import re
import os
import sys
import pysam
import pandas as pd

def extract_sample_and_caller(filename):
    base = os.path.basename(filename)
    parts = base.split('.')
    sample = parts[0]
    caller = parts[1]
    return sample, caller

def parse_format_and_samples(format_field, sample_field):
    return dict(zip(format_field.split(':'), sample_field.split(':')))

def extract_delly(record, format_field, sample_field):
    data = parse_format_and_samples(format_field, sample_field)
    return {
        'PE': record.info.get('PE', 'NA'),
        'SR': record.info.get('SR', 'NA'),
        'CT': record.info.get('CT', 'NA'),
        'DP': data.get('RC', 'NA')
    }

def extract_gridss(record, format_field, sample_field):
    return {
        'PE': record.info.get('ASRP', 'NA'),
        'SR': record.info.get('SR', 'NA'),
        'CT': record.info.get('INSRMRC', 'NA'),
        'DP': record.info.get('ASSR', 'NA')
    }

def extract_manta(record, format_field, sample_field):
    alt = record.alts[0] if record.alts else ''
    connection_type = '-' if '[' in alt or ']' in alt else '+'
    data = parse_format_and_samples(format_field, sample_field)
    pe = data.get('PR', 'NA').split(',')[1] if 'PR' in data else 'NA'
    return {
        'PE': pe,
        'SR': 'NA',
        'CT': connection_type,
        'DP': record.info.get('BND_DEPTH', 'NA')
    }

def extract_svaba(record, format_field, sample_field):
    sctg = record.info.get('SCTG', 'NA')
    ct = 'NA'
    if sctg != 'NA' and '__' in sctg:
        ct_parts = sctg.split('__')[0].split('-')
        ct = f"{ct_parts[0][-1]}/{ct_parts[1][-1]}" if len(ct_parts) > 1 else 'NA'
    data = parse_format_and_samples(format_field, sample_field)
    return {
        'PE': data.get('DR', 'NA'),
        'SR': data.get('SR', 'NA'),
        'CT': ct,
        'DP': data.get('DP', 'NA')
    }

def extract_tiddit(record, format_field, sample_field):
    alt = record.alts[0] if record.alts else ''
    connection_type = '-' if '[' in alt or ']' in alt else '+'
    data = parse_format_and_samples(format_field, sample_field)
    return {
        'PE': data.get('DV', 'NA'),
        'SR': data.get('RV', 'NA'),
        'CT': connection_type,
        'DP': record.info.get('LTE', 'NA')
    }

CALLER_EXTRACTORS = {
    'delly': extract_delly,
    'gridss': extract_gridss,
    'manta': extract_manta,
    'svaba': extract_svaba,
    'tiddit': extract_tiddit,
}


def extract_chr2_from_alt(alt, chrom):
    if not alt or ('[' not in alt and ']' not in alt):
        return chrom 
    match = re.search(r'[\[\]]([A-Za-z]*\d+|[XYMT])[:]', alt)
    if match:
        return match.group(1)
    return chrom 

def process_record(record, caller):
    chrom = record.chrom
    pos = record.pos
    ref = record.ref
    alt = ','.join(record.alts) if record.alts else '.'
    qual = record.qual if record.qual is not None else '.'
    filt = ';'.join(record.filter.keys()) if record.filter else '.'

    info = record.info
    svtype = info.get('SVTYPE', 'NA')
    end = info.get('END', pos + 1)
    svlen = info.get('SVLEN', [end - pos])
    if isinstance(svlen, (list, tuple)):
        svlen = svlen[0] if svlen else end - pos

    ciend = ','.join(map(str, info['CIEND'])) if 'CIEND' in info else 'NA'
    cipos = ','.join(map(str, info['CIPOS'])) if 'CIPOS' in info else 'NA'
    supp_vec = info['SUPP_VEC'] if 'SUPP_VEC' in info else 'NA'
    supp = info['SUPP'] if 'SUPP' in info else 'NA'
    svmethod = info['SVMETHOD'] if 'SVMETHOD' in info else 'NA'
    strands = info['STRANDS'] if 'STRANDS' in info else 'NA'

    chr2 = extract_chr2_from_alt(alt, chrom)

    sample = list(record.samples)[0]
    sample_data = record.samples[sample]
    format_field = ':'.join(sample_data.keys())
    sample_field = ':'.join([
        ','.join(map(str, sample_data[key])) if isinstance(sample_data[key], (list, tuple)) 
        else str(sample_data[key])
        for key in sample_data
    ])

    if caller in CALLER_EXTRACTORS:
        support_data = CALLER_EXTRACTORS[caller](record, format_field, sample_field)
    else:
        support_data = {'PE': 'NA', 'SR': 'NA', 'CT': 'NA', 'DP': 'NA'}

    return {
        'CHROM': chrom,
        'POS': pos,
        'ID': record.id if record.id else '.',
        'REF': ref,
        'ALT': alt,
        'QUAL': qual,
        'FILTER': filt,
        'CIEND': ','.join(map(str, ciend)),
        'CIPOS': ','.join(map(str, cipos)),
        'CHR2': chr2,
        'END': end,
        'SVLEN': svlen,
        'SVMETHOD': svmethod,
        'SVTYPE': svtype,
        'SUPP_VEC': supp_vec,
        'SUPP': supp,
        'STRANDS': strands,
        **support_data
    }, record

def parse_and_standardize_vcf(input_vcf):
    sample, caller = extract_sample_and_caller(input_vcf)

    standardized_vcf = f"{sample}.{caller}.standard.vcf"
    summary_tsv = f"{sample}.{caller}.summary.tsv"

    variants = []

    with pysam.VariantFile(input_vcf) as vcf_in, \
         pysam.VariantFile(standardized_vcf, 'w', header=vcf_in.header) as vcf_out:

        for record in vcf_in:
            variant_data, standardized_record = process_record(record, caller)
            variants.append(variant_data)
            vcf_out.write(standardized_record)

    df = pd.DataFrame(variants)
    df.to_csv(summary_tsv, sep='\t', index=False)

    print(f"Processed: {input_vcf} → {standardized_vcf}, {summary_tsv}")

def main():
    if len(sys.argv) < 2:
        print("Usage: python standardize_svs.py file1.vcf [file2.vcf ...]")
        sys.exit(1)

    for vcf_file in sys.argv[1:]:
        parse_and_standardize_vcf(vcf_file)

if __name__ == '__main__':
    main()

