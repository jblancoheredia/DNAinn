#!/usr/bin/env python3

import argparse
import pandas as pd
import sys
from collections import defaultdict

def safe_float(value):
    if isinstance(value, str) and ',' in value:
        try:
            parts = [float(x) for x in value.split(',') if x.strip() != '']
            return sum(parts) / len(parts) if parts else None
        except ValueError:
            return None
    try:
        return float(value)
    except (ValueError, TypeError):
        return None

def parse_freebayes(file_path):
    records = []
    with open(file_path) as f:
        for line in f:
            if line.startswith("#"): continue
            fields = line.strip().split('\t')
            chrom, pos, vid, ref, alt, qual, fltr, info, fmt, sample = fields[:10]
            info_dict = dict(kv.split('=') for kv in info.split(';') if '=' in kv)
            fmt_dict = dict(zip(fmt.split(':'), sample.split(':')))
            af = fmt_dict.get("AF") or info_dict.get("AF")
            vd = fmt_dict.get("AO")
            dp = fmt_dict.get("DP")
            if (not af or af == ".") and vd and dp:
                try: af = str(round(float(vd) / float(dp), 4))
                except ZeroDivisionError: af = "0"
            records.append({"KEY": f"{chrom}_{pos}_{ref}_{alt}", "CHROM": chrom, "POS": int(pos), "ID": vid,
                            "REF": ref, "ALT": alt, "QUAL": qual, "FILTER": fltr, "AF": af,
                            "VD": vd, "DP": dp, "Caller": "freebayes"})
    return pd.DataFrame(records)

def parse_lofreq(file_path):
    records = []
    with open(file_path) as f:
        for line in f:
            if line.startswith("#"): continue
            fields = line.strip().split('\t')
            chrom, pos, vid, ref, alt, qual, fltr, info = fields[:8]
            info_dict = dict(kv.split('=') for kv in info.split(';') if '=' in kv)
            dp = info_dict.get("DP")
            dp4 = info_dict.get("DP4")
            af = info_dict.get("AF")
            vd = None
            if dp4:
                parts = list(map(int, dp4.split(',')))
                vd = parts[2] + parts[3]
            if (not af or af == ".") and vd and dp:
                try: af = str(round(float(vd) / float(dp), 4))
                except ZeroDivisionError: af = "0"
            records.append({"KEY": f"{chrom}_{pos}_{ref}_{alt}", "CHROM": chrom, "POS": int(pos), "ID": vid,
                            "REF": ref, "ALT": alt, "QUAL": qual, "FILTER": fltr, "AF": af,
                            "VD": str(vd) if vd is not None else None, "DP": dp, "Caller": "lofreq"})
    return pd.DataFrame(records)

def parse_mutect2(file_path):
    records = []
    with open(file_path) as f:
        for line in f:
            if line.startswith("#"): continue
            fields = line.strip().split('\t')
            chrom, pos, vid, ref, alt, qual, fltr, info, *ids, match = fields
            fmt, sample = ids if len(ids) >= 2 else ("", "")
            info_dict = dict(kv.split('=') for kv in info.split(';') if '=' in kv)
            af = info_dict.get("AF")
            dp = info_dict.get("DP")
            vd = info_dict.get("AS_UNIQ_ALT_READ_COUNT")
            if not vd and "AS_SB_TABLE" in info_dict:
                try:
                    alt_counts = info_dict["AS_SB_TABLE"].split('|')[1]
                    alt_forward, alt_reverse = map(int, alt_counts.split(','))
                    vd = str(alt_forward + alt_reverse)
                except Exception:
                    vd = None
            if (not af or af == ".") and vd and dp:
                try: af = str(round(float(vd) / float(dp), 4))
                except ZeroDivisionError: af = "0"
            records.append({"KEY": f"{chrom}_{pos}_{ref}_{alt}", "CHROM": chrom, "POS": int(pos), "ID": vid,
                            "REF": ref, "ALT": alt, "QUAL": qual, "FILTER": fltr, "AF": af,
                            "VD": vd, "DP": dp, "Caller": "mutect2"})
    return pd.DataFrame(records)

def parse_vardict(file_path):
    records = []
    with open(file_path) as f:
        for line in f:
            if line.startswith("#"): continue
            fields = line.strip().split('\t')
            chrom, pos, vid, ref, alt, qual, fltr, info, fmt, sample, match = fields[:11]
            fmt_dict = dict(zip(fmt.split(':'), sample.split(':')))
            af = fmt_dict.get("AF")
            vd = fmt_dict.get("VD")
            dp = fmt_dict.get("DP")
            if (not af or af == ".") and vd and dp:
                try: af = str(round(float(vd) / float(dp), 4))
                except ZeroDivisionError: af = "0"
            records.append({"KEY": f"{chrom}_{pos}_{ref}_{alt}", "CHROM": chrom, "POS": int(pos), "ID": vid,
                            "REF": ref, "ALT": alt, "QUAL": qual, "FILTER": fltr, "AF": af,
                            "VD": vd, "DP": dp, "Caller": "vardir"})
    return pd.DataFrame(records)

def merge_callers(dfs, caller_order):
    grouped = defaultdict(list)
    for df in dfs:
        for _, row in df.iterrows():
            grouped[row['KEY']].append(row)

    merged = []
    for key, records in grouped.items():
        base = records[0]
        supp_vec = ''.join(['1' if any(r['Caller'] == c for r in records) else '0' for c in caller_order])
        af_values = [safe_float(r['AF']) for r in records if safe_float(r['AF']) is not None]
        vd_values = [int(r['VD']) for r in records if r['VD'] and r['VD'].isdigit()]
        dp_values = [int(r['DP']) for r in records if r['DP'] and r['DP'].isdigit()]
        merged_row = {
            'CHROM': base['CHROM'],
            'POS': base['POS'],
            'ID': base['ID'],
            'REF': base['REF'],
            'ALT': base['ALT'],
            'FILTER': 'PASS' if any(r['FILTER'] == 'PASS' for r in records) else records[0]['FILTER'],
            'AF': round(sum(af_values) / len(af_values), 4) if af_values else None,
            'VD': int(round(sum(vd_values) / len(vd_values))) if vd_values else None,
            'DP': int(round(sum(dp_values) / len(dp_values))) if dp_values else None,
            'SUPP_VEC': supp_vec
        }
        merged.append(merged_row)
    return pd.DataFrame(merged)

def main():
    parser = argparse.ArgumentParser(
        description="Merge VCFs from FreeBayes, LoFreq_SNV, LoFreq_Indels, Mutect2, and VarDict into a unified TSV with SUPP_VEC.",
        epilog="Expected order of files: freebayes.vcf lofreq_snvs.vcf lofreq_indels.vcf mutect2.vcf vardir.vcf"
    )
    parser.add_argument("vcfs", nargs=5, metavar="VCF", help="5 VCF files from the 4 callers (in specific order)")
    parser.add_argument("-o", "--output", default="merged_variants.tsv", help="Output TSV file")
    args = parser.parse_args()

    parsers = [parse_freebayes, parse_lofreq, parse_lofreq, parse_mutect2, parse_vardict]
    caller_order = ['freebayes', 'lofreq_snv', 'lofreq_indels', 'mutect2', 'vardir']

    dfs = []
    for func, file_path in zip(parsers, args.vcfs):
        try:
            df = func(file_path)
            dfs.append(df)
        except Exception as e:
            print(f"Error parsing {file_path} with {func.__name__}: {e}", file=sys.stderr)
            sys.exit(1)

    merged_df = merge_callers(dfs, caller_order)
    merged_df.to_csv(args.output, sep='\t', index=False)
    print(f"Merged TSV written to {args.output}")

if __name__ == "__main__":
    main()
