#!/usr/bin/env python3

import re
import sys
import gzip
import argparse
import pandas as pd
from collections import defaultdict

VERSION = "2.2.1"

CORE_COLUMNS = {
    "KEY", "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "AF", "VD", "DP",
    "Caller", "CALLER_INFO", "CALLER_FORMAT", "CALLER_SAMPLE_NAMES", "CALLER_SAMPLES",
    "SELECTED_SAMPLE", "SNPEFF_ANN", "VEP_CSQ",
}

ANNOTATION_PRIORITY_FIELDS = [
    "VEP_CONSEQUENCE", "VEP_SYMBOL", "VEP_HGVSC", "VEP_HGVSP",
    "SNPEFF_ANNOTATION", "SNPEFF_GENE_NAME", "SNPEFF_HGVS_C", "SNPEFF_HGVS_P",
]

SNPEFF_ANN_FIELDS = [
    "Allele", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID",
    "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank", "HGVS.c", "HGVS.p",
    "cDNA.pos / cDNA.length", "CDS.pos / CDS.length", "AA.pos / AA.length",
    "Distance", "ERRORS / WARNINGS / INFO",
]

def get_version():
    return VERSION

def open_vcf(file_path):
    if file_path.endswith(".gz"):
        return gzip.open(file_path, "rt")
    return open(file_path)

def parse_info_field(info):
    info_dict = {}
    if not info or info == ".":
        return info_dict
    for item in info.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            info_dict[key] = value
        else:
            info_dict[item] = "true"
    return info_dict

def strip_annotation_tags(info):
    """Drop ANN/CSQ from the raw INFO blob; they are stored as dedicated columns."""
    if not info or info == ".":
        return info
    kept = []
    for item in info.split(";"):
        if not item:
            continue
        key = item.split("=", 1)[0]
        if key not in {"ANN", "CSQ"}:
            kept.append(item)
    return ";".join(kept) if kept else "."

def sanitize_column_name(value):
    value = re.sub(r"[^0-9A-Za-z]+", "_", value.strip())
    value = re.sub(r"_+", "_", value).strip("_")
    return value.upper() if value else "FIELD"

def parse_vcf_header(file_path):
    samples = []
    csq_fields = []
    with open_vcf(file_path) as f:
        for line in f:
            if line.startswith("##INFO=<ID=CSQ"):
                match = re.search(r"Format: ([^\"]+)", line)
                if match:
                    csq_fields = [x.strip() for x in match.group(1).split("|")]
            elif line.startswith("#CHROM"):
                fields = line.rstrip("\n").split("\t")
                samples = fields[9:]
                break
    return samples, csq_fields

def split_alt_values(value, alt_idx):
    if value is None:
        return None
    values = value.split(",")
    if len(values) == 1:
        return values[0]
    if alt_idx < len(values):
        return values[alt_idx]
    return values[0]

def matching_annotation_blob(annotation_blob, alt):
    if not annotation_blob:
        return None
    entries = [x for x in annotation_blob.split(",") if x]
    if not entries:
        return None
    matches = [entry for entry in entries if entry.split("|", 1)[0] == alt]
    return ",".join(matches) if matches else entries[0]

def parse_pipe_annotation(annotation_blob, alt, fields, prefix):
    chosen_blob = matching_annotation_blob(annotation_blob, alt)
    if not chosen_blob:
        return {}
    chosen = None
    entries = [x for x in chosen_blob.split(",") if x]
    for entry in entries:
        parts = entry.split("|")
        if parts and parts[0] == alt:
            chosen = parts
            break
    if chosen is None:
        chosen = entries[0].split("|") if entries else []
    parsed = {}
    for i, field_name in enumerate(fields):
        value = chosen[i] if i < len(chosen) and chosen[i] != "" else None
        parsed[f"{prefix}_{sanitize_column_name(field_name)}"] = value
    return parsed

def parse_annotations(info_dict, alt, csq_fields):
    ann = matching_annotation_blob(info_dict.get("ANN"), alt)
    csq = matching_annotation_blob(info_dict.get("CSQ"), alt)
    parsed = {
        "SNPEFF_ANN": ann,
        "VEP_CSQ": csq,
    }
    parsed.update(parse_pipe_annotation(info_dict.get("ANN"), alt, SNPEFF_ANN_FIELDS, "SNPEFF"))
    parsed.update(parse_pipe_annotation(info_dict.get("CSQ"), alt, csq_fields, "VEP"))
    return parsed

def add_info_columns(record, info_dict, alt_idx):
    for key, value in info_dict.items():
        if key in {"ANN", "CSQ"}:
            continue
        record[f"INFO_{sanitize_column_name(key)}"] = split_alt_values(value, alt_idx)

def add_format_columns(record, fmt_dict):
    for key, value in fmt_dict.items():
        record[f"FORMAT_{sanitize_column_name(key)}"] = value

def format_sample_payload(sample_names, sample_values):
    if not sample_values:
        return None
    if sample_names and len(sample_names) == len(sample_values):
        return ";".join(f"{name}={value}" for name, value in zip(sample_names, sample_values))
    return ";".join(sample_values)

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
    sample_names, csq_fields = parse_vcf_header(file_path)
    with open_vcf(file_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 10:
                continue
            chrom, pos, vid, ref, alt, qual, fltr, info, fmt = fields[:9]
            sample_values = fields[9:]
            sample = sample_values[0]
            info_dict = parse_info_field(info)
            fmt_dict = dict(zip(fmt.split(':'), sample.split(':')))
            alts = alt.split(",")
            for alt_idx, alt_allele in enumerate(alts):
                af = split_alt_values(fmt_dict.get("AF") or info_dict.get("AF"), alt_idx)
                vd = split_alt_values(fmt_dict.get("AO") or info_dict.get("AO"), alt_idx)
                dp = fmt_dict.get("DP") or info_dict.get("DP")
                if (not af or af == ".") and vd and dp:
                    try:
                        af = str(round(float(vd) / float(dp), 4))
                    except ZeroDivisionError:
                        af = "0"
                record = {
                    "KEY": f"{chrom}_{pos}_{ref}_{alt_allele}",
                    "CHROM": chrom,
                    "POS": int(pos),
                    "ID": vid,
                    "REF": ref,
                    "ALT": alt_allele,
                    "QUAL": qual,
                    "FILTER": fltr,
                    "AF": af,
                    "VD": vd,
                    "DP": dp,
                    "Caller": "freebayes",
                    "CALLER_INFO": strip_annotation_tags(info),
                    "CALLER_FORMAT": fmt,
                    "CALLER_SAMPLE_NAMES": ";".join(sample_names) if sample_names else None,
                    "CALLER_SAMPLES": format_sample_payload(sample_names, sample_values),
                    "SELECTED_SAMPLE": sample_names[0] if sample_names else None,
                }
                record.update(parse_annotations(info_dict, alt_allele, csq_fields))
                add_info_columns(record, info_dict, alt_idx)
                add_format_columns(record, fmt_dict)
                records.append(record)
    return pd.DataFrame(records)

def parse_lofreq(file_path, caller_label):
    records = []
    sample_names, csq_fields = parse_vcf_header(file_path)
    with open_vcf(file_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue
            chrom, pos, vid, ref, alt, qual, fltr, info = fields[:8]
            info_dict = parse_info_field(info)
            dp = info_dict.get("DP")
            dp4 = info_dict.get("DP4")
            alts = alt.split(",")
            for alt_idx, alt_allele in enumerate(alts):
                af = split_alt_values(info_dict.get("AF"), alt_idx)
                vd = None
                if dp4 and len(alts) == 1:
                    parts = list(map(int, dp4.split(',')))
                    vd = parts[2] + parts[3]
                if (not af or af == ".") and vd and dp:
                    try:
                        af = str(round(float(vd) / float(dp), 4))
                    except ZeroDivisionError:
                        af = "0"
                record = {
                    "KEY": f"{chrom}_{pos}_{ref}_{alt_allele}",
                    "CHROM": chrom,
                    "POS": int(pos),
                    "ID": vid,
                    "REF": ref,
                    "ALT": alt_allele,
                    "QUAL": qual,
                    "FILTER": fltr,
                    "AF": af,
                    "VD": str(vd) if vd is not None else None,
                    "DP": dp,
                    "Caller": caller_label,
                    "CALLER_INFO": strip_annotation_tags(info),
                    "CALLER_FORMAT": None,
                    "CALLER_SAMPLE_NAMES": ";".join(sample_names) if sample_names else None,
                    "CALLER_SAMPLES": None,
                    "SELECTED_SAMPLE": None,
                }
                record.update(parse_annotations(info_dict, alt_allele, csq_fields))
                add_info_columns(record, info_dict, alt_idx)
                records.append(record)
    return pd.DataFrame(records)

def parse_mutect2(file_path):
    records = []
    sample_names, csq_fields = parse_vcf_header(file_path)
    with open_vcf(file_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 10:
                continue
            chrom, pos, vid, ref, alt, qual, fltr, info, fmt = fields[:9]
            sample_fields = fields[9:]
            sample = sample_fields[-1]
            selected_sample = sample_names[-1] if sample_names else None
            fmt_dict = dict(zip(fmt.split(":"), sample.split(":")))
            info_dict = parse_info_field(info)
            alts = alt.split(",")
            for alt_idx, alt_allele in enumerate(alts):
                af = split_alt_values(fmt_dict.get("AF") or info_dict.get("AF"), alt_idx)
                dp = fmt_dict.get("DP") or info_dict.get("DP")
                vd = None
                ad = fmt_dict.get("AD")
                if ad:
                    ad_parts = ad.split(",")
                    if alt_idx + 1 < len(ad_parts):
                        vd = ad_parts[alt_idx + 1]
                if not vd:
                    vd = split_alt_values(info_dict.get("AS_UNIQ_ALT_READ_COUNT"), alt_idx)
                if not vd and "AS_SB_TABLE" in info_dict:
                    try:
                        per_allele = info_dict["AS_SB_TABLE"].split("|")
                        if alt_idx + 1 < len(per_allele):
                            alt_forward, alt_reverse = map(int, per_allele[alt_idx + 1].split(","))
                            vd = str(alt_forward + alt_reverse)
                    except Exception:
                        vd = None
                if (not af or af == ".") and vd and dp:
                    try:
                        af = str(round(float(vd) / float(dp), 4))
                    except ZeroDivisionError:
                        af = "0"
                record = {
                    "KEY": f"{chrom}_{pos}_{ref}_{alt_allele}",
                    "CHROM": chrom,
                    "POS": int(pos),
                    "ID": vid,
                    "REF": ref,
                    "ALT": alt_allele,
                    "QUAL": qual,
                    "FILTER": fltr,
                    "AF": af,
                    "VD": vd,
                    "DP": dp,
                    "Caller": "mutect2",
                    "CALLER_INFO": strip_annotation_tags(info),
                    "CALLER_FORMAT": fmt,
                    "CALLER_SAMPLE_NAMES": ";".join(sample_names) if sample_names else None,
                    "CALLER_SAMPLES": format_sample_payload(sample_names, sample_fields),
                    "SELECTED_SAMPLE": selected_sample,
                }
                record.update(parse_annotations(info_dict, alt_allele, csq_fields))
                add_info_columns(record, info_dict, alt_idx)
                add_format_columns(record, fmt_dict)
                records.append(record)
    return pd.DataFrame(records)

def parse_vardict(file_path):
    records = []
    sample_names, csq_fields = parse_vcf_header(file_path)
    with open_vcf(file_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 10:
                continue
            chrom, pos, vid, ref, alt, qual, fltr, info, fmt = fields[:9]
            sample_values = fields[9:]
            sample = sample_values[0]
            info_dict = parse_info_field(info)
            fmt_dict = dict(zip(fmt.split(':'), sample.split(':')))
            alts = alt.split(",")
            for alt_idx, alt_allele in enumerate(alts):
                af = split_alt_values(fmt_dict.get("AF") or info_dict.get("AF"), alt_idx)
                vd = split_alt_values(fmt_dict.get("VD") or info_dict.get("VD"), alt_idx)
                dp = fmt_dict.get("DP") or info_dict.get("DP")
                if (not af or af == ".") and vd and dp:
                    try:
                        af = str(round(float(vd) / float(dp), 4))
                    except ZeroDivisionError:
                        af = "0"
                record = {
                    "KEY": f"{chrom}_{pos}_{ref}_{alt_allele}",
                    "CHROM": chrom,
                    "POS": int(pos),
                    "ID": vid,
                    "REF": ref,
                    "ALT": alt_allele,
                    "QUAL": qual,
                    "FILTER": fltr,
                    "AF": af,
                    "VD": vd,
                    "DP": dp,
                    "Caller": "vardict",
                    "CALLER_INFO": strip_annotation_tags(info),
                    "CALLER_FORMAT": fmt,
                    "CALLER_SAMPLE_NAMES": ";".join(sample_names) if sample_names else None,
                    "CALLER_SAMPLES": format_sample_payload(sample_names, sample_values),
                    "SELECTED_SAMPLE": sample_names[0] if sample_names else None,
                }
                record.update(parse_annotations(info_dict, alt_allele, csq_fields))
                add_info_columns(record, info_dict, alt_idx)
                add_format_columns(record, fmt_dict)
                records.append(record)
    return pd.DataFrame(records)

def merged_af(records_by_caller):
    """Merged AF: VarDict when present, otherwise mean of Mutect2/LoFreq (never
    FreeBayes)."""
    vardict = records_by_caller.get("vardict")
    if vardict is not None:
        af = safe_float(vardict.get("AF"))
        if af is not None:
            return round(af, 4)
    af_values = []
    for caller in ("mutect2", "lofreq_snv", "lofreq_indels"):
        row = records_by_caller.get(caller)
        if row is not None:
            af = safe_float(row.get("AF"))
            if af is not None:
                af_values.append(af)
    if af_values:
        return round(sum(af_values) / len(af_values), 4)
    return None

def merge_callers(dfs, caller_order):
    grouped = defaultdict(list)
    for df in dfs:
        for row in df.to_dict("records"):
            grouped[row['KEY']].append(row)
    merged = []
    for key, records in grouped.items():
        base = records[0]
        records_by_caller = {r['Caller']: r for r in records}
        supp_vec = ''.join(['1' if any(r['Caller'] == c for r in records) else '0' for c in caller_order])
        vd_values = [int(r['VD']) for r in records if r['VD'] and r['VD'].isdigit()]
        dp_values = [int(r['DP']) for r in records if r['DP'] and r['DP'].isdigit()]
        ann_source = None
        for caller in caller_order:
            row = records_by_caller.get(caller)
            if row is not None and any(row.get(k) for k in ANNOTATION_PRIORITY_FIELDS):
                ann_source = row
                break
        if ann_source is None:
            ann_source = base
        merged_row = {
            'CHROM': base['CHROM'],
            'POS': base['POS'],
            'ID': base['ID'],
            'REF': base['REF'],
            'ALT': base['ALT'],
            'FILTER': 'PASS' if any(r['FILTER'] == 'PASS' for r in records) else records[0]['FILTER'],
            'AF': merged_af(records_by_caller),
            'VD': int(round(sum(vd_values) / len(vd_values))) if vd_values else None,
            'DP': int(round(sum(dp_values) / len(dp_values))) if dp_values else None,
            'NUM_CALLERS': len(records_by_caller),
            'CALLERS': ",".join([c for c in caller_order if c in records_by_caller]),
            'SUPP_VEC': supp_vec,
        }
        for col in ann_source.keys():
            if col.startswith("VEP_") or col.startswith("SNPEFF_"):
                merged_row[col] = ann_source.get(col)
        for caller in caller_order:
            row = records_by_caller.get(caller)
            prefix = sanitize_column_name(caller)
            merged_row[f"{prefix}_PRESENT"] = 1 if row is not None else 0
            if row is None:
                continue
            for col in row.keys():
                if col in {"KEY", "CHROM", "POS", "REF", "ALT", "Caller"}:
                    continue
                merged_row[f"{prefix}_{col}"] = row.get(col)
        merged.append(merged_row)
    return pd.DataFrame(merged)

def main():
    parser = argparse.ArgumentParser(
        description="Merge annotated VCFs from FreeBayes, LoFreq_SNV, LoFreq_Indels, Mutect2, and VarDict into a unified TSV with SUPP_VEC and annotation fields.",
        epilog="Expected order of files: freebayes.vcf(.gz) lofreq_snvs.vcf(.gz) lofreq_indels.vcf(.gz) mutect2.vcf(.gz) vardict.vcf(.gz)"
    )
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {get_version()}")
    parser.add_argument("vcfs", nargs=5, metavar="VCF", help="5 VCF files from the 4 callers (in specific order)")
    parser.add_argument("-o", "--output", default="merged_variants.tsv", help="Output TSV file")
    args = parser.parse_args()
    parsers = [
        parse_freebayes,
        lambda p: parse_lofreq(p, "lofreq_snv"),
        lambda p: parse_lofreq(p, "lofreq_indels"),
        parse_mutect2,
        parse_vardict
    ]
    caller_order = ['freebayes', 'lofreq_snv', 'lofreq_indels', 'mutect2', 'vardict']
    dfs = []
    for func, file_path in zip(parsers, args.vcfs):
        try:
            df = func(file_path)
            dfs.append(df)
        except Exception as e:
            print(f"Error parsing {file_path} with {func.__name__}: {e}", file=sys.stderr)
            sys.exit(1)
    merged_df = merge_callers(dfs, caller_order).fillna(".").replace("", ".")
    merged_df.to_csv(args.output, sep='\t', index=False, na_rep=".")
    print(f"Merged TSV written to {args.output}")

if __name__ == "__main__":
    main()
