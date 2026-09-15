#!/usr/bin/env python3

import argparse
import csv
import re
import sys

VERSION = "1.3.1"

def _raise_csv_field_size_limit():
    """Merged VCFs can carry VEP/SnpEff blobs larger than csv's 128 KiB default."""
    limit = sys.maxsize
    while True:
        try:
            csv.field_size_limit(limit)
            return
        except OverflowError:
            limit = limit // 10

_raise_csv_field_size_limit()

OUTPUT_COLUMNS = [
    "Sample", "NormalUsed", "Chrom", "Start", "Ref", "Alt", "VariantClass",
    "Gene", "Exon", "Call_Confidence", "Comments", "TranscriptID",
    "cDNAchange", "AAchange", "dbSNP_ID", "Cosmic_ID", "1000G_MAF",
    "FailureReason", "CallMethod", "COSMIC_site", "N_TotalDepth",
    "N_RefCount", "N_AltCount", "N_AltFreq", "T_TotalDepth", "T_RefCount",
    "T_AltCount", "T_AltFreq", "T_Ref+", "T_Ref-", "T_Alt+", "T_Alt-",
    "Strand_Bias", "All_N_Aggregate_AlleleDepth", "All_N_Median_AlleleFreq",
    "T_freq/All_N_Freq", "Occurence_in_Normals",
    "Germline", "Germline_Info", "Hotspot", "Hotspot_Info",
]

CALLER_ORDER = ["VARDICT", "MUTECT2", "FREEBAYES", "LOFREQ_SNV", "LOFREQ_INDELS"]

AA3 = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
    "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
    "*": "*",
}


def clean(value):
    if value in (None, "", "."):
        return "."
    return str(value)


def first_value(row, fields):
    for field in fields:
        value = row.get(field)
        if value not in (None, "", "."):
            return value
    return "."


def as_int(value):
    try:
        if value in (None, "", "."):
            return None
        return int(round(float(str(value).split(",")[0])))
    except (TypeError, ValueError):
        return None


def as_float(value):
    try:
        if value in (None, "", "."):
            return None
        return float(str(value).split(",")[0])
    except (TypeError, ValueError):
        return None


def split_ints(value):
    if value in (None, "", "."):
        return []
    out = []
    for part in str(value).split(","):
        try:
            out.append(int(round(float(part))))
        except ValueError:
            out.append(None)
    return out


def fmt_number(value):
    if value is None:
        return "."
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def sample_names(row, caller):
    names = row.get(f"{caller}_CALLER_SAMPLE_NAMES")
    if names in (None, "", "."):
        return []
    return names.split(";")


def sample_maps(row, caller):
    payload = row.get(f"{caller}_CALLER_SAMPLES")
    fmt = row.get(f"{caller}_CALLER_FORMAT")
    if payload in (None, "", ".") or fmt in (None, "", "."):
        return []
    keys = fmt.split(":")
    maps = []
    for sample_payload in payload.split(";"):
        if "=" in sample_payload:
            name, values = sample_payload.split("=", 1)
        else:
            name, values = ".", sample_payload
        maps.append((name, dict(zip(keys, values.split(":")))))
    return maps


def preferred_sample_names(row):
    for caller in ["VARDICT", "MUTECT2", "FREEBAYES", "LOFREQ_SNV", "LOFREQ_INDELS"]:
        names = sample_names(row, caller)
        if names:
            if len(names) > 1:
                return names[-2], names[-1]
            return names[0], "."
    return ".", "."


def vardict_counts(row):
    maps = sample_maps(row, "VARDICT")
    if not maps:
        return None

    def counts_from_map(sample):
        dp = as_int(sample.get("DP"))
        vd = as_int(sample.get("VD"))
        af = as_float(sample.get("AF"))
        rd = split_ints(sample.get("RD"))
        ald = split_ints(sample.get("ALD"))
        ref_plus = rd[0] if len(rd) > 0 else None
        ref_minus = rd[1] if len(rd) > 1 else None
        alt_plus = ald[0] if len(ald) > 0 else None
        alt_minus = ald[1] if len(ald) > 1 else None
        ref = None
        if ref_plus is not None and ref_minus is not None:
            ref = ref_plus + ref_minus
        elif dp is not None and vd is not None:
            ref = dp - vd
        return {
            "dp": dp, "ref": ref, "alt": vd, "af": af,
            "ref_plus": ref_plus, "ref_minus": ref_minus,
            "alt_plus": alt_plus, "alt_minus": alt_minus,
        }

    tumor = counts_from_map(maps[0][1])
    normal = counts_from_map(maps[1][1]) if len(maps) > 1 else {}
    return tumor, normal


def mutect_counts(row):
    maps = sample_maps(row, "MUTECT2")
    if not maps:
        return None

    def counts_from_map(sample):
        dp = as_int(sample.get("DP"))
        ad = split_ints(sample.get("AD"))
        af = as_float(sample.get("AF"))
        sb = split_ints(sample.get("SB"))
        ref = ad[0] if len(ad) > 0 else None
        alt = sum(x for x in ad[1:] if x is not None) if len(ad) > 1 else None
        ref_plus = sb[0] if len(sb) > 0 else None
        ref_minus = sb[1] if len(sb) > 1 else None
        alt_plus = sb[2] if len(sb) > 2 else None
        alt_minus = sb[3] if len(sb) > 3 else None
        return {
            "dp": dp, "ref": ref, "alt": alt, "af": af,
            "ref_plus": ref_plus, "ref_minus": ref_minus,
            "alt_plus": alt_plus, "alt_minus": alt_minus,
        }

    normal = counts_from_map(maps[0][1])
    tumor = counts_from_map(maps[1][1]) if len(maps) > 1 else {}
    return tumor, normal


def freebayes_counts(row):
    dp = as_int(row.get("FREEBAYES_DP"))
    alt = as_int(row.get("FREEBAYES_VD"))
    af = as_float(row.get("FREEBAYES_AF"))
    ref = dp - alt if dp is not None and alt is not None else None
    return {
        "dp": dp,
        "ref": ref,
        "alt": alt,
        "af": af,
        "ref_plus": as_int(row.get("FREEBAYES_INFO_SRF")),
        "ref_minus": as_int(row.get("FREEBAYES_INFO_SRR")),
        "alt_plus": as_int(row.get("FREEBAYES_INFO_SAF")),
        "alt_minus": as_int(row.get("FREEBAYES_INFO_SAR")),
    }, {}


def lofreq_counts(row):
    caller = "LOFREQ_SNV" if row.get("LOFREQ_SNV_PRESENT") == "1" else "LOFREQ_INDELS"
    dp4 = split_ints(row.get(f"{caller}_INFO_DP4"))
    dp = as_int(row.get(f"{caller}_DP"))
    alt = as_int(row.get(f"{caller}_VD"))
    af = as_float(row.get(f"{caller}_AF"))
    ref_plus = dp4[0] if len(dp4) > 0 else None
    ref_minus = dp4[1] if len(dp4) > 1 else None
    alt_plus = dp4[2] if len(dp4) > 2 else None
    alt_minus = dp4[3] if len(dp4) > 3 else None
    ref = None
    if ref_plus is not None and ref_minus is not None:
        ref = ref_plus + ref_minus
    elif dp is not None and alt is not None:
        ref = dp - alt
    return {
        "dp": dp, "ref": ref, "alt": alt, "af": af,
        "ref_plus": ref_plus, "ref_minus": ref_minus,
        "alt_plus": alt_plus, "alt_minus": alt_minus,
    }, {}


def best_counts(row):
    for caller, fn in [
        ("VARDICT", vardict_counts),
        ("MUTECT2", mutect_counts),
        ("FREEBAYES", freebayes_counts),
        ("LOFREQ_SNV", lofreq_counts),
        ("LOFREQ_INDELS", lofreq_counts),
    ]:
        if row.get(f"{caller}_PRESENT") == "1":
            counts = fn(row)
            if counts:
                return counts
    return {}, {}


def extract_existing_ids(*values):
    """Categorise variant IDs from one or more sources (VEP Existing_variation
    and/or the hotspot-list match) into dbSNP (rs*) and COSMIC (COSM*/COSV*).
    De-duplicated, first-seen order preserved. Other IDs (e.g. HGMD CM*) are
    ignored since the clinical MAF has no column for them."""
    dbsnp = []
    cosmic = []
    for value in values:
        if value in (None, "", "."):
            continue
        for item in re.split(r"[,&;|]", value):
            item = item.strip()
            if item.startswith("rs"):
                if item not in dbsnp:
                    dbsnp.append(item)
            elif item.upper().startswith(("COSM", "COSV")):
                if item not in cosmic:
                    cosmic.append(item)
    return ",".join(dbsnp) if dbsnp else ".", ",".join(cosmic) if cosmic else "."


def variant_class(row):
    consequence = first_value(row, ["VEP_CONSEQUENCE", "SNPEFF_ANNOTATION"])
    if consequence != ".":
        return consequence
    ref = row.get("REF", "")
    alt = row.get("ALT", "")
    if len(ref) == len(alt):
        return "SNV"
    if len(ref) > len(alt):
        return "DEL"
    return "INS"


def snpeff_for_vep_transcript(row):
    vep_feature = row.get("VEP_FEATURE")
    ann = row.get("SNPEFF_ANN")
    if vep_feature in (None, "", ".") or ann in (None, "", "."):
        return {}
    vep_stem = vep_feature.split(".", 1)[0]
    for entry in ann.split(","):
        parts = entry.split("|")
        if len(parts) < 11:
            continue
        feature_id = parts[6]
        if feature_id.split(".", 1)[0] == vep_stem:
            return {
                "feature_id": feature_id,
                "rank": parts[8],
                "hgvsc": parts[9],
                "hgvsp": parts[10],
            }
    return {}


def cdna_change(row):
    hgvsc = row.get("VEP_HGVSC")
    if hgvsc not in (None, "", "."):
        return hgvsc
    matched = snpeff_for_vep_transcript(row)
    if matched.get("hgvsc"):
        return matched["hgvsc"]
    return first_value(row, ["SNPEFF_HGVS_C"])


def aa_change(row):
    hgvsp = row.get("VEP_HGVSP")
    if hgvsp not in (None, "", "."):
        return hgvsp

    amino_acids = row.get("VEP_AMINO_ACIDS")
    protein_position = row.get("VEP_PROTEIN_POSITION")
    if amino_acids not in (None, "", ".") and protein_position not in (None, "", "."):
        if "/" in amino_acids and "-" not in protein_position:
            ref_aa, alt_aa = amino_acids.split("/", 1)
            ref_aa = AA3.get(ref_aa, ref_aa)
            alt_aa = AA3.get(alt_aa, alt_aa)
            return f"p.{ref_aa}{protein_position}{alt_aa}"

    matched = snpeff_for_vep_transcript(row)
    if matched.get("hgvsp"):
        return matched["hgvsp"]

    return first_value(row, ["SNPEFF_HGVS_P"])


def failure_reason(row):
    if row.get("NOISE_TIER") in {"HIGH_CONFIDENCE", "SINGLE_CALLER", "REVIEW"}:
        return "."
    flags = []
    for field in ["HARD_FILTER_FLAGS", "NOISE_FLAGS"]:
        value = row.get(field)
        if value not in (None, "", ".", "NONE"):
            flags.append(value)
    return ";".join(flags) if flags else "."


def convert_row(row):
    tumor, normal = best_counts(row)
    sample, normal_used = preferred_sample_names(row)
    # Merge annotator IDs (empty until VEP/SnpEff carry a known-variants DB)
    # with IDs from the hotspot-list match, so COSMIC/dbSNP populate now via the
    # hotspot list and automatically gain annotator IDs once those are present.
    dbsnp, cosmic = extract_existing_ids(
        row.get("VEP_EXISTING_VARIATION"), row.get("HOTSPOT_IDS")
    )
    comments = []
    for field in ["RESCUE_FLAGS", "NOISE_FLAGS"]:
        value = row.get(field)
        if value not in (None, "", ".", "NONE"):
            comments.append(f"{field}={value}")

    out = {
        "Sample": clean(sample),
        "NormalUsed": clean(normal_used),
        "Chrom": clean(row.get("CHROM")),
        "Start": clean(row.get("POS")),
        "Ref": clean(row.get("REF")),
        "Alt": clean(row.get("ALT")),
        "VariantClass": variant_class(row),
        "Gene": first_value(row, ["VEP_SYMBOL", "SNPEFF_GENE_NAME"]),
        "Exon": first_value(row, ["VEP_EXON", "SNPEFF_RANK"]),
        "Call_Confidence": clean(row.get("NOISE_TIER")),
        "Comments": ";".join(comments) if comments else ".",
        "TranscriptID": first_value(row, ["VEP_FEATURE", "SNPEFF_FEATURE_ID"]),
        "cDNAchange": cdna_change(row),
        "AAchange": aa_change(row),
        "dbSNP_ID": dbsnp,
        "Cosmic_ID": cosmic,
        "1000G_MAF": ".",
        "FailureReason": failure_reason(row),
        "CallMethod": clean(row.get("EFFECTIVE_CALLERS") or row.get("CALLERS")),
        "COSMIC_site": cosmic,
        "N_TotalDepth": fmt_number(normal.get("dp")),
        "N_RefCount": fmt_number(normal.get("ref")),
        "N_AltCount": fmt_number(normal.get("alt")),
        "N_AltFreq": fmt_number(normal.get("af")),
        "T_TotalDepth": fmt_number(tumor.get("dp")),
        "T_RefCount": fmt_number(tumor.get("ref")),
        "T_AltCount": fmt_number(tumor.get("alt")),
        "T_AltFreq": fmt_number(tumor.get("af")),
        "T_Ref+": fmt_number(tumor.get("ref_plus")),
        "T_Ref-": fmt_number(tumor.get("ref_minus")),
        "T_Alt+": fmt_number(tumor.get("alt_plus")),
        "T_Alt-": fmt_number(tumor.get("alt_minus")),
        "Strand_Bias": first_value(row, ["VARDICT_FORMAT_SBF", "LOFREQ_SNV_INFO_SB"]),
        "All_N_Aggregate_AlleleDepth": ".",
        "All_N_Median_AlleleFreq": ".",
        "T_freq/All_N_Freq": ".",
        "Occurence_in_Normals": ".",
        "Germline": clean(row.get("GERMLINE")),
        "Germline_Info": clean(row.get("GERMLINE_INFO")),
        "Hotspot": clean(row.get("HOTSPOT")),
        "Hotspot_Info": clean(row.get("HOTSPOT_INFO")),
    }
    return {column: clean(out.get(column)) for column in OUTPUT_COLUMNS}


def convert(input_path, output_path):
    with open(input_path, newline="") as in_f, open(output_path, "w", newline="") as out_f:
        reader = csv.DictReader(in_f, delimiter="\t")
        writer = csv.DictWriter(out_f, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        count = 0
        for row in reader:
            writer.writerow(convert_row(row))
            count += 1
    return count


def main():
    parser = argparse.ArgumentParser(
        description="Create a polished clinical-team TSV from the high-confidence variant table."
    )
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {VERSION}")
    parser.add_argument("input", help="High-confidence TSV from tier_merged_calls.py")
    parser.add_argument("-o", "--output", required=True, help="Polished clinical TSV output")
    args = parser.parse_args()
    try:
        count = convert(args.input, args.output)
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)
    print(f"Polished clinical TSV written to {args.output}")
    print(f"Rows: {count}")


if __name__ == "__main__":
    main()
