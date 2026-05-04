#!/usr/bin/env python3
"""
MAPK signaling pathway alteration filter and summarizer.

Purpose:
  - Use a small MAPK gene set, e.g. 26 genes.
  - Keep only alterations whose gene symbol intersects the gene set.
  - Support SNV/InDel, CNV, and SV records.
  - Ignore/drop the SV `fusion` column by design.
  - Emit alteration-level, sample-level, and gene-level TSVs.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

VERSION = "0.2.0"

GENE_COLUMN_ALIASES = [
    "gene_symbol", "gene", "Gene", "GENE", "SYMBOL", "Hugo_Symbol", "HUGO_SYMBOL",
    "gene_name", "Gene.refGene", "Gene_Name", "symbol", "ANN_SYMBOL"
]
SAMPLE_COLUMN_ALIASES = [
    "sample", "Sample", "SAMPLE", "sample_id", "Sample_ID", "tumor_sample", "Tumor_Sample_Barcode",
    "Tumor_Sample", "sample_name"
]
CONSEQUENCE_COLUMN_ALIASES = [
    "consequence", "Consequence", "effect", "Effect", "variant_classification",
    "Variant_Classification", "Func.refGene", "ExonicFunc.refGene", "BIOTYPE", "ANN_CONSEQUENCE"
]
PROTEIN_COLUMN_ALIASES = [
    "protein_change", "Protein_Change", "HGVSp", "hgvsp", "HGVSp_Short",
    "AAChange", "AAChange.refGene", "amino_acid_change", "Protein_position", "ANN_HGVSP"
]
VAF_COLUMN_ALIASES = [
    "vaf", "VAF", "tumor_vaf", "Tumor_VAF", "AF", "allele_fraction", "Allele_Fraction",
    "t_alt_freq", "TUMOR_AF"
]
DEPTH_COLUMN_ALIASES = [
    "depth", "Depth", "DP", "tumor_depth", "Tumor_Depth", "t_depth", "TUMOR_DP"
]
ALT_DEPTH_COLUMN_ALIASES = ["alt_depth", "Alt_Depth", "VD", "AD", "tumor_alt_count", "t_alt_count"]
REF_COLUMN_ALIASES = ["ref", "REF", "Reference_Allele", "Reference", "ReferenceAllele"]
ALT_COLUMN_ALIASES = ["alt", "ALT", "Tumor_Seq_Allele2", "Alternate_Allele", "Alternate", "Allele"]
CHROM_COLUMN_ALIASES = ["chrom", "chr", "CHROM", "Chromosome", "chromosome"]
POS_COLUMN_ALIASES = ["pos", "POS", "start", "Start_Position", "position", "Position"]
CALLER_COLUMN_ALIASES = ["caller", "Caller", "variant_caller", "Variant_Caller", "SUPP_VEC"]

# CNV aliases matching DMND*_CNV_MERGED.tsv, while staying flexible.
CNV_GENE_COLUMN_ALIASES = ["gene", "Gene", "GENE", "gene_symbol", "ANN_SYMBOL"]
CNV_SAMPLE_COLUMN_ALIASES = SAMPLE_COLUMN_ALIASES
CNV_CALL_COLUMN_ALIASES = ["call", "cnv_call", "CALL", "status", "seg_status"]
CNV_CN_COLUMN_ALIASES = ["cn", "CN", "copy_number", "Copy_Number", "tcn"]
CNV_LOG2_COLUMN_ALIASES = ["log2", "LOG2", "log2_ratio", "cnlr_median"]
CNV_TOOL_COLUMN_ALIASES = ["tool", "caller", "source", "source_file"]
CNV_CHROM_COLUMN_ALIASES = ["chrom", "chr", "CHROM", "Chromosome", "chromosome"]
CNV_START_COLUMN_ALIASES = ["start", "START", "pos", "POS", "Start_Position"]
CNV_END_COLUMN_ALIASES = ["end", "END", "End_Position"]

# SV aliases matching DMND*_SOMTIC_SV_ANN.tsv.
SV_GENE1_COLUMN_ALIASES = ["gene1", "Gene1", "GENE1", "gene_1"]
SV_GENE2_COLUMN_ALIASES = ["gene2", "Gene2", "GENE2", "gene_2"]
SV_SAMPLE_COLUMN_ALIASES = SAMPLE_COLUMN_ALIASES
SVTYPE_COLUMN_ALIASES = ["SVTYPE", "svtype", "sv_type", "type"]
SV_CALLER_COLUMN_ALIASES = ["SUPP_VEC", "caller", "Caller", "SVMETHOD"]
SV_FILTER_COLUMN_ALIASES = ["FILTER", "filter"]
SV_SUPP_COLUMN_ALIASES = ["SUPP", "supp"]
SV_PE_COLUMN_ALIASES = ["PE", "pe"]
SV_SR_COLUMN_ALIASES = ["SR", "sr"]
SV_DP_COLUMN_ALIASES = ["DP", "depth"]
SV_SITE1_COLUMN_ALIASES = ["site1", "Site1", "SITE1"]
SV_SITE2_COLUMN_ALIASES = ["site2", "Site2", "SITE2"]
SV_CHR1_COLUMN_ALIASES = ["chr1", "CHROM", "chrom1"]
SV_POS1_COLUMN_ALIASES = ["pos1", "POS", "start1"]
SV_CHR2_COLUMN_ALIASES = ["chr2", "CHR2", "chrom2"]
SV_POS2_COLUMN_ALIASES = ["pos2", "END", "end2"]

# Cancer-relevant MAPK hotspot pattern, intentionally conservative.
HOTSPOT_PATTERNS = {
    "KRAS": [r"G12", r"G13", r"Q61", r"A146"],
    "NRAS": [r"G12", r"G13", r"Q61"],
    "HRAS": [r"G12", r"G13", r"Q61"],
    "BRAF": [r"V600", r"K601", r"G469", r"D594", r"L597"],
    "MAP2K1": [r"K57", r"Q56", r"I99", r"P124", r"C121"],
    "MAP2K2": [r"K57", r"Q60", r"I103", r"P128"],
    "EGFR": [r"L858", r"T790", r"G719", r"L861", r"S768", r"E746", r"L747", r"E709"],
    "FGFR1": [r"N546", r"K656"],
    "FGFR2": [r"S252", r"P253", r"N549", r"K659"],
    "FGFR3": [r"S249", r"Y373", r"R248", r"K650", r"G370"],
    "FGFR4": [r"N535", r"V550"],
    "NTRK1": [r"G595", r"G667"],
    "NTRK2": [r"G639"],
    "NTRK3": [r"G623"],
    "PDGFRA": [r"D842", r"V561", r"N659"],
    "RAC1": [r"P29"],
    "AKT1": [r"E17"],
}

TUMOR_SUPPRESSOR_OR_NEGATIVE_REGULATOR = {"NF1", "TP53", "TGFBR1", "TGFBR2"}
ONCOGENE_LIKE = {
    "AKT1", "BRAF", "EGFR", "FGFR1", "FGFR2", "FGFR3", "FGFR4", "HRAS", "KRAS",
    "MAP2K1", "MAP2K2", "MYC", "NTRK1", "NTRK2", "NTRK3", "PDGFRA", "RAC1", "RAF1", "RRAS2", "SOS1"
}

LOF_TERMS = [
    "frameshift", "frame_shift", "stop_gained", "nonsense", "stopgain",
    "splice_acceptor", "splice_donor", "splice_site", "start_lost", "startloss", "nonstop",
    "stop_lost"
]
MISSENSE_TERMS = ["missense", "nonsynonymous", "non_synonymous", "missense_variant"]
INFRAME_TERMS = ["inframe", "in-frame", "in_frame", "nonframeshift"]
SILENT_TERMS = ["synonymous", "silent", "utr", "intron", "intergenic", "upstream", "downstream"]

NEUTRAL_CNV_CALLS = {"", ".", "NA", "NAN", "NONE", "NULL", "NEUTRAL", "NEUTR", "DIPLOID"}
AMP_CNV_CALLS = {"AMP", "AMPLIFICATION", "AMPLIFIED", "HIGH_GAIN", "GAIN_HIGH"}
GAIN_CNV_CALLS = {"GAIN", "LOW_GAIN", "COPY_GAIN"}
LOSS_CNV_CALLS = {"LOSS", "DEL", "DELETION", "COPY_LOSS"}
HOMDEL_CNV_CALLS = {"HOMDEL", "HOM_DEL", "HOMODELETION", "DEEP_DEL", "DEEP_DELETION"}

SVTYPE_DETAIL = {
    "TRA": "translocation",
    "BND": "breakend",
    "INV": "inversion",
    "DEL": "deletion",
    "DUP": "duplication",
    "INS": "insertion",
}

ALTERATION_COLUMNS = [
    "sample", "alteration_type", "gene", "partner_gene", "event_class", "confidence", "detail",
    "chrom", "start", "end", "ref", "alt", "consequence", "protein_change",
    "vaf", "alt_depth", "depth", "cnv_call", "copy_number", "log2", "svtype", "site1", "site2",
    "caller", "filter", "support", "pathway", "dataset", "note",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Filter and summarize MAPK SNV/InDel, CNV, and SV calls in a restricted gene set.")
    parser.add_argument("--sample-id", required=True, help="Sample ID used if an input table lacks a sample column.")
    parser.add_argument("--variants-tsv", "--snvindel-tsv", dest="variants_tsv", required=True, type=Path, help="Input SNV/InDel variant TSV.")
    parser.add_argument("--cnv-tsv", default=None, type=Path, help="Optional CNV TSV, e.g. DMND*_CNV_MERGED.tsv.")
    parser.add_argument("--sv-tsv", default=None, type=Path, help="Optional SV annotation TSV, e.g. DMND*_SOMTIC_SV_ANN.tsv.")
    parser.add_argument("--mapk-genes", required=True, type=Path, help="MAPK gene-set TSV with gene_symbol or gene column.")
    parser.add_argument("--out-prefix", required=True, help="Output prefix.")
    parser.add_argument("--min-vaf", type=float, default=None, help="Optional minimum VAF filter for SNV/InDel rows.")
    parser.add_argument("--min-depth", type=float, default=None, help="Optional minimum depth filter for SNV/InDel rows.")
    parser.add_argument("--include-neutral-cnv", action="store_true", help="Emit NEUTRAL CNV rows. Default: exclude them.")
    parser.add_argument("--version", action="version", version=f"mapk_sig_path {VERSION}")
    return parser.parse_args()


def valid_optional_path(path: Optional[Path]) -> bool:
    if path is None:
        return False
    text = str(path).strip()
    if text == "" or text.lower() in {"null", "none", "na", "."}:
        return False
    return path.exists() and path.is_file()


def read_tsv(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"No header found in {path}")
        return list(reader.fieldnames), [dict(row) for row in reader]


def find_col(columns: Sequence[str], aliases: Sequence[str]) -> Optional[str]:
    exact = {c: c for c in columns}
    lower = {c.lower(): c for c in columns}
    for alias in aliases:
        if alias in exact:
            return exact[alias]
        if alias.lower() in lower:
            return lower[alias.lower()]
    return None


def clean_gene(value: str) -> str:
    value = (value or "").strip()
    if not value or value in {".", "-", "NA", "N/A", "nan", "NaN"}:
        return ""
    value = value.replace(" ", "")
    return value.upper()


def expand_gene_field(value: str) -> List[str]:
    value = (value or "").strip()
    if not value or value in {".", "-"}:
        return []
    parts = re.split(r"[;,|/]", value)
    seen = set()
    genes = []
    for part in parts:
        gene = clean_gene(part)
        if gene and gene not in seen:
            seen.add(gene)
            genes.append(gene)
    return genes


def to_float(value: str) -> Optional[float]:
    if value is None:
        return None
    value = str(value).strip()
    if value == "" or value.upper() in {"NA", "NAN", ".", "NONE", "NULL"}:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def read_gene_set(path: Path) -> Dict[str, Dict[str, str]]:
    columns, rows = read_tsv(path)
    gene_col = find_col(columns, ["gene_symbol", "gene", "symbol", "Gene", "SYMBOL"])
    if gene_col is None:
        raise ValueError(f"Could not find a gene column in {path}. Columns: {', '.join(columns)}")

    gene_info: Dict[str, Dict[str, str]] = {}
    for row in rows:
        gene = clean_gene(row.get(gene_col, ""))
        if not gene:
            continue
        gene_info[gene] = {
            "gene": gene,
            "pathway": row.get("pathway", "MAPK signaling pathway"),
            "dataset": row.get("dataset", "custom_MAPK_gene_set"),
            "threshold_value": row.get("threshold_value", ""),
        }
    return gene_info


def infer_variant_type(row: Dict[str, str], ref_col: Optional[str], alt_col: Optional[str], consequence: str) -> str:
    ref = row.get(ref_col, "") if ref_col else ""
    alt = row.get(alt_col, "") if alt_col else ""
    ref = ref.strip()
    alt = alt.strip()

    if ref and alt and ref not in {".", "-"} and alt not in {".", "-"}:
        if len(ref) == 1 and len(alt) == 1:
            return "SNV"
        return "InDel"

    text = consequence.lower()
    if any(term in text for term in ["del", "ins", "indel", "frameshift", "inframe", "nonframeshift"]):
        return "InDel"
    return "SNV/InDel"


def is_hotspot(gene: str, protein_change: str) -> bool:
    protein_change = protein_change or ""
    patterns = HOTSPOT_PATTERNS.get(gene, [])
    return any(re.search(pattern, protein_change, flags=re.IGNORECASE) for pattern in patterns)


def classify_snvindel_event(gene: str, consequence: str, protein_change: str) -> Tuple[str, str, str, str]:
    text = f"{consequence} {protein_change}".lower()

    if is_hotspot(gene, protein_change):
        return "known_mapk_hotspot", "high", "hotspot", "Known recurrent MAPK/cancer hotspot pattern"

    if any(term in text for term in LOF_TERMS):
        if gene in TUMOR_SUPPRESSOR_OR_NEGATIVE_REGULATOR:
            return "loss_of_function", "possible", "lof", "LoF event in tumor suppressor / negative regulator"
        return "loss_of_function", "uncertain", "lof", "LoF event; interpret based on gene biology"

    if any(term in text for term in INFRAME_TERMS):
        if gene in ONCOGENE_LIKE:
            return "inframe_indel", "possible", "inframe_indel", "In-frame event in MAPK-associated oncogene-like gene"
        return "inframe_indel", "uncertain", "inframe_indel", "In-frame event in MAPK gene"

    if any(term in text for term in MISSENSE_TERMS):
        if gene in ONCOGENE_LIKE:
            return "missense_non_hotspot", "uncertain", "missense", "Non-hotspot missense in MAPK-associated oncogene-like gene"
        if gene in TUMOR_SUPPRESSOR_OR_NEGATIVE_REGULATOR:
            return "missense_non_hotspot", "uncertain", "missense", "Non-hotspot missense in tumor suppressor / negative regulator"
        return "missense_non_hotspot", "uncertain", "missense", "Non-hotspot missense"

    if any(term in text for term in SILENT_TERMS):
        return "silent_or_non_coding", "low", "silent_or_non_coding", "Synonymous/non-coding annotation"

    return "other_snv_indel", "uncertain", "other", "SNV/InDel event with unclassified consequence"


def classify_cnv_event(gene: str, call: str) -> Tuple[str, str, str, str]:
    normalized = (call or "").strip().upper()
    if normalized in AMP_CNV_CALLS:
        conf = "high" if gene in ONCOGENE_LIKE else "possible"
        return "amplification", conf, "amplification", "CNV amplification in MAPK gene"
    if normalized in GAIN_CNV_CALLS:
        conf = "possible" if gene in ONCOGENE_LIKE else "uncertain"
        return "copy_number_gain", conf, "gain", "CNV gain in MAPK gene"
    if normalized in HOMDEL_CNV_CALLS:
        conf = "high" if gene in TUMOR_SUPPRESSOR_OR_NEGATIVE_REGULATOR else "possible"
        return "deep_deletion", conf, "deep_deletion", "Deep deletion in MAPK gene"
    if normalized in LOSS_CNV_CALLS:
        conf = "possible" if gene in TUMOR_SUPPRESSOR_OR_NEGATIVE_REGULATOR else "uncertain"
        return "copy_number_loss", conf, "loss", "CNV loss in MAPK gene"
    return "copy_number_alteration", "uncertain", normalized.lower() or "cnv", "CNV call in MAPK gene"


def classify_sv_event(gene: str, partner_gene: str, svtype: str, site1: str, site2: str) -> Tuple[str, str, str, str]:
    svtype_norm = (svtype or "").strip().upper()
    detail = SVTYPE_DETAIL.get(svtype_norm, svtype_norm.lower() or "structural_variant")
    partner = clean_gene(partner_gene)

    if partner and partner in ONCOGENE_LIKE | TUMOR_SUPPRESSOR_OR_NEGATIVE_REGULATOR:
        event_class = "mapk_pathway_sv"
    elif partner and partner == gene:
        event_class = "mapk_intragenic_sv"
    else:
        event_class = "mapk_partner_sv"

    # Conservative: this is not fusion logic and does not use the fusion column.
    if svtype_norm in {"TRA", "BND", "INV", "DEL", "DUP"}:
        confidence = "possible"
    else:
        confidence = "uncertain"

    note = "MAPK-associated structural variant; fusion annotation intentionally ignored"
    if site1 or site2:
        note += f"; breakpoint context: {site1 or 'NA'} / {site2 or 'NA'}"
    return event_class, confidence, detail, note


def empty_alteration_row() -> Dict[str, str]:
    return {col: "" for col in ALTERATION_COLUMNS}


def parse_snvindel_rows(args: argparse.Namespace, gene_info: Dict[str, Dict[str, str]]) -> List[Dict[str, str]]:
    mapk_genes = set(gene_info)
    columns, rows = read_tsv(args.variants_tsv)

    gene_col = find_col(columns, GENE_COLUMN_ALIASES)
    if gene_col is None:
        raise ValueError(
            "Could not find a gene-symbol column in SNV/InDel TSV. "
            f"Supported aliases: {', '.join(GENE_COLUMN_ALIASES)}. "
            f"Observed columns: {', '.join(columns)}"
        )

    sample_col = find_col(columns, SAMPLE_COLUMN_ALIASES)
    consequence_col = find_col(columns, CONSEQUENCE_COLUMN_ALIASES)
    protein_col = find_col(columns, PROTEIN_COLUMN_ALIASES)
    vaf_col = find_col(columns, VAF_COLUMN_ALIASES)
    depth_col = find_col(columns, DEPTH_COLUMN_ALIASES)
    alt_depth_col = find_col(columns, ALT_DEPTH_COLUMN_ALIASES)
    ref_col = find_col(columns, REF_COLUMN_ALIASES)
    alt_col = find_col(columns, ALT_COLUMN_ALIASES)
    chrom_col = find_col(columns, CHROM_COLUMN_ALIASES)
    pos_col = find_col(columns, POS_COLUMN_ALIASES)
    caller_col = find_col(columns, CALLER_COLUMN_ALIASES)
    filter_col = find_col(columns, ["FILTER", "filter"])

    alteration_rows: List[Dict[str, str]] = []
    for row in rows:
        genes = expand_gene_field(row.get(gene_col, ""))
        hits = [g for g in genes if g in mapk_genes]
        if not hits:
            continue

        vaf = to_float(row.get(vaf_col, "")) if vaf_col else None
        depth = to_float(row.get(depth_col, "")) if depth_col else None

        if args.min_vaf is not None and vaf is not None and vaf < args.min_vaf:
            continue
        if args.min_depth is not None and depth is not None and depth < args.min_depth:
            continue

        sample = row.get(sample_col, args.sample_id) if sample_col else args.sample_id
        consequence = row.get(consequence_col, "") if consequence_col else ""
        protein_change = row.get(protein_col, "") if protein_col else ""
        variant_type = infer_variant_type(row, ref_col, alt_col, consequence)

        for gene in hits:
            event_class, confidence, detail, note = classify_snvindel_event(gene, consequence, protein_change)
            info = gene_info.get(gene, {})
            out = empty_alteration_row()
            out.update({
                "sample": sample,
                "alteration_type": variant_type,
                "gene": gene,
                "partner_gene": "",
                "event_class": event_class,
                "confidence": confidence,
                "detail": detail,
                "chrom": row.get(chrom_col, "") if chrom_col else "",
                "start": row.get(pos_col, "") if pos_col else "",
                "end": row.get(pos_col, "") if pos_col else "",
                "ref": row.get(ref_col, "") if ref_col else "",
                "alt": row.get(alt_col, "") if alt_col else "",
                "consequence": consequence,
                "protein_change": protein_change,
                "vaf": "" if vaf is None else str(vaf),
                "alt_depth": row.get(alt_depth_col, "") if alt_depth_col else "",
                "depth": "" if depth is None else str(depth),
                "caller": row.get(caller_col, "") if caller_col else "",
                "filter": row.get(filter_col, "") if filter_col else "",
                "pathway": info.get("pathway", "MAPK signaling pathway"),
                "dataset": info.get("dataset", "custom_MAPK_gene_set"),
                "note": note,
            })
            alteration_rows.append(out)
    return alteration_rows


def parse_cnv_rows(args: argparse.Namespace, gene_info: Dict[str, Dict[str, str]]) -> List[Dict[str, str]]:
    if not valid_optional_path(args.cnv_tsv):
        return []

    mapk_genes = set(gene_info)
    columns, rows = read_tsv(args.cnv_tsv)
    gene_col = find_col(columns, CNV_GENE_COLUMN_ALIASES)
    if gene_col is None:
        raise ValueError(
            "Could not find a gene column in CNV TSV. "
            f"Supported aliases: {', '.join(CNV_GENE_COLUMN_ALIASES)}. "
            f"Observed columns: {', '.join(columns)}"
        )

    sample_col = find_col(columns, CNV_SAMPLE_COLUMN_ALIASES)
    call_col = find_col(columns, CNV_CALL_COLUMN_ALIASES)
    cn_col = find_col(columns, CNV_CN_COLUMN_ALIASES)
    log2_col = find_col(columns, CNV_LOG2_COLUMN_ALIASES)
    tool_col = find_col(columns, CNV_TOOL_COLUMN_ALIASES)
    chrom_col = find_col(columns, CNV_CHROM_COLUMN_ALIASES)
    start_col = find_col(columns, CNV_START_COLUMN_ALIASES)
    end_col = find_col(columns, CNV_END_COLUMN_ALIASES)

    # Collapse repeated gene-level rows from callers such as OncoCNV so one gene call
    # is counted once even when the file contains multiple exon/probe intervals.
    collapsed: Dict[Tuple[str, str, str, str, str, str, str], Dict[str, object]] = {}

    for row in rows:
        genes = expand_gene_field(row.get(gene_col, ""))
        hits = [g for g in genes if g in mapk_genes]
        if not hits:
            continue

        call = row.get(call_col, "") if call_col else ""
        normalized_call = (call or "").strip().upper()
        if not args.include_neutral_cnv and normalized_call in NEUTRAL_CNV_CALLS:
            continue

        sample = row.get(sample_col, args.sample_id) if sample_col else args.sample_id
        caller = row.get(tool_col, "") if tool_col else ""
        chrom = row.get(chrom_col, "") if chrom_col else ""
        cn = row.get(cn_col, "") if cn_col else ""
        log2 = row.get(log2_col, "") if log2_col else ""
        start_val = to_float(row.get(start_col, "")) if start_col else None
        end_val = to_float(row.get(end_col, "")) if end_col else None

        for gene in hits:
            key = (sample, gene, normalized_call, caller, chrom, cn, log2)
            if key not in collapsed:
                collapsed[key] = {
                    "sample": sample,
                    "gene": gene,
                    "call": call,
                    "caller": caller,
                    "chrom": chrom,
                    "cn": cn,
                    "log2": log2,
                    "starts": [],
                    "ends": [],
                    "n_regions": 0,
                }
            item = collapsed[key]
            item["n_regions"] = int(item["n_regions"]) + 1
            if start_val is not None:
                item["starts"].append(start_val)  # type: ignore[index, union-attr]
            if end_val is not None:
                item["ends"].append(end_val)  # type: ignore[index, union-attr]

    alteration_rows: List[Dict[str, str]] = []
    for item in collapsed.values():
        gene = str(item["gene"])
        call = str(item["call"])
        event_class, confidence, detail, note = classify_cnv_event(gene, call)
        n_regions = int(item["n_regions"])
        if n_regions > 1:
            note = f"{note}; collapsed {n_regions} CNV rows for this sample/gene/caller/call"
        starts = item["starts"]  # type: ignore[assignment]
        ends = item["ends"]  # type: ignore[assignment]
        start = str(int(min(starts))) if starts else ""
        end = str(int(max(ends))) if ends else ""
        info = gene_info.get(gene, {})
        out = empty_alteration_row()
        out.update({
            "sample": str(item["sample"]),
            "alteration_type": "CNV",
            "gene": gene,
            "partner_gene": "",
            "event_class": event_class,
            "confidence": confidence,
            "detail": detail,
            "chrom": str(item["chrom"]),
            "start": start,
            "end": end,
            "cnv_call": call,
            "copy_number": str(item["cn"]),
            "log2": str(item["log2"]),
            "caller": str(item["caller"]),
            "support": f"n_regions={n_regions}",
            "pathway": info.get("pathway", "MAPK signaling pathway"),
            "dataset": info.get("dataset", "custom_MAPK_gene_set"),
            "note": note,
        })
        alteration_rows.append(out)
    return alteration_rows


def norm_chr(value: str) -> str:
    value = (value or "").strip()
    value = re.sub(r"^chr", "", value, flags=re.IGNORECASE)
    return value


def sv_endpoint(row: Dict[str, str], chr_col: Optional[str], pos_col: Optional[str]) -> str:
    chrom = norm_chr(row.get(chr_col, "") if chr_col else "")
    pos = (row.get(pos_col, "") if pos_col else "").strip()
    return f"{chrom}:{pos}"


def sv_dedup_key(row: Dict[str, str], chr1_col: Optional[str], pos1_col: Optional[str], chr2_col: Optional[str], pos2_col: Optional[str], svtype_col: Optional[str]) -> Tuple[str, str, str]:
    ep1 = sv_endpoint(row, chr1_col, pos1_col)
    ep2 = sv_endpoint(row, chr2_col, pos2_col)
    svtype = (row.get(svtype_col, "") if svtype_col else "").strip().upper()
    a, b = sorted([ep1, ep2])
    return (a, b, svtype)


def parse_sv_rows(args: argparse.Namespace, gene_info: Dict[str, Dict[str, str]]) -> List[Dict[str, str]]:
    if not valid_optional_path(args.sv_tsv):
        return []

    mapk_genes = set(gene_info)
    columns, rows = read_tsv(args.sv_tsv)

    gene1_col = find_col(columns, SV_GENE1_COLUMN_ALIASES)
    gene2_col = find_col(columns, SV_GENE2_COLUMN_ALIASES)
    if gene1_col is None or gene2_col is None:
        raise ValueError(
            "Could not find gene1/gene2 columns in SV TSV. "
            f"Observed columns: {', '.join(columns)}"
        )

    sample_col = find_col(columns, SV_SAMPLE_COLUMN_ALIASES)
    svtype_col = find_col(columns, SVTYPE_COLUMN_ALIASES)
    caller_col = find_col(columns, SV_CALLER_COLUMN_ALIASES)
    filter_col = find_col(columns, SV_FILTER_COLUMN_ALIASES)
    supp_col = find_col(columns, SV_SUPP_COLUMN_ALIASES)
    pe_col = find_col(columns, SV_PE_COLUMN_ALIASES)
    sr_col = find_col(columns, SV_SR_COLUMN_ALIASES)
    dp_col = find_col(columns, SV_DP_COLUMN_ALIASES)
    site1_col = find_col(columns, SV_SITE1_COLUMN_ALIASES)
    site2_col = find_col(columns, SV_SITE2_COLUMN_ALIASES)
    chr1_col = find_col(columns, SV_CHR1_COLUMN_ALIASES)
    pos1_col = find_col(columns, SV_POS1_COLUMN_ALIASES)
    chr2_col = find_col(columns, SV_CHR2_COLUMN_ALIASES)
    pos2_col = find_col(columns, SV_POS2_COLUMN_ALIASES)

    seen_keys = set()
    alteration_rows: List[Dict[str, str]] = []

    for row in rows:
        genes1 = expand_gene_field(row.get(gene1_col, ""))
        genes2 = expand_gene_field(row.get(gene2_col, ""))
        hits1 = [g for g in genes1 if g in mapk_genes]
        hits2 = [g for g in genes2 if g in mapk_genes]
        hits = []
        for g in hits1 + hits2:
            if g not in hits:
                hits.append(g)
        if not hits:
            continue

        dedup_key = sv_dedup_key(row, chr1_col, pos1_col, chr2_col, pos2_col, svtype_col)
        if dedup_key in seen_keys:
            continue
        seen_keys.add(dedup_key)

        sample = row.get(sample_col, args.sample_id) if sample_col else args.sample_id
        svtype = row.get(svtype_col, "") if svtype_col else ""
        site1 = row.get(site1_col, "") if site1_col else ""
        site2 = row.get(site2_col, "") if site2_col else ""
        chrom = norm_chr(row.get(chr1_col, "") if chr1_col else "")
        start = row.get(pos1_col, "") if pos1_col else ""
        end = row.get(pos2_col, "") if pos2_col else ""
        caller = row.get(caller_col, "") if caller_col else ""
        filt = row.get(filter_col, "") if filter_col else ""
        support_parts = []
        for label, col in [("SUPP", supp_col), ("PE", pe_col), ("SR", sr_col), ("DP", dp_col)]:
            if col:
                value = row.get(col, "")
                if value not in {"", ".", "NA"}:
                    support_parts.append(f"{label}={value}")
        support = ";".join(support_parts)

        for gene in hits:
            if gene in hits1:
                partner_candidates = genes2 or [""]
                gene_site = site1
                partner_site = site2
            else:
                partner_candidates = genes1 or [""]
                gene_site = site2
                partner_site = site1
            partner_gene = ",".join([g for g in partner_candidates if g and g != gene])

            event_class, confidence, detail, note = classify_sv_event(gene, partner_gene, svtype, gene_site, partner_site)
            info = gene_info.get(gene, {})
            out = empty_alteration_row()
            out.update({
                "sample": sample,
                "alteration_type": "SV",
                "gene": gene,
                "partner_gene": partner_gene,
                "event_class": event_class,
                "confidence": confidence,
                "detail": detail,
                "chrom": chrom,
                "start": start,
                "end": end,
                "svtype": svtype,
                "site1": gene_site,
                "site2": partner_site,
                "caller": caller,
                "filter": filt,
                "support": support,
                "pathway": info.get("pathway", "MAPK signaling pathway"),
                "dataset": info.get("dataset", "custom_MAPK_gene_set"),
                "note": note,
            })
            alteration_rows.append(out)
    return alteration_rows


def write_tsv(path: Path, rows: List[Dict[str, str]], columns: Sequence[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(columns), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def summarize(alteration_rows: List[Dict[str, str]], sample_id: str, out_prefix: str) -> None:
    by_sample: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    by_sample_gene: Dict[Tuple[str, str], List[Dict[str, str]]] = defaultdict(list)
    for row in alteration_rows:
        by_sample[row["sample"]].append(row)
        by_sample_gene[(row["sample"], row["gene"])].append(row)

    samples = sorted(set(by_sample) | {sample_id})
    sample_summary_rows: List[Dict[str, str]] = []
    confidence_rank = {"high": 3, "possible": 2, "uncertain": 1, "low": 0, "": -1}

    for sample in samples:
        events = by_sample.get(sample, [])
        n_high = sum(1 for e in events if e["confidence"] == "high")
        n_possible = sum(1 for e in events if e["confidence"] == "possible")
        n_uncertain = sum(1 for e in events if e["confidence"] == "uncertain")
        n_low = sum(1 for e in events if e["confidence"] == "low")
        genes_hit = sorted({e["gene"] for e in events})
        alteration_counts = Counter(e["alteration_type"] for e in events)

        if n_high > 0:
            status = "High-confidence MAPK alteration"
        elif n_possible > 0:
            status = "Possible MAPK alteration"
        elif n_uncertain > 0 or n_low > 0:
            status = "MAPK variant(s), uncertain significance"
        else:
            status = "No MAPK alteration detected"

        best = ""
        if events:
            def best_key(e: Dict[str, str]) -> Tuple[int, float]:
                vaf = to_float(e.get("vaf", ""))
                log2 = to_float(e.get("log2", ""))
                magnitude = vaf if vaf is not None else (abs(log2) if log2 is not None else 0.0)
                return (confidence_rank.get(e["confidence"], -1), magnitude)

            best_event = sorted(events, key=best_key, reverse=True)[0]
            if best_event["alteration_type"] in {"SNV", "InDel", "SNV/InDel"}:
                suffix = best_event["protein_change"] or best_event["event_class"]
            elif best_event["alteration_type"] == "CNV":
                suffix = best_event["cnv_call"] or best_event["event_class"]
            else:
                suffix = f"{best_event['svtype']} {best_event['partner_gene']}".strip()
            best = f"{best_event['gene']} {suffix}".strip()

        sample_summary_rows.append({
            "sample": sample,
            "n_mapk_events": str(len(events)),
            "n_mapk_genes_hit": str(len(genes_hit)),
            "n_snvindel_events": str(sum(v for k, v in alteration_counts.items() if k in {"SNV", "InDel", "SNV/InDel"})),
            "n_cnv_events": str(alteration_counts.get("CNV", 0)),
            "n_sv_events": str(alteration_counts.get("SV", 0)),
            "n_high_confidence": str(n_high),
            "n_possible": str(n_possible),
            "n_uncertain": str(n_uncertain),
            "n_low": str(n_low),
            "mapk_status": status,
            "top_candidate": best,
            "genes_hit": ",".join(genes_hit),
        })

    sample_summary_columns = [
        "sample", "n_mapk_events", "n_mapk_genes_hit", "n_snvindel_events", "n_cnv_events", "n_sv_events",
        "n_high_confidence", "n_possible", "n_uncertain", "n_low", "mapk_status", "top_candidate", "genes_hit"
    ]
    write_tsv(Path(f"{out_prefix}.mapk_sig_path.sample_summary.tsv"), sample_summary_rows, sample_summary_columns)

    gene_summary_rows: List[Dict[str, str]] = []
    for (sample, gene), events in sorted(by_sample_gene.items()):
        event_counts = Counter(e["event_class"] for e in events)
        alteration_counts = Counter(e["alteration_type"] for e in events)
        best_event = sorted(
            events,
            key=lambda e: confidence_rank.get(e["confidence"], -1),
            reverse=True,
        )[0]
        gene_summary_rows.append({
            "sample": sample,
            "gene": gene,
            "n_events": str(len(events)),
            "alteration_types": ",".join(f"{k}:{v}" for k, v in sorted(alteration_counts.items())),
            "event_classes": ",".join(f"{k}:{v}" for k, v in sorted(event_counts.items())),
            "best_event_class": best_event["event_class"],
            "best_confidence": best_event["confidence"],
            "best_detail": best_event["detail"],
            "best_protein_change": best_event["protein_change"],
            "best_cnv_call": best_event["cnv_call"],
            "best_svtype": best_event["svtype"],
        })

    gene_summary_columns = [
        "sample", "gene", "n_events", "alteration_types", "event_classes", "best_event_class",
        "best_confidence", "best_detail", "best_protein_change", "best_cnv_call", "best_svtype"
    ]
    write_tsv(Path(f"{out_prefix}.mapk_sig_path.gene_summary.tsv"), gene_summary_rows, gene_summary_columns)


def main() -> int:
    args = parse_args()

    gene_info = read_gene_set(args.mapk_genes)
    if not gene_info:
        raise ValueError(f"No genes found in {args.mapk_genes}")

    alteration_rows: List[Dict[str, str]] = []
    alteration_rows.extend(parse_snvindel_rows(args, gene_info))
    alteration_rows.extend(parse_cnv_rows(args, gene_info))
    alteration_rows.extend(parse_sv_rows(args, gene_info))

    write_tsv(Path(f"{args.out_prefix}.mapk_sig_path.alterations.tsv"), alteration_rows, ALTERATION_COLUMNS)
    summarize(alteration_rows, args.sample_id, args.out_prefix)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
