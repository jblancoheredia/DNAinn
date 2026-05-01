#!/usr/bin/env python3
"""
MAPK SNV/InDel filter and summarizer.

Purpose:
  - Use a small MAPK gene set, e.g. 26 genes.
  - Keep only SNV/InDel calls whose gene symbol intersects the gene set.
  - Classify variant consequence into simple, transparent buckets.
  - Emit alteration-level, sample-level, and gene-level TSVs.

This intentionally avoids CNV/SV/fusion logic.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

VERSION = "0.1.0"

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
REF_COLUMN_ALIASES = ["ref", "REF", "Reference_Allele", "Reference", "ReferenceAllele"]
ALT_COLUMN_ALIASES = ["alt", "ALT", "Tumor_Seq_Allele2", "Alternate_Allele", "Alternate", "Allele"]
CHROM_COLUMN_ALIASES = ["chrom", "chr", "CHROM", "Chromosome", "chromosome"]
POS_COLUMN_ALIASES = ["pos", "POS", "start", "Start_Position", "position", "Position"]
CALLER_COLUMN_ALIASES = ["caller", "Caller", "variant_caller", "Variant_Caller", "SUPP_VEC"]

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

TUMOR_SUPPRESSOR_OR_NEGATIVE_REGULATOR = {"NF1", "TP53"}
ONCOGENE_LIKE = {
    "AKT1", "BRAF", "EGFR", "FGFR1", "FGFR2", "FGFR3", "FGFR4", "HRAS", "KRAS",
    "MAP2K1", "MAP2K2", "MYC", "NTRK1", "NTRK2", "PDGFRA", "RAC1", "RAF1", "RRAS2", "SOS1"
}

LOF_TERMS = [
    "frameshift", "frame_shift", "stop_gained", "nonsense", "stopgain", "stop_gained",
    "splice_acceptor", "splice_donor", "splice_site", "start_lost", "startloss", "nonstop",
    "stop_lost"
]
MISSENSE_TERMS = ["missense", "nonsynonymous", "non_synonymous", "missense_variant"]
INFRAME_TERMS = ["inframe", "in-frame", "in_frame", "nonframeshift"]
SILENT_TERMS = ["synonymous", "silent", "utr", "intron", "intergenic", "upstream", "downstream"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Filter and summarize SNV/InDel calls in a restricted MAPK gene set.")
    parser.add_argument("--sample-id", required=True, help="Sample ID used if the variant table lacks a sample column.")
    parser.add_argument("--variants-tsv", required=True, type=Path, help="Input SNV/InDel variant TSV.")
    parser.add_argument("--mapk-genes", required=True, type=Path, help="MAPK gene-set TSV with gene_symbol or gene column.")
    parser.add_argument("--out-prefix", required=True, help="Output prefix.")
    parser.add_argument("--min-vaf", type=float, default=None, help="Optional minimum VAF filter.")
    parser.add_argument("--min-depth", type=float, default=None, help="Optional minimum depth filter.")
    parser.add_argument("--version", action="version", version=f"mapk_snvindel {VERSION}")
    return parser.parse_args()


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
    if not value:
        return ""
    # Some annotations contain multiple genes separated by punctuation.
    # Keep exact symbol cleaning conservative here; expansion happens separately.
    value = value.replace(" ", "")
    return value.upper()


def expand_gene_field(value: str) -> List[str]:
    value = (value or "").strip()
    if not value:
        return []
    parts = re.split(r"[;,|/]", value)
    return [clean_gene(p) for p in parts if clean_gene(p)]


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


def classify_event(gene: str, consequence: str, protein_change: str) -> Tuple[str, str, str]:
    text = f"{consequence} {protein_change}".lower()

    if is_hotspot(gene, protein_change):
        return "known_mapk_hotspot", "high", "Known recurrent MAPK/cancer hotspot pattern"

    if any(term in text for term in LOF_TERMS):
        if gene in TUMOR_SUPPRESSOR_OR_NEGATIVE_REGULATOR:
            return "loss_of_function", "possible", "LoF event in tumor suppressor / negative regulator"
        return "loss_of_function", "uncertain", "LoF event; interpret based on gene biology"

    if any(term in text for term in INFRAME_TERMS):
        if gene in ONCOGENE_LIKE:
            return "inframe_indel", "possible", "In-frame event in MAPK-associated oncogene-like gene"
        return "inframe_indel", "uncertain", "In-frame event in MAPK gene"

    if any(term in text for term in MISSENSE_TERMS):
        if gene in ONCOGENE_LIKE:
            return "missense_non_hotspot", "uncertain", "Non-hotspot missense in MAPK-associated oncogene-like gene"
        if gene in TUMOR_SUPPRESSOR_OR_NEGATIVE_REGULATOR:
            return "missense_non_hotspot", "uncertain", "Non-hotspot missense in tumor suppressor / negative regulator"
        return "missense_non_hotspot", "uncertain", "Non-hotspot missense"

    if any(term in text for term in SILENT_TERMS):
        return "silent_or_non_coding", "low", "Synonymous/non-coding annotation"

    return "other_snv_indel", "uncertain", "SNV/InDel event with unclassified consequence"


def write_tsv(path: Path, rows: List[Dict[str, str]], columns: Sequence[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(columns), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> int:
    args = parse_args()

    gene_info = read_gene_set(args.mapk_genes)
    mapk_genes = set(gene_info)
    if not mapk_genes:
        raise ValueError(f"No genes found in {args.mapk_genes}")

    columns, rows = read_tsv(args.variants_tsv)

    gene_col = find_col(columns, GENE_COLUMN_ALIASES)
    if gene_col is None:
        raise ValueError(
            "Could not find a gene-symbol column in variants TSV. "
            f"Supported aliases: {', '.join(GENE_COLUMN_ALIASES)}. "
            f"Observed columns: {', '.join(columns)}"
        )

    sample_col = find_col(columns, SAMPLE_COLUMN_ALIASES)
    consequence_col = find_col(columns, CONSEQUENCE_COLUMN_ALIASES)
    protein_col = find_col(columns, PROTEIN_COLUMN_ALIASES)
    vaf_col = find_col(columns, VAF_COLUMN_ALIASES)
    depth_col = find_col(columns, DEPTH_COLUMN_ALIASES)
    ref_col = find_col(columns, REF_COLUMN_ALIASES)
    alt_col = find_col(columns, ALT_COLUMN_ALIASES)
    chrom_col = find_col(columns, CHROM_COLUMN_ALIASES)
    pos_col = find_col(columns, POS_COLUMN_ALIASES)
    caller_col = find_col(columns, CALLER_COLUMN_ALIASES)

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
            event_class, confidence, note = classify_event(gene, consequence, protein_change)
            info = gene_info.get(gene, {})
            alteration_rows.append({
                "sample": sample,
                "gene": gene,
                "variant_type": variant_type,
                "event_class": event_class,
                "confidence": confidence,
                "consequence": consequence,
                "protein_change": protein_change,
                "vaf": "" if vaf is None else str(vaf),
                "depth": "" if depth is None else str(depth),
                "chrom": row.get(chrom_col, "") if chrom_col else "",
                "pos": row.get(pos_col, "") if pos_col else "",
                "ref": row.get(ref_col, "") if ref_col else "",
                "alt": row.get(alt_col, "") if alt_col else "",
                "caller": row.get(caller_col, "") if caller_col else "",
                "pathway": info.get("pathway", "MAPK signaling pathway"),
                "dataset": info.get("dataset", "custom_MAPK_gene_set"),
                "note": note,
            })

    alteration_columns = [
        "sample", "gene", "variant_type", "event_class", "confidence", "consequence",
        "protein_change", "vaf", "depth", "chrom", "pos", "ref", "alt", "caller",
        "pathway", "dataset", "note"
    ]
    write_tsv(Path(f"{args.out_prefix}.mapk_snvindel.alterations.tsv"), alteration_rows, alteration_columns)

    by_sample: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    by_sample_gene: Dict[Tuple[str, str], List[Dict[str, str]]] = defaultdict(list)
    for row in alteration_rows:
        by_sample[row["sample"]].append(row)
        by_sample_gene[(row["sample"], row["gene"])].append(row)

    # Always emit at least one summary row for the requested sample, even if no MAPK events are found.
    samples = sorted(set(by_sample) | {args.sample_id})
    sample_summary_rows: List[Dict[str, str]] = []
    confidence_rank = {"high": 3, "possible": 2, "uncertain": 1, "low": 0, "": -1}

    for sample in samples:
        events = by_sample.get(sample, [])
        n_high = sum(1 for e in events if e["confidence"] == "high")
        n_possible = sum(1 for e in events if e["confidence"] == "possible")
        n_uncertain = sum(1 for e in events if e["confidence"] == "uncertain")
        n_low = sum(1 for e in events if e["confidence"] == "low")
        genes_hit = sorted({e["gene"] for e in events})

        if n_high > 0:
            status = "High-confidence MAPK SNV/InDel alteration"
        elif n_possible > 0:
            status = "Possible MAPK SNV/InDel alteration"
        elif n_uncertain > 0 or n_low > 0:
            status = "MAPK SNV/InDel variant(s), uncertain significance"
        else:
            status = "No MAPK SNV/InDel alteration detected"

        best = ""
        if events:
            best_event = sorted(
                events,
                key=lambda e: (confidence_rank.get(e["confidence"], -1), e.get("vaf", "")),
                reverse=True,
            )[0]
            best = f"{best_event['gene']} {best_event['protein_change'] or best_event['event_class']}".strip()

        sample_summary_rows.append({
            "sample": sample,
            "n_mapk_snvindel_events": str(len(events)),
            "n_mapk_genes_hit": str(len(genes_hit)),
            "n_high_confidence": str(n_high),
            "n_possible": str(n_possible),
            "n_uncertain": str(n_uncertain),
            "n_low": str(n_low),
            "mapk_snvindel_status": status,
            "top_candidate": best,
            "genes_hit": ",".join(genes_hit),
        })

    sample_summary_columns = [
        "sample", "n_mapk_snvindel_events", "n_mapk_genes_hit", "n_high_confidence",
        "n_possible", "n_uncertain", "n_low", "mapk_snvindel_status", "top_candidate", "genes_hit"
    ]
    write_tsv(Path(f"{args.out_prefix}.mapk_snvindel.sample_summary.tsv"), sample_summary_rows, sample_summary_columns)

    gene_summary_rows: List[Dict[str, str]] = []
    for (sample, gene), events in sorted(by_sample_gene.items()):
        event_counts = Counter(e["event_class"] for e in events)
        best_event = sorted(
            events,
            key=lambda e: confidence_rank.get(e["confidence"], -1),
            reverse=True,
        )[0]
        gene_summary_rows.append({
            "sample": sample,
            "gene": gene,
            "n_events": str(len(events)),
            "event_classes": ",".join(f"{k}:{v}" for k, v in sorted(event_counts.items())),
            "best_event_class": best_event["event_class"],
            "best_confidence": best_event["confidence"],
            "best_protein_change": best_event["protein_change"],
        })

    gene_summary_columns = [
        "sample", "gene", "n_events", "event_classes", "best_event_class", "best_confidence", "best_protein_change"
    ]
    write_tsv(Path(f"{args.out_prefix}.mapk_snvindel.gene_summary.tsv"), gene_summary_rows, gene_summary_columns)

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
