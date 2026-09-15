#!/usr/bin/env python3

import csv
import re
import sys
import argparse

VERSION = "1.16.1"

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

CALLERS = ["FREEBAYES", "LOFREQ_SNV", "LOFREQ_INDELS", "MUTECT2", "VARDICT"]
CLINIC_TN_RATIO_THRESHOLD = 5.0
FREEBAYES_MIN_MQM = 40.0
FREEBAYES_MAX_BIAS_PHRED = 30.0
LOFREQ_MAX_SB = 30.0
LOFREQ_MIN_UQ = 50.0
DMP_HOTSPOT_MIN_VAF = 0.02
DMP_STANDARD_MIN_VAF = 0.05
DMP_FLT3_ITD_MIN_VAF = 0.01

# --- Germline classification thresholds ---
# Mutect2 POPAF is Phred-scaled -log10(population allele frequency); a small
# POPAF means the allele is common in the population (i.e. germline-like).
# We convert POPAF back to a population AF and compare against these floors.
GERMLINE_POPAF_COMMON_AF = 0.01     # pop AF >= 1%   -> strong germline evidence
GERMLINE_POPAF_RARE_AF = 0.001      # pop AF >= 0.1% -> mild germline evidence
# Mutect2 GERMQ is a Phred-scaled probability that the site is NOT germline;
# a low value means germline cannot be confidently ruled out.
GERMLINE_GERMQ_MAX = 30.0
# Mutect2 NLOD is the normal log-odds that the site is hom-ref in the normal;
# a value at/below this suggests the alt allele is present in the normal.
GERMLINE_NLOD_MAX = 0.0
# Matched-normal VAF at/above this is treated as germline presence in the normal.
GERMLINE_NORMAL_VAF_MIN = 0.25
# Evidence-weighted score cutoffs for the 3-state GERMLINE label.
GERMLINE_SCORE_LIKELY = 3
GERMLINE_SCORE_UNLIKELY = 0

GERMLINE_LIKELY = "Likely Germline"
GERMLINE_INCONCLUSIVE = "Non conclusive"
GERMLINE_UNLIKELY = "Likely non-germline"

# VAF-bin overlay (CN-free). --vaf-somatic-max is the somatic-leaning ceiling;
# bins above that are labeled and capped at REVIEW. Not proof of germline.
VAF_SOMATIC_MAX_DEFAULT = 0.39
VAF_HET_MAX = 0.60
VAF_HOM_MIN = 0.90
VAF_BIN_SOMATIC_LEANING = "SOMATIC_LEANING"
VAF_BIN_HET_GERMLINE_LIKE = "HET_GERMLINE_LIKE"
VAF_BIN_ATYPICAL_DIPLOID = "ATYPICAL_DIPLOID"
VAF_BIN_HOM_GERMLINE_LIKE = "HOM_GERMLINE_LIKE"
VAF_BIN_UNKNOWN = "UNKNOWN"
VAF_BINS_CAP_REVIEW = {
    VAF_BIN_HET_GERMLINE_LIKE,
    VAF_BIN_ATYPICAL_DIPLOID,
    VAF_BIN_HOM_GERMLINE_LIKE,
}

HOTSPOT_LIKELY = "Likely Hotspot"
HOTSPOT_INCONCLUSIVE = "Non conclusive"
HOTSPOT_UNLIKELY = "Likely non-hotspot"

# Three-letter to one-letter amino-acid codes, used to normalise VEP/SnpEff
# HGVSp (three-letter) into the one-letter form used by the hotspot list.
AA_THREE_TO_ONE = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*", "Xaa": "X", "Sec": "U", "Pyl": "O",
}

def normalize_output_row(row, fieldnames):
    return {
        field: row.get(field) if row.get(field) not in (None, "") else "."
        for field in fieldnames
    }

def as_int(value):
    try:
        if value in (None, "", "."):
            return None
        return int(round(float(value)))
    except (TypeError, ValueError):
        return None

def as_float(value):
    try:
        if value in (None, "", "."):
            return None
        return float(value)
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

def truthy(value):
    return str(value).lower() in {"1", "true", "yes", "pass", "somatic", "strongsomatic"}

def called_by(row, caller):
    return row.get(f"{caller}_PRESENT") == "1"

def freebayes_gt_is_ref(row):
    gt = row.get("FREEBAYES_FORMAT_GT")
    return gt in {"0/0", "0|0", "0"}

def hard_failed_callers(row):
    failed = {}
    if called_by(row, "FREEBAYES") and freebayes_gt_is_ref(row):
        failed["FREEBAYES"] = "FREEBAYES_GT_0_0"

    mutect_filter = row.get("MUTECT2_FILTER")
    if called_by(row, "MUTECT2") and mutect_filter != "PASS":
        failed["MUTECT2"] = "MUTECT2_FILTER_NOT_PASS"

    vardict_filter = row.get("VARDICT_FILTER")
    if called_by(row, "VARDICT") and vardict_filter != "PASS":
        failed["VARDICT"] = "VARDICT_FILTER_NOT_PASS"

    return failed

def effective_callers(row, failed_callers):
    return [caller for caller in CALLERS if called_by(row, caller) and caller not in failed_callers]

def caller_filter_ok(row, caller):
    flt = row.get(f"{caller}_FILTER")
    status = row.get(f"{caller}_INFO_STATUS")
    somatic = row.get(f"{caller}_INFO_SOMATIC")
    if flt in {"PASS", "."}:
        return True
    if status == "StrongSomatic":
        return True
    if truthy(somatic):
        return True
    return False

def any_filter_ok(row):
    return any(called_by(row, caller) and caller_filter_ok(row, caller) for caller in CALLERS)

def any_effective_filter_ok(row, active_callers):
    return any(caller_filter_ok(row, caller) for caller in active_callers)

def normal_alt_support(row):
    flags = []
    mutect_samples = row.get("MUTECT2_CALLER_SAMPLES", "")
    if mutect_samples:
        normal = mutect_samples.split(";")[0]
        fields = normal.split("=", 1)[-1].split(":")
        fmt = row.get("MUTECT2_CALLER_FORMAT", "").split(":")
        fmt_map = dict(zip(fmt, fields))
        ad = split_ints(fmt_map.get("AD"))
        if len(ad) > 1 and any(x and x > 0 for x in ad[1:]):
            flags.append("MUTECT2_NORMAL_ALT_SUPPORT")
    vardict_samples = row.get("VARDICT_CALLER_SAMPLES", "")
    if vardict_samples:
        parts = vardict_samples.split(";")
        if len(parts) > 1:
            normal = parts[1]
            fields = normal.split("=", 1)[-1].split(":")
            fmt = row.get("VARDICT_CALLER_FORMAT", "").split(":")
            fmt_map = dict(zip(fmt, fields))
            vd = as_int(fmt_map.get("VD"))
            ad = split_ints(fmt_map.get("AD"))
            if (vd and vd > 0) or (len(ad) > 1 and any(x and x > 0 for x in ad[1:])):
                flags.append("VARDICT_MATCH_ALT_SUPPORT")
    return flags

def one_sided_alt_support(row):
    flags = []
    if called_by(row, "FREEBAYES"):
        saf = as_int(row.get("FREEBAYES_INFO_SAF"))
        sar = as_int(row.get("FREEBAYES_INFO_SAR"))
        if saf is not None and sar is not None and saf + sar > 0 and min(saf, sar) == 0:
            flags.append("FREEBAYES_ONE_SIDED_ALT")
    if called_by(row, "LOFREQ_SNV"):
        dp4 = split_ints(row.get("LOFREQ_SNV_INFO_DP4"))
        if len(dp4) == 4:
            alt_forward, alt_reverse = dp4[2], dp4[3]
            if alt_forward is not None and alt_reverse is not None:
                if alt_forward + alt_reverse > 0 and min(alt_forward, alt_reverse) == 0:
                    flags.append("LOFREQ_SNV_ONE_SIDED_ALT")
    if called_by(row, "MUTECT2"):
        table = row.get("MUTECT2_INFO_AS_SB_TABLE", "")
        alleles = table.split("|")
        if len(alleles) > 1:
            for alt_counts in alleles[1:]:
                counts = split_ints(alt_counts)
                if len(counts) == 2 and None not in counts and sum(counts) > 0 and min(counts) == 0:
                    flags.append("MUTECT2_ONE_SIDED_ALT")
                    break
    if called_by(row, "VARDICT"):
        ald = split_ints(row.get("VARDICT_FORMAT_ALD"))
        if len(ald) == 2 and None not in ald and sum(ald) > 0 and min(ald) == 0:
            flags.append("VARDICT_ONE_SIDED_ALT")
    return flags

def caller_artifact_flags(row):
    flags = []

    if called_by(row, "FREEBAYES"):
        mqm = as_float(row.get("FREEBAYES_INFO_MQM"))
        if mqm is not None and mqm < FREEBAYES_MIN_MQM:
            flags.append(f"FREEBAYES_MQM_LT_{int(FREEBAYES_MIN_MQM)}")

        for field in ["SAP", "RPP", "EPP"]:
            value = as_float(row.get(f"FREEBAYES_INFO_{field}"))
            if value is not None and value > FREEBAYES_MAX_BIAS_PHRED:
                flags.append(f"FREEBAYES_{field}_GT_{int(FREEBAYES_MAX_BIAS_PHRED)}")

    if called_by(row, "LOFREQ_SNV"):
        dp4 = split_ints(row.get("LOFREQ_SNV_INFO_DP4"))
        if len(dp4) == 4:
            alt_forward, alt_reverse = dp4[2], dp4[3]
            if alt_forward is not None and alt_reverse is not None:
                alt_total = alt_forward + alt_reverse
                if alt_total > 0:
                    alt_strand_fraction = min(alt_forward, alt_reverse) / alt_total
                    if alt_strand_fraction < 0.1:
                        flags.append("LOFREQ_SNV_ALT_STRAND_FRAC_LT_0.1")
                    elif alt_strand_fraction < 0.2:
                        flags.append("LOFREQ_SNV_ALT_STRAND_FRAC_LT_0.2")

        sb = as_float(row.get("LOFREQ_SNV_INFO_SB"))
        if sb is not None and sb >= LOFREQ_MAX_SB:
            flags.append(f"LOFREQ_SNV_SB_GE_{int(LOFREQ_MAX_SB)}")

        uq = as_float(row.get("LOFREQ_SNV_INFO_UQ"))
        if uq is not None and uq < LOFREQ_MIN_UQ:
            flags.append(f"LOFREQ_SNV_UQ_LT_{int(LOFREQ_MIN_UQ)}")

    return flags

def selected_alt_reads(row):
    vd = as_int(row.get("VD"))
    if vd is not None:
        return vd
    values = []
    for caller in CALLERS:
        if called_by(row, caller):
            caller_vd = as_int(row.get(f"{caller}_VD"))
            if caller_vd is not None:
                values.append(caller_vd)
    return min(values) if values else None

def selected_depth(row):
    dp = as_int(row.get("DP"))
    if dp is not None:
        return dp
    values = []
    for caller in CALLERS:
        if called_by(row, caller):
            caller_dp = as_int(row.get(f"{caller}_DP"))
            if caller_dp is not None:
                values.append(caller_dp)
    return min(values) if values else None

VAF_CALLERS_NO_FB = ["MUTECT2", "LOFREQ_SNV", "LOFREQ_INDELS"]

def selected_vaf(row):
    af = as_float(row.get("AF"))
    if af is not None:
        return af
    if called_by(row, "VARDICT"):
        vardict_af = as_float(row.get("VARDICT_AF"))
        if vardict_af is not None:
            return vardict_af
    values = []
    for caller in VAF_CALLERS_NO_FB:
        if called_by(row, caller):
            caller_af = as_float(row.get(f"{caller}_AF"))
            if caller_af is not None:
                values.append(caller_af)
    return sum(values) / len(values) if values else None

def is_indel(row):
    return len(row.get("REF", "")) != len(row.get("ALT", ""))

def gene_symbol(row):
    gene = row.get("VEP_SYMBOL") or row.get("SNPEFF_GENE_NAME") or ""
    return gene.strip("()").upper()

def first_token(value):
    if value in (None, "", "."):
        return None
    return str(value).split(",")[0]

def normalize_hgvsp(value):
    """Normalise a VEP/SnpEff HGVSp string to the one-letter form used by the
    hotspot list (e.g. 'ENSP..:p.Cys746Tyr' -> 'p.C746Y', 'p.Cys746Cys' ->
    'p.C746='). Returns None for anything that is not a protein change."""
    if value in (None, "", "."):
        return None
    text = str(value).strip()
    if ":" in text:
        text = text.split(":")[-1]
    text = text.replace("%3D", "=")
    if not text.startswith("p."):
        return None
    for three, one in AA_THREE_TO_ONE.items():
        text = text.replace(three, one)
    synonymous = re.fullmatch(r"p\.([A-Z*])(\d+)\1", text)
    if synonymous:
        text = f"p.{synonymous.group(1)}{synonymous.group(2)}="
    return text

def normalize_chrom(value):
    if value in (None, "", "."):
        return None
    chrom = str(value).strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    return chrom.upper()

def minimal_variant_key(chrom, pos, ref, alt):
    """Reference-free minimal representation of a variant so that MAF-style
    (trimmed / '-') and VCF-style (anchor-base) alleles compare equal. Trims the
    shared suffix then the shared prefix, adjusting POS. Returns None if
    unusable. Note: exact left-alignment across repeats needs the reference
    FASTA, which is unavailable here, so repeat-region indels may not reconcile;
    protein-change matching backstops those coding cases."""
    chrom = normalize_chrom(chrom)
    pos = as_int(pos)
    if chrom is None or pos is None or ref is None or alt is None:
        return None
    ref = "" if ref == "-" else str(ref).upper()
    alt = "" if alt == "-" else str(alt).upper()
    while ref and alt and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while ref and alt and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return (chrom, pos, ref, alt)

class HotspotIndex:
    def __init__(self):
        self.evidence = {}   # (gene, normalized HGVSp) -> set(ids)  [protein]
        self.coords = {}     # minimal (chrom,pos,ref,alt) -> set(ids)  [coordinate]
        self.genes = set()

    def add_protein(self, gene, hgvsp, existing):
        self.genes.add(gene)
        key = (gene, hgvsp)
        self.evidence.setdefault(key, set())
        if existing and existing not in ("", "."):
            self.evidence[key].add(existing)

    def add_coord(self, key, existing):
        self.coords.setdefault(key, set())
        if existing and existing not in ("", "."):
            self.coords[key].add(existing)

    def __contains__(self, key):
        return key in self.evidence

    def ids(self, key):
        return self.evidence.get(key, set())

    def coord_ids(self, key):
        return self.coords.get(key, set())

def _hotspot_coord_key(row):
    """Minimal-representation key for a hotspot MAF row, reconciling MAF's '-'
    and anchor-free indel positions with VCF-style variants."""
    chrom = row.get("Chromosome")
    start = as_int(row.get("Start_Position"))
    ref = row.get("Reference_Allele")
    alt = row.get("Tumor_Seq_Allele2")
    if start is None:
        return None
    if alt == "-":                       # deletion
        return minimal_variant_key(chrom, start, ref, "")
    if ref == "-":                       # insertion (MAF Start = base before)
        return minimal_variant_key(chrom, start + 1, "", alt)
    return minimal_variant_key(chrom, start, ref, alt)

def load_hotspots(path):
    """Load the hg38 hotspot list into both a protein-change index
    (gene + normalized HGVSp) and a coordinate index (minimal chrom/pos/ref/alt).
    Matching is the union of the two (additive). Every row is loaded."""
    index = HotspotIndex()
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"Hugo_Symbol", "Chromosome", "Start_Position",
                    "Reference_Allele", "Tumor_Seq_Allele2", "HGVSp_Short"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(
                "Hotspot file must be tab-delimited with columns: "
                + ", ".join(sorted(required))
            )
        for row in reader:
            existing = (row.get("Existing_variation") or "").strip()
            gene = (row.get("Hugo_Symbol") or "").strip("()").upper()
            hgvsp = normalize_hgvsp(row.get("HGVSp_Short"))
            if gene and hgvsp:
                index.add_protein(gene, hgvsp, existing)
            coord_key = _hotspot_coord_key(row)
            if coord_key is not None:
                index.add_coord(coord_key, existing)
    return index

def split_variant_ids(id_strings):
    """Flatten a set of (possibly comma-joined) Existing_variation strings from
    the hotspot list into a set of individual variant IDs (rs*, COSM*, etc.)."""
    out = set()
    for s in id_strings:
        for token in re.split(r"[,&;|]", s or ""):
            token = token.strip()
            if token and token not in (".",):
                out.add(token)
    return out

def hotspot_candidates(row):
    """(annotator, gene, normalized HGVSp) pairs available for this row."""
    candidates = []
    vep_gene = (row.get("VEP_SYMBOL") or "").strip("()").upper()
    vep_p = normalize_hgvsp(row.get("VEP_HGVSP"))
    if vep_gene and vep_p:
        candidates.append(("VEP", vep_gene, vep_p))
    eff_gene = (row.get("SNPEFF_GENE_NAME") or "").strip("()").upper()
    eff_p = normalize_hgvsp(row.get("SNPEFF_HGVS_P"))
    if eff_gene and eff_p:
        candidates.append(("SNPEFF", eff_gene, eff_p))
    return candidates

def hotspot_evaluate(row, hotspots):
    """Evaluate hotspot status via protein-change AND coordinate matching (union,
    additive). Returns a dict with the label/info plus booleans used by tiering."""
    candidates = hotspot_candidates(row)
    matched_ids = set()
    protein_details = []
    for source, gene, hgvsp in candidates:
        key = (gene, hgvsp)
        if key in hotspots:
            matched_ids |= split_variant_ids(hotspots.ids(key))
            ids = ",".join(sorted(hotspots.ids(key)))
            protein_details.append(f"{source}:{gene}:{hgvsp}" + (f"({ids})" if ids else ""))

    coord_key = minimal_variant_key(
        row.get("CHROM"), row.get("POS"), row.get("REF"), row.get("ALT")
    )
    coord_detail = None
    if coord_key is not None and coord_key in hotspots.coords:
        matched_ids |= split_variant_ids(hotspots.coord_ids(coord_key))
        ids = ",".join(sorted(hotspots.coord_ids(coord_key)))
        c, p, r, a = coord_key
        coord_detail = f"POS:{c}:{p}:{r or '-'}:{a or '-'}" + (f"({ids})" if ids else "")

    ids_str = ",".join(sorted(matched_ids))
    protein_match = bool(protein_details)
    coord_match = coord_detail is not None
    list_match = protein_match or coord_match
    # Coordinate matched but the annotators' protein change did not -> the
    # variant is the same locus on a different transcript than the list uses.
    transcript_discrepancy = coord_match and bool(candidates) and not protein_match

    if list_match:
        details = list(protein_details)
        if coord_detail:
            details.append(coord_detail)
        if transcript_discrepancy:
            details.append("TRANSCRIPT_DISCREPANCY")
        return {
            "label": HOTSPOT_LIKELY, "info": "; ".join(details), "ids": ids_str,
            "list_match": True, "coord_match": coord_match,
            "protein_match": protein_match, "transcript_discrepancy": transcript_discrepancy,
        }

    exception_flags = dmp_exception_flags(row)
    if exception_flags:
        return {
            "label": HOTSPOT_LIKELY, "info": "DMP_EXCEPTION:" + ",".join(exception_flags),
            "ids": "", "list_match": False, "coord_match": False,
            "protein_match": False, "transcript_discrepancy": False,
        }
    if candidates:
        shown = "; ".join(f"{s}:{g}:{p}" for s, g, p in candidates)
        info = f"evaluated, not in hotspot list ({shown})"
        label = HOTSPOT_UNLIKELY
    else:
        info = "no coding protein change available to evaluate"
        label = HOTSPOT_UNLIKELY if coord_key is not None else HOTSPOT_INCONCLUSIVE
    return {
        "label": label, "info": info, "ids": "", "list_match": False,
        "coord_match": False, "protein_match": False, "transcript_discrepancy": False,
    }

def germline_classify(row):
    """Evidence-weighted 3-state germline call from Mutect2 (POPAF/GERMQ/NLOD/
    PON/FILTER), VarDict STATUS, and matched-normal VAF."""
    score = 0
    evidence = []
    mfilter_tags = set((row.get("MUTECT2_FILTER") or "").split(";"))

    if called_by(row, "MUTECT2") and "germline" in mfilter_tags:
        score += 3
        evidence.append("MUTECT2_FILTER=germline")

    popaf = as_float(first_token(row.get("MUTECT2_INFO_POPAF")))
    if popaf is not None:
        pop_af = 10 ** (-popaf)
        if pop_af >= GERMLINE_POPAF_COMMON_AF:
            score += 3
            evidence.append(f"POPAF={popaf:g}(pop_af~{pop_af:.3g})")
        elif pop_af >= GERMLINE_POPAF_RARE_AF:
            score += 1
            evidence.append(f"POPAF={popaf:g}(pop_af~{pop_af:.3g})")
        else:
            score -= 1
            evidence.append(f"POPAF={popaf:g}(rare,somatic-leaning)")

    germq = as_float(first_token(row.get("MUTECT2_INFO_GERMQ")))
    if germq is not None:
        if germq < GERMLINE_GERMQ_MAX:
            score += 1
            evidence.append(f"GERMQ={germq:g}(<{GERMLINE_GERMQ_MAX:g})")
        else:
            score -= 1

    nlod = as_float(first_token(row.get("MUTECT2_INFO_NLOD")))
    if nlod is not None and nlod <= GERMLINE_NLOD_MAX:
        score += 1
        evidence.append(f"NLOD={nlod:g}(alt present in normal)")

    if row.get("MUTECT2_INFO_PON") not in (None, "", "."):
        score += 1
        evidence.append("MUTECT2_PON")
    if "panel_of_normals" in mfilter_tags:
        score += 1
        evidence.append("MUTECT2_FILTER=panel_of_normals")

    status = row.get("VARDICT_INFO_STATUS") or ""
    if status == "Germline":
        score += 3
        evidence.append("VARDICT_STATUS=Germline")
    elif status in {"StrongLOH", "LikelyLOH"}:
        score += 1
        evidence.append(f"VARDICT_STATUS={status}")
    elif status in {"StrongSomatic", "LikelySomatic"}:
        score -= 2
        evidence.append(f"VARDICT_STATUS={status}")

    normal_vaf = matched_normal_vaf(row)
    if normal_vaf >= GERMLINE_NORMAL_VAF_MIN:
        score += 3
        evidence.append(f"NORMAL_VAF={normal_vaf:.3g}(>={GERMLINE_NORMAL_VAF_MIN:g})")
    elif normal_vaf > 0:
        evidence.append(f"NORMAL_VAF={normal_vaf:.3g}")

    if row.get("LOFREQ_SNV_INFO_SOMATIC") not in (None, "", "."):
        score -= 1
        evidence.append("LOFREQ_SOMATIC")

    if not evidence:
        return GERMLINE_INCONCLUSIVE, "no germline/somatic signals available"
    if score >= GERMLINE_SCORE_LIKELY:
        label = GERMLINE_LIKELY
    elif score <= GERMLINE_SCORE_UNLIKELY:
        label = GERMLINE_UNLIKELY
    else:
        label = GERMLINE_INCONCLUSIVE
    return label, "; ".join(evidence)

def vaf_bin_label(vaf, vaf_somatic_max):
    """CN-free VAF bin. Soft labels only; does not drop variants."""
    if vaf is None:
        return VAF_BIN_UNKNOWN
    if vaf <= vaf_somatic_max:
        return VAF_BIN_SOMATIC_LEANING
    if vaf <= VAF_HET_MAX:
        return VAF_BIN_HET_GERMLINE_LIKE
    if vaf < VAF_HOM_MIN:
        return VAF_BIN_ATYPICAL_DIPLOID
    return VAF_BIN_HOM_GERMLINE_LIKE

def vaf_germline_overlay(vaf, vaf_somatic_max):
    """Append-only VAF tags for GERMLINE_INFO. Suspicious bins cap HC/SINGLE_CALLER at REVIEW."""
    vaf_bin = vaf_bin_label(vaf, vaf_somatic_max)
    tags = [f"VAF_BIN={vaf_bin}"]
    if vaf is not None:
        tags.append(f"TUMOR_VAF={vaf:.3g}")
    if vaf_bin == VAF_BIN_ATYPICAL_DIPLOID:
        tags.append("VAF_CONTEXT=REQUIRES_CN_FOR_MECHANISM")
    elif vaf_bin == VAF_BIN_HOM_GERMLINE_LIKE:
        tags.append("VAF_CONTEXT=REQUIRES_CN_FOR_MECHANISM")
    return vaf_bin, "; ".join(tags)

def row_pos(row):
    return as_int(row.get("POS"))

def row_end(row):
    pos = row_pos(row)
    if pos is None:
        return None
    ref = row.get("REF", "")
    return pos + max(len(ref), 1) - 1

def in_interval(value, start, end):
    return value is not None and start < value < end

def overlaps_interval(row, start, end):
    pos = row_pos(row)
    end_pos = row_end(row)
    return (
        pos is not None and end_pos is not None
        and pos <= end
        and end_pos >= start
    )

def annotation_text(row):
    fields = [
        "VEP_CONSEQUENCE", "SNPEFF_ANNOTATION", "VEP_IMPACT", "SNPEFF_ANNOTATION_IMPACT",
        "VEP_HGVSP", "SNPEFF_HGVS_P", "VEP_AMINO_ACIDS",
    ]
    return "|".join(row.get(field, "") for field in fields).lower()

def is_insertion(row):
    return len(row.get("REF", "")) < len(row.get("ALT", ""))

def has_any_annotation(row, terms):
    text = annotation_text(row)
    return any(term in text for term in terms)

def dmp_exception_flags(row):
    flags = []
    gene = gene_symbol(row)
    pos = row_pos(row)
    text = annotation_text(row)

    # GRCh38/hg38 coordinates lifted from the legacy DMP hg19 exception windows.
    if gene == "MET" and (
        overlaps_interval(row, 116771799, 116771848)
        or overlaps_interval(row, 116771990, 116772039)
    ):
        flags.append("DMP_EXCEPTION_MET_INTRON_13_14_SPLICE_REGION")

    # GRCh38/hg38 positions for the canonical TERT promoter hotspots
    # commonly referred to as C228T and C250T in older literature.
    if gene == "TERT" and pos in {1295113, 1295135}:
        flags.append("DMP_EXCEPTION_TERT_PROMOTER_HOTSPOT_GRCH38")

    if gene == "FLT3" and is_insertion(row) and in_interval(pos, 28028288, 28034301):
        flags.append("DMP_EXCEPTION_FLT3_ITD_VAF_0.01")

    if gene == "CALR" and in_interval(pos, 12943712, 12943830) and has_any_annotation(row, ["frameshift", "frame_shift"]):
        flags.append("DMP_EXCEPTION_CALR_EXON9_FRAMESHIFT")

    # This hg19 interval lifts to the same GRCh38 numeric coordinates.
    if gene == "JAK2" and in_interval(pos, 5070016, 5070053):
        if is_indel(row) and has_any_annotation(row, ["inframe", "in_frame"]):
            flags.append("DMP_EXCEPTION_JAK2_EXON12_INFRAME_INDEL")
        elif has_any_annotation(row, ["missense"]):
            flags.append("DMP_EXCEPTION_JAK2_EXON12_MISSENSE")

    if gene == "CEBPA" and in_interval(pos, 33301379, 33301584):
        if is_indel(row) and has_any_annotation(row, ["inframe", "in_frame"]):
            flags.append("DMP_EXCEPTION_CEBPA_INFRAME_INDEL")
        elif has_any_annotation(row, ["missense"]):
            flags.append("DMP_EXCEPTION_CEBPA_MISSENSE")

    if gene == "FLT3" and is_insertion(row) and has_any_annotation(row, ["inframe", "in_frame"]):
        flags.append("DMP_EXCEPTION_FLT3_INFRAME_INSERTION")

    if gene == "CSF3R" and in_interval(pos, 36466359, 36467566):
        if has_any_annotation(row, ["frameshift", "frame_shift"]):
            flags.append("DMP_EXCEPTION_CSF3R_FRAMESHIFT")
        elif has_any_annotation(row, ["nonsense", "stop_gained"]):
            flags.append("DMP_EXCEPTION_CSF3R_NONSENSE")

    return flags

def is_dmp_exception(row):
    return bool(dmp_exception_flags(row))

def dmp_min_vaf(row, hotspot_like):
    if "DMP_EXCEPTION_FLT3_ITD_VAF_0.01" in dmp_exception_flags(row):
        return DMP_FLT3_ITD_MIN_VAF
    return DMP_HOTSPOT_MIN_VAF if hotspot_like else DMP_STANDARD_MIN_VAF

def sample_field_map(sample_payload, sample_index, format_value):
    if not sample_payload or sample_payload == ".":
        return {}
    parts = sample_payload.split(";")
    if sample_index >= len(parts):
        return {}
    sample_values = parts[sample_index].split("=", 1)[-1].split(":")
    fmt = format_value.split(":") if format_value and format_value != "." else []
    return dict(zip(fmt, sample_values))

def matched_normal_vaf(row):
    values = []

    mutect_map = sample_field_map(
        row.get("MUTECT2_CALLER_SAMPLES", ""),
        0,
        row.get("MUTECT2_CALLER_FORMAT", ""),
    )
    mutect_af = mutect_map.get("AF")
    if mutect_af:
        values.extend(x for x in mutect_af.split(",") if x not in {"", "."})

    vardict_map = sample_field_map(
        row.get("VARDICT_CALLER_SAMPLES", ""),
        1,
        row.get("VARDICT_CALLER_FORMAT", ""),
    )
    vardict_af = vardict_map.get("AF")
    if vardict_af:
        values.extend(x for x in vardict_af.split(",") if x not in {"", "."})

    floats = [as_float(x) for x in values]
    floats = [x for x in floats if x is not None]
    return max(floats) if floats else 0.0

def clinic_filter_flags(row, dp, alt_reads, vaf, hotspot_list_match=False, apply_clinic_floors=False):
    """Clinic support flags. Default (cfDNA): no VAF/AD/DP floors — only missing-data
    and T/N ratio. Pass apply_clinic_floors=True to enable tissue-style thresholds."""
    flags = []
    rescue = []
    exception_flags = dmp_exception_flags(row)
    hotspot_like = is_high_impact(row) or bool(exception_flags) or hotspot_list_match
    normal_vaf = matched_normal_vaf(row)
    tn_ratio = CLINIC_TN_RATIO_THRESHOLD if normal_vaf == 0 else (vaf or 0) / normal_vaf
    variant_kind = "INDEL" if is_indel(row) else "SNV"

    if apply_clinic_floors:
        min_dp = 20
        min_ad = 8 if hotspot_like else 10
        min_vaf = dmp_min_vaf(row, hotspot_like)
        if dp is None:
            flags.append("CLINIC_NO_DP")
        elif dp < min_dp:
            flags.append(f"CLINIC_{variant_kind}_DP_LT_{min_dp}")
        if alt_reads is None:
            flags.append("CLINIC_NO_AD")
        elif alt_reads < min_ad:
            flags.append(f"CLINIC_{variant_kind}_AD_LT_{min_ad}")
        if vaf is None:
            flags.append("CLINIC_NO_VAF")
        elif vaf < min_vaf:
            flags.append(f"CLINIC_{variant_kind}_VAF_LT_{min_vaf:g}")
    else:
        if dp is None:
            flags.append("CLINIC_NO_DP")
        if alt_reads is None:
            flags.append("CLINIC_NO_AD")
        if vaf is None:
            flags.append("CLINIC_NO_VAF")

    if tn_ratio < CLINIC_TN_RATIO_THRESHOLD:
        flags.append(f"CLINIC_TN_RATIO_LT_{CLINIC_TN_RATIO_THRESHOLD:g}")

    if hotspot_like:
        if is_high_impact(row):
            rescue.append("CLINIC_HOTSPOT_LIKE_HIGH_IMPACT")
        if hotspot_list_match:
            rescue.append("CLINIC_HOTSPOT_LIST_MATCH")
        rescue.extend(exception_flags)
    if not any(flag.startswith("CLINIC_") and "_LT_" in flag for flag in flags):
        rescue.append("CLINIC_THRESHOLD_PASS")

    return flags, rescue

def is_high_impact(row):
    consequence = "|".join([
        row.get("VEP_CONSEQUENCE", ""),
        row.get("SNPEFF_ANNOTATION", ""),
    ]).lower()
    impact = "|".join([
        row.get("VEP_IMPACT", ""),
        row.get("SNPEFF_ANNOTATION_IMPACT", ""),
    ]).upper()
    high_terms = [
        "stop_gained", "stop_lost", "start_lost", "frameshift",
        "splice_acceptor", "splice_donor", "transcript_ablation",
        "inframe_insertion", "inframe_deletion",
    ]
    return impact == "HIGH" or any(term in consequence for term in high_terms)

def is_clinic_rescue_annotation(row):
    return is_high_impact(row) or is_dmp_exception(row)

def is_high_or_moderate_impact(row):
    impact = "|".join([
        row.get("VEP_IMPACT", ""),
        row.get("SNPEFF_ANNOTATION_IMPACT", ""),
    ]).upper()
    return "HIGH" in impact.split("|") or "MODERATE" in impact.split("|")

def has_priority_coding_annotation(row):
    consequence = "|".join([
        row.get("VEP_CONSEQUENCE", ""),
        row.get("SNPEFF_ANNOTATION", ""),
    ]).lower()
    priority_terms = [
        "missense", "protein_altering", "frameshift", "stop_gained", "stop_lost",
        "start_lost", "inframe", "splice",
    ]
    return any(term in consequence for term in priority_terms)

def eligible_for_high_confidence(row):
    return (
        has_priority_coding_annotation(row)
        or is_high_or_moderate_impact(row)
        or is_dmp_exception(row)
    )

def low_impact_or_noncoding_flags(row):
    if is_clinic_rescue_annotation(row):
        return []

    consequence = "|".join([
        row.get("VEP_CONSEQUENCE", ""),
        row.get("SNPEFF_ANNOTATION", ""),
    ]).lower()
    impact = "|".join([
        row.get("VEP_IMPACT", ""),
        row.get("SNPEFF_ANNOTATION_IMPACT", ""),
    ]).upper()
    biotype = "|".join([
        row.get("VEP_BIOTYPE", ""),
        row.get("SNPEFF_TRANSCRIPT_BIOTYPE", ""),
    ]).lower()

    coding_terms = [
        "missense", "protein_altering", "frameshift", "stop_gained", "stop_lost",
        "start_lost", "inframe", "splice_acceptor", "splice_donor",
    ]
    if any(term in consequence for term in coding_terms):
        return []

    noncoding_terms = [
        "synonymous", "intron", "intergenic", "upstream", "downstream", "utr",
        "non_coding", "ncrna", "regulatory", "tf_binding", "mature_mirna",
    ]
    flags = []
    if any(term in consequence for term in noncoding_terms):
        flags.append("LOW_IMPACT_OR_NONCODING_ANNOTATION")
    elif "MODIFIER" in impact or impact == "LOW" or "|LOW" in impact or "LOW|" in impact:
        flags.append("LOW_IMPACT_OR_NONCODING_ANNOTATION")

    if biotype and "protein_coding" not in biotype and biotype != ".|.":
        flags.append("NON_PROTEIN_CODING_TRANSCRIPT")

    return sorted(set(flags))

def tier_row(row, min_alt_reads_review, min_alt_reads_strong, hotspots, apply_clinic_floors=False, vaf_somatic_max=VAF_SOMATIC_MAX_DEFAULT):
    flags = []
    rescue = []
    failed_callers = hard_failed_callers(row)
    active_callers = effective_callers(row, failed_callers)
    num_callers = len(active_callers)
    alt_reads = selected_alt_reads(row)
    dp = selected_depth(row)
    vaf = selected_vaf(row)
    filter_ok = any_effective_filter_ok(row, active_callers)
    hotspot = hotspot_evaluate(row, hotspots)
    hotspot_list_match = hotspot["list_match"]
    hotspot_label = hotspot["label"]
    hotspot_info = hotspot["info"]
    hotspot_ids = hotspot["ids"]
    germline_label, germline_info = germline_classify(row)
    is_likely_germline = germline_label == GERMLINE_LIKELY
    vaf_bin, vaf_overlay = vaf_germline_overlay(vaf, vaf_somatic_max)
    if germline_info:
        germline_info = f"{germline_info}; {vaf_overlay}"
    else:
        germline_info = vaf_overlay
    # A hotspot-list hit is treated as at least as protective as the hardcoded
    # exceptions/high-impact annotation (additive, never replacing them).
    clinic_rescue_anno = is_clinic_rescue_annotation(row) or hotspot_list_match

    flags.extend(failed_callers.values())
    if num_callers <= 1:
        flags.append("SINGLE_CALLER")
    if alt_reads is None:
        flags.append("NO_ALT_READ_SUMMARY")
    elif alt_reads <= 1:
        flags.append("ALT_READS_LE_1")
    elif alt_reads < min_alt_reads_review:
        flags.append(f"ALT_READS_LT_{min_alt_reads_review}")
    if not filter_ok:
        flags.append("NO_PASS_OR_STRONG_CALLER_STATUS")
    flags.extend(normal_alt_support(row))
    flags.extend(one_sided_alt_support(row))
    flags.extend(caller_artifact_flags(row))
    flags.extend(low_impact_or_noncoding_flags(row))
    clinic_flags, clinic_rescue = clinic_filter_flags(
        row, dp, alt_reads, vaf, hotspot_list_match, apply_clinic_floors
    )
    flags.extend(clinic_flags)
    rescue.extend(clinic_rescue)
    if num_callers >= 2:
        rescue.append("MULTI_CALLER")
    if is_high_impact(row):
        rescue.append("HIGH_IMPACT_ANNOTATION")
    if is_dmp_exception(row):
        rescue.append("DMP_HARD_CODED_EXCEPTION")
    if hotspot_list_match:
        rescue.append("HOTSPOT_LIST_MATCH")
    if hotspot["coord_match"]:
        rescue.append("HOTSPOT_COORD_MATCH")
    if hotspot["transcript_discrepancy"]:
        rescue.append("HOTSPOT_TRANSCRIPT_DISCREPANCY")
    if is_likely_germline:
        flags.append("LIKELY_GERMLINE")
    if hotspot_list_match and is_likely_germline:
        rescue.append("HOTSPOT_GERMLINE_CONFLICT")
    if filter_ok:
        rescue.append("CALLER_STATUS_OK")
    if alt_reads is not None and alt_reads >= min_alt_reads_strong:
        rescue.append(f"ALT_READS_GE_{min_alt_reads_strong}")
    hard_filter = bool(failed_callers) and not active_callers
    clinic_tn_fail = any(flag.startswith("CLINIC_TN_RATIO_LT_") for flag in flags)
    clinic_depth_ad_fail = apply_clinic_floors and any(
        flag.startswith("CLINIC_") and ("_DP_LT_" in flag or "_AD_LT_" in flag or "_VAF_LT_" in flag)
        for flag in flags
    )
    normal_flag = any(flag.endswith("_ALT_SUPPORT") for flag in flags)
    one_sided = any(flag.endswith("_ONE_SIDED_ALT") for flag in flags)
    low_impact_or_noncoding = any(
        flag in {"LOW_IMPACT_OR_NONCODING_ANNOTATION", "NON_PROTEIN_CODING_TRANSCRIPT"}
        for flag in flags
    )
    very_low_alt = alt_reads is not None and alt_reads <= 1
    low_alt = alt_reads is not None and alt_reads < min_alt_reads_review
    high_conf_alt_ok = alt_reads is None or alt_reads >= min_alt_reads_review
    hotspot_list_hc = (
        hotspot_list_match
        and num_callers >= 1
        and filter_ok
        and not normal_flag
        and not clinic_tn_fail
        and alt_reads is not None
        and alt_reads >= 1
    )
    standard_hc = (
        (eligible_for_high_confidence(row) or hotspot_list_match)
        and num_callers >= 2
        and not normal_flag
        and not clinic_tn_fail
        and filter_ok
        and not one_sided
        and high_conf_alt_ok
    )
    # Same bar as standard_hc, but only one caller → dedicated tier (not HC).
    standard_hc_single = (
        num_callers == 1
        and (eligible_for_high_confidence(row) or hotspot_list_match)
        and not normal_flag
        and not clinic_tn_fail
        and filter_ok
        and not one_sided
        and high_conf_alt_ok
    )
    if hard_filter:
        tier = "LIKELY_NOISE"
    elif standard_hc or hotspot_list_hc:
        tier = "HIGH_CONFIDENCE"
        if hotspot_list_hc and not standard_hc:
            rescue.append("HOTSPOT_LIST_SINGLE_CALLER_HC")
    elif standard_hc_single:
        tier = "SINGLE_CALLER"
    elif clinic_tn_fail and num_callers <= 1:
        tier = "LIKELY_NOISE"
    elif clinic_depth_ad_fail and num_callers <= 1 and not clinic_rescue_anno:
        tier = "LIKELY_NOISE"
    elif normal_flag and num_callers <= 1:
        tier = "LIKELY_NOISE"
    elif low_impact_or_noncoding and num_callers <= 1 and not clinic_rescue_anno:
        tier = "LIKELY_NOISE"
    elif num_callers <= 1 and very_low_alt and not clinic_rescue_anno:
        tier = "LIKELY_NOISE"
    elif num_callers <= 1 and low_alt and (not filter_ok or one_sided) and not clinic_rescue_anno:
        tier = "LIKELY_NOISE"
    else:
        tier = "REVIEW"
    # Likely-germline caps HC/SINGLE_CALLER at REVIEW; never pushes to LIKELY_NOISE.
    if is_likely_germline and tier in {"HIGH_CONFIDENCE", "SINGLE_CALLER"}:
        tier = "REVIEW"
        flags.append("GERMLINE_DEMOTED_FROM_HIGH_CONFIDENCE")
    # VAF-bin overlay: het/atypical/hom-like VAFs cannot be HC or SINGLE_CALLER.
    # Existing germline evidence still wins; this only adds a REVIEW cap.
    if vaf_bin in VAF_BINS_CAP_REVIEW and tier in {"HIGH_CONFIDENCE", "SINGLE_CALLER"}:
        tier = "REVIEW"
        flags.append(f"VAF_BIN_CAP_REVIEW={vaf_bin}")
    score = 0
    score += 2 if num_callers <= 1 else -2
    score += 3 if very_low_alt else 1 if low_alt else -1
    score += 3 if normal_flag else 0
    score += 1 if one_sided else 0
    score += 2 if not filter_ok else -1
    score += 4 if hard_filter else 0
    score += 2 if clinic_tn_fail else 0
    score += 1 if clinic_depth_ad_fail else 0
    score += 1 if low_impact_or_noncoding else 0
    score -= 2 if is_high_impact(row) else 0
    score -= 2 if is_dmp_exception(row) else 0
    score -= 2 if hotspot_list_match else 0
    score += 2 if is_likely_germline else 0
    score += 1 if vaf_bin in VAF_BINS_CAP_REVIEW else 0
    return {
        "EFFECTIVE_NUM_CALLERS": str(num_callers),
        "EFFECTIVE_CALLERS": ",".join(active_callers) if active_callers else "NONE",
        "HARD_FILTER_FLAGS": ",".join(failed_callers.values()) if failed_callers else "NONE",
        "NOISE_TIER": tier,
        "NOISE_SCORE": str(score),
        "NOISE_FLAGS": ",".join(flags) if flags else "NONE",
        "RESCUE_FLAGS": ",".join(rescue) if rescue else "NONE",
        "KEEP_FOR_REVIEW": "0" if tier == "LIKELY_NOISE" else "1",
        "GERMLINE": germline_label,
        "GERMLINE_INFO": germline_info,
        "HOTSPOT": hotspot_label,
        "HOTSPOT_INFO": hotspot_info,
        "HOTSPOT_IDS": hotspot_ids if hotspot_ids else ".",
    }

def write_filtered(
    input_path,
    output_path,
    high_confidence_output,
    review_output,
    single_caller_output,
    noise_output,
    min_alt_reads_review,
    min_alt_reads_strong,
    hotspots,
    apply_clinic_floors=False,
    vaf_somatic_max=VAF_SOMATIC_MAX_DEFAULT,
):
    with open(input_path, newline="") as in_f, open(output_path, "w", newline="") as out_f:
        reader = csv.DictReader(in_f, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Input TSV has no header")
        added = [
            "EFFECTIVE_NUM_CALLERS", "EFFECTIVE_CALLERS", "HARD_FILTER_FLAGS",
            "NOISE_TIER", "NOISE_SCORE", "NOISE_FLAGS", "RESCUE_FLAGS", "KEEP_FOR_REVIEW",
            "GERMLINE", "GERMLINE_INFO", "HOTSPOT", "HOTSPOT_INFO", "HOTSPOT_IDS",
        ]
        fieldnames = reader.fieldnames + [x for x in added if x not in reader.fieldnames]
        writer = csv.DictWriter(out_f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        high_f = open(high_confidence_output, "w", newline="") if high_confidence_output else None
        review_f = open(review_output, "w", newline="") if review_output else None
        single_f = open(single_caller_output, "w", newline="") if single_caller_output else None
        noise_f = open(noise_output, "w", newline="") if noise_output else None
        high_writer = review_writer = single_writer = noise_writer = None
        if high_f:
            high_writer = csv.DictWriter(high_f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
            high_writer.writeheader()
        if review_f:
            review_writer = csv.DictWriter(review_f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
            review_writer.writeheader()
        if single_f:
            single_writer = csv.DictWriter(single_f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
            single_writer.writeheader()
        if noise_f:
            noise_writer = csv.DictWriter(noise_f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
            noise_writer.writeheader()
        counts = {"HIGH_CONFIDENCE": 0, "SINGLE_CALLER": 0, "REVIEW": 0, "LIKELY_NOISE": 0}
        try:
            for row in reader:
                row.update(tier_row(
                    row, min_alt_reads_review, min_alt_reads_strong, hotspots,
                    apply_clinic_floors, vaf_somatic_max,
                ))
                counts[row["NOISE_TIER"]] = counts.get(row["NOISE_TIER"], 0) + 1
                out_row = normalize_output_row(row, fieldnames)
                writer.writerow(out_row)
                if high_writer and row["NOISE_TIER"] == "HIGH_CONFIDENCE":
                    high_writer.writerow(out_row)
                if single_writer and row["NOISE_TIER"] == "SINGLE_CALLER":
                    single_writer.writerow(out_row)
                if review_writer and row["NOISE_TIER"] == "REVIEW":
                    review_writer.writerow(out_row)
                if noise_writer and row["NOISE_TIER"] == "LIKELY_NOISE":
                    noise_writer.writerow(out_row)
        finally:
            if high_f:
                high_f.close()
            if review_f:
                review_f.close()
            if single_f:
                single_f.close()
            if noise_f:
                noise_f.close()
    return counts

def main():
    parser = argparse.ArgumentParser(
        description="Tier merged cfDNA caller TSVs for likely noise without filtering by low VAF."
    )
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {VERSION}")
    parser.add_argument("input", help="Merged TSV from vcfcalls2tsv.py")
    parser.add_argument("-o", "--output", required=True, help="Output TSV with tier/flag columns")
    parser.add_argument("--hotspots", required=True,
                        help="hg38 hotspot list (tab-delimited); matched additively by "
                             "gene + protein change AND by minimal genomic coordinates")
    parser.add_argument("--high-confidence-output", help="Optional TSV containing only HIGH_CONFIDENCE rows")
    parser.add_argument("--review-output", help="Optional TSV containing only REVIEW rows")
    parser.add_argument("--single-caller-output",
                        help="Optional TSV containing only SINGLE_CALLER rows (HC bar, one caller)")
    parser.add_argument("--noise-output", help="Optional TSV containing LIKELY_NOISE rows")
    parser.add_argument("--min-alt-reads-review", type=int, default=3,
                        help="Alt-read threshold used for low-evidence flags; does not use AF")
    parser.add_argument("--min-alt-reads-strong", type=int, default=5,
                        help="Alt-read threshold used as a rescue/evidence flag")
    parser.add_argument("--apply-clinic-floors", action="store_true",
                        help="Enable tissue-style clinic VAF/AD/DP floor flags and allow "
                             "them to demote single-caller rows; default is off (cfDNA)")
    parser.add_argument("--vaf-somatic-max", type=float, default=VAF_SOMATIC_MAX_DEFAULT,
                        help="Tumor VAF ceiling treated as somatic-leaning and HC-eligible "
                             f"(default: {VAF_SOMATIC_MAX_DEFAULT}). VAFs above this are labeled "
                             "het/atypical/hom-like and capped at REVIEW. Does not drop variants.")
    args = parser.parse_args()
    try:
        hotspots = load_hotspots(args.hotspots)
        counts = write_filtered(
            args.input,
            args.output,
            args.high_confidence_output,
            args.review_output,
            args.single_caller_output,
            args.noise_output,
            args.min_alt_reads_review,
            args.min_alt_reads_strong,
            hotspots,
            args.apply_clinic_floors,
            args.vaf_somatic_max,
        )
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)
    print(f"Tiered TSV written to {args.output}")
    for tier in ["HIGH_CONFIDENCE", "SINGLE_CALLER", "REVIEW", "LIKELY_NOISE"]:
        print(f"{tier}: {counts.get(tier, 0)}")

if __name__ == "__main__":
    main()
