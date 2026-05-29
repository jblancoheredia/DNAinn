#!/usr/bin/env bash
# Discordant anchors: overlap-label reads, then group by exact anchor+label (chr,pos,label).
# Usage: line_probe_candidates.sh <bam> <line_l1.bed> <prefix> <bam_layer> <suffix>
set -euo pipefail

BAM="$1"
LINE_ANNOT="$2"
PREFIX="$3"
BAM_LAYER="$4"
SUFFIX="$5"

OUT="${PREFIX}.line_probe.candidates.${SUFFIX}.tsv"

HEADER='sample	bam_layer	read_id	chr	pos	n_discordant_reads	n_lp_alignment	n_hg38_mate_on_lp	n_secondary	n_supplementary	n_sa_tag	mean_mapq	max_mapq	anchor_mean_depth	anchor_max_depth	nearest_l1_name	n_l1_overlaps	nearest_l1_strand	candidate_score'

samtools view -F 4 "$BAM" | awk -v OFS='	' '
{
    disc = ($3 != "LP" && $7 == "LP") || ($3 == "LP" && $7 != "LP" && $7 != "=" && $7 != "*")
    if (!disc) next
    if ($3 == "LP") {
        chr = $7
        pos = $8
        lp_align = 1
        hg38_mate = 0
    } else {
        chr = $3
        pos = $4
        lp_align = 0
        hg38_mate = 1
    }
    if (chr == "*" || chr == "LP" || pos < 1) next
    sec = (int($2 / 256) % 2 == 1) ? 1 : 0
    sup = (int($2 / 2048) % 2 == 1) ? 1 : 0
    sa = ($0 ~ /SA:Z:/) ? 1 : 0
    print $1, chr, pos, $5, lp_align, hg38_mate, sec, sup, sa
}' > "${PREFIX}.discordant_reads.tsv"

if [[ ! -s "${PREFIX}.discordant_reads.tsv" ]]; then
    printf '%s
' "$HEADER" > "$OUT"
    exit 0
fi

awk -v OFS='	' '{ print $2, $3 - 1, $3, $1 }' "${PREFIX}.discordant_reads.tsv" > "${PREFIX}.anchors.bed"

LINE_ANNOT_SORTED="${PREFIX}.line_l1.sorted.bed"
bedtools sort -i "$LINE_ANNOT" > "$LINE_ANNOT_SORTED"
bedtools sort -i "${PREFIX}.anchors.bed" > "${PREFIX}.anchors.sorted.bed"

bedtools intersect     -a "${PREFIX}.anchors.sorted.bed"     -b "$LINE_ANNOT_SORTED"     -wa -wb     > "${PREFIX}.l1_overlaps.tsv" || true

# read_id -> label (not_in_bed | ambiguous | single L1 name); n_overlaps; strand
awk -F'	' -v OFS='	' '
{
    rid = $4
    l1 = $8
    strand = $10
    key = rid SUBSEP l1
    if (key in seen) next
    seen[key] = 1
    n_ov[rid]++
    if (n_names[rid] == 0) {
        only[rid] = l1
        only_strand[rid] = strand
    } else if (only[rid] != l1) {
        only[rid] = ""
        only_strand[rid] = "."
    }
    n_names[rid]++
}
END {
    for (rid in n_ov) {
        if (only[rid] == "") label = "ambiguous"
        else label = only[rid]
        print rid, label, n_ov[rid], only_strand[rid]
    }
}
' "${PREFIX}.l1_overlaps.tsv" > "${PREFIX}.read_labels.tsv"

# Group by exact anchor+label: chr, pos, nearest_l1_name
awk -F'	' -v OFS='	' -v labels_file="${PREFIX}.read_labels.tsv" '
BEGIN {
    while ((getline < labels_file) > 0) {
        rid = $1
        label[rid] = $2
        nov[rid] = $3
        strand[rid] = $4
    }
    close(labels_file)
}
{
    rid = $1; chr = $2; pos = $3; mapq = $4; lp = $5; hg38 = $6; sec = $7; sup = $8; sa = $9
    lab = (rid in label) ? label[rid] : "not_in_bed"
    n_ov = (rid in nov) ? nov[rid] : 0
    st = (rid in strand) ? strand[rid] : "."

    key = chr SUBSEP pos SUBSEP lab
    n[key]++
    n_lp[key] += lp
    n_hg38[key] += hg38
    n_sec[key] += sec
    n_sup[key] += sup
    n_sa[key] += sa
    sum_mapq[key] += mapq
    if (!(key in max_mapq) || mapq > max_mapq[key]) max_mapq[key] = mapq
    if (!(key in max_ov) || n_ov > max_ov[key]) max_ov[key] = n_ov
    if (!(key in strand_out) || strand_out[key] == ".") strand_out[key] = st
    if (lab == "ambiguous" || lab == "not_in_bed") strand_out[key] = "."
}
END {
    for (k in n) {
        split(k, a, SUBSEP)
        chr = a[1]; pos = a[2]; lab = a[3]
        mean_mapq = (n[k] > 0) ? (sum_mapq[k] / n[k]) : 0
        event_id = chr ":" pos ":" lab
        print event_id, chr, pos, n[k], n_lp[k], n_hg38[k], n_sec[k], n_sup[k], n_sa[k], mean_mapq, max_mapq[k], max_ov[k], strand_out[k], lab
    }
}
' "${PREFIX}.discordant_reads.tsv" > "${PREFIX}.grouped.tsv"

printf '%s
' "$HEADER" > "$OUT"

while IFS=$'	' read -r event_id chr pos n n_lp n_hg38 n_sec n_sup n_sa mean_mapq max_mapq n_overlaps l1_strand l1_name; do
    region_start=$((pos - 100))
    [[ $region_start -lt 0 ]] && region_start=0
    region_end=$((pos + 100))

    depth_stats=$(samtools depth -aa -r "${chr}:${region_start}-${region_end}" "$BAM" | awk '
        { s += $3; if ($3 > m) m = $3; c++ }
        END {
            if (c == 0) print "0	0"
            else printf "%.2f	%d", s / c, m
        }
    ')
    anchor_mean_depth=$(echo "$depth_stats" | cut -f1)
    anchor_max_depth=$(echo "$depth_stats" | cut -f2)

    score=$(awk -v label="$l1_name" -v n="$n" -v n_sa="$n_sa" -v n_sup="$n_sup"         -v mean_mapq="$mean_mapq" -v max_depth="$anchor_max_depth" 'BEGIN {
        if (label != "not_in_bed" && label != "ambiguous") l1_pts = 40
        else l1_pts = 0
        read_pts = (n >= 10) ? 30 : int(n * 3)
        split_pts = (n_sa > 0 || n_sup > 0) ? 15 : 0
        mapq_pts = int(mean_mapq * 15 / 60)
        if (mapq_pts > 15) mapq_pts = 15
        depth_pts = int(max_depth)
        if (depth_pts > 15) depth_pts = 15
        printf "%d", l1_pts + read_pts + split_pts + mapq_pts + depth_pts
    }')

    printf '%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s
'         "${PREFIX}" "${BAM_LAYER}" "${event_id}" "${chr}" "${pos}" "${n}" "${n_lp}" "${n_hg38}"         "${n_sec}" "${n_sup}" "${n_sa}" "${mean_mapq}" "${max_mapq}"         "${anchor_mean_depth}" "${anchor_max_depth}" "${l1_name}" "${n_overlaps}" "${l1_strand}" "${score}"         >> "$OUT"
done < "${PREFIX}.grouped.tsv"

rm -f "${PREFIX}.discordant_reads.tsv" "${PREFIX}.anchors.bed" "${PREFIX}.anchors.sorted.bed"     "${PREFIX}.line_l1.sorted.bed" "${PREFIX}.l1_overlaps.tsv" "${PREFIX}.read_labels.tsv" "${PREFIX}.grouped.tsv"
