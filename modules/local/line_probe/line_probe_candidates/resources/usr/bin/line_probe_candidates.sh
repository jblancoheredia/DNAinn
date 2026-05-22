#!/usr/bin/env bash
# Per-read discordant anchors: label by BED overlap (not_in_bed / ambiguous / L1 name).
# Usage: line_probe_candidates.sh <bam> <line_l1.bed> <prefix> <bam_layer> <suffix>
set -euo pipefail

BAM="$1"
LINE_ANNOT="$2"
PREFIX="$3"
BAM_LAYER="$4"
SUFFIX="$5"

OUT="${PREFIX}.line_probe.candidates.${SUFFIX}.tsv"

HEADER='sample	bam_layer	read_id	chr	pos	n_discordant_reads	n_lp_alignment	n_hg38_mate_on_lp	n_secondary	n_supplementary	n_sa_tag	mean_mapq	max_mapq	anchor_mean_depth	anchor_max_depth	nearest_l1_name	n_l1_overlaps	nearest_l1_strand	candidate_score'

samtools view -F 4 "$BAM" | awk -v OFS='\t' '
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
    printf '%s\n' "$HEADER" > "$OUT"
    exit 0
fi

awk -v OFS='\t' '{ print $2, $3 - 1, $3, $1 }' "${PREFIX}.discordant_reads.tsv" > "${PREFIX}.anchors.bed"

LINE_ANNOT_SORTED="${PREFIX}.line_l1.sorted.bed"
bedtools sort -i "$LINE_ANNOT" > "$LINE_ANNOT_SORTED"
bedtools sort -i "${PREFIX}.anchors.bed" > "${PREFIX}.anchors.sorted.bed"

bedtools intersect \
    -a "${PREFIX}.anchors.sorted.bed" \
    -b "$LINE_ANNOT_SORTED" \
    -wa -wb \
    > "${PREFIX}.l1_overlaps.tsv" || true

# read_id -> label (not_in_bed | ambiguous | single L1 name); n_overlaps; strand
awk -F'\t' -v OFS='\t' '
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

printf '%s\n' "$HEADER" > "$OUT"

declare -A READ_L1_NAME READ_N_OVERLAPS READ_L1_STRAND
if [[ -s "${PREFIX}.read_labels.tsv" ]]; then
    while IFS=$'\t' read -r rid lab nov st; do
        READ_L1_NAME["$rid"]="$lab"
        READ_N_OVERLAPS["$rid"]="$nov"
        READ_L1_STRAND["$rid"]="$st"
    done < "${PREFIX}.read_labels.tsv"
fi

while IFS=$'\t' read -r read_id chr pos mapq lp hg38 sec sup sa; do
    l1_name="${READ_L1_NAME[$read_id]:-not_in_bed}"
    n_overlaps="${READ_N_OVERLAPS[$read_id]:-0}"
    l1_strand="${READ_L1_STRAND[$read_id]:-.}"

    region_start=$((pos - 100))
    [[ $region_start -lt 0 ]] && region_start=0
    region_end=$((pos + 100))

    depth_stats=$(samtools depth -aa -r "${chr}:${region_start}-${region_end}" "$BAM" | awk '
        { s += $3; if ($3 > m) m = $3; c++ }
        END {
            if (c == 0) print "0\t0"
            else printf "%.2f\t%d", s / c, m
        }
    ')
    anchor_mean_depth=$(echo "$depth_stats" | cut -f1)
    anchor_max_depth=$(echo "$depth_stats" | cut -f2)

    score=$(awk -v label="$l1_name" -v n_sa="$sa" -v n_sup="$sup" \
        -v mapq="$mapq" -v max_depth="$anchor_max_depth" 'BEGIN {
        if (label != "not_in_bed" && label != "ambiguous") l1_pts = 40
        else l1_pts = 0
        read_pts = 3
        split_pts = (n_sa > 0 || n_sup > 0) ? 15 : 0
        mapq_pts = int(mapq * 15 / 60)
        if (mapq_pts > 15) mapq_pts = 15
        depth_pts = int(max_depth)
        if (depth_pts > 15) depth_pts = 15
        printf "%d", l1_pts + read_pts + split_pts + mapq_pts + depth_pts
    }')

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${PREFIX}" "${BAM_LAYER}" "${read_id}" "${chr}" "${pos}" "1" "${lp}" "${hg38}" \
        "${sec}" "${sup}" "${sa}" "${mapq}" "${mapq}" \
        "${anchor_mean_depth}" "${anchor_max_depth}" "${l1_name}" "${n_overlaps}" "${l1_strand}" "${score}" \
        >> "$OUT"
done < "${PREFIX}.discordant_reads.tsv"

rm -f "${PREFIX}.discordant_reads.tsv" "${PREFIX}.anchors.bed" "${PREFIX}.anchors.sorted.bed" \
    "${PREFIX}.line_l1.sorted.bed" "${PREFIX}.l1_overlaps.tsv" "${PREFIX}.read_labels.tsv"
