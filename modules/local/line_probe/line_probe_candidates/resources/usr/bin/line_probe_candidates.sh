#!/usr/bin/env bash
# Cluster lp_hg38_discordant anchors (150 bp), annotate nearest L1 with bedtools (500 bp), score.
# Usage: line_probe_candidates.sh <bam> <line_l1.bed> <prefix> <bam_layer> <suffix>
set -euo pipefail

BAM="$1"
LINE_ANNOT="$2"
PREFIX="$3"
BAM_LAYER="$4"
SUFFIX="$5"

L1_MAX_DIST=500
CLUSTER_DIST=150
OUT="${PREFIX}.line_probe.candidates.${SUFFIX}.tsv"

HEADER='sample	bam_layer	cluster_id	chr	pos	n_discordant_reads	n_lp_alignment	n_hg38_mate_on_lp	n_secondary	n_supplementary	n_sa_tag	mean_mapq	max_mapq	anchor_mean_depth	anchor_max_depth	nearest_l1_name	nearest_l1_dist_bp	nearest_l1_strand	candidate_score'

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

sort -k2,2V -k3,3n "${PREFIX}.discordant_reads.tsv" | awk -v OFS='\t' -v cluster_dist="${CLUSTER_DIST}" '
function emit() {
    mean_mapq = (n > 0) ? sum_mapq / n : 0
    print cid, c_chr, c_pos, n, n_lp, n_hg38, n_sec, n_sup, n_sa, mean_mapq, max_mapq
}
BEGIN { cid = 0 }
{
    chr = $2; pos = $3; mapq = $4; lp = $5; hg38 = $6; sec = $7; sup = $8; sa = $9
    if (cid == 0 || chr != c_chr || pos - c_pos > cluster_dist) {
        if (cid > 0) emit()
        cid++
        c_chr = chr
        c_pos = pos
        n = 0; n_lp = 0; n_hg38 = 0; n_sec = 0; n_sup = 0; n_sa = 0
        sum_mapq = 0; max_mapq = 0
    }
    n++
    n_lp += lp
    n_hg38 += hg38
    n_sec += sec
    n_sup += sup
    n_sa += sa
    sum_mapq += mapq
    if (mapq > max_mapq) max_mapq = mapq
}
END { if (cid > 0) emit() }
' > "${PREFIX}.clusters.tsv"

awk -v OFS='\t' '{ print $2, $3 - 1, $3, "cluster_" $1 }' "${PREFIX}.clusters.tsv" > "${PREFIX}.clusters.bed"

LINE_ANNOT_SORTED="${PREFIX}.line_l1.sorted.bed"
bedtools sort -i "$LINE_ANNOT" > "$LINE_ANNOT_SORTED"
bedtools sort -i "${PREFIX}.clusters.bed" > "${PREFIX}.clusters.sorted.bed"

bedtools closest \
    -a "${PREFIX}.clusters.sorted.bed" \
    -b "$LINE_ANNOT_SORTED" \
    -d \
    -t first \
    -k 1 \
    > "${PREFIX}.l1_closest.tsv"

n_clusters=$(wc -l < "${PREFIX}.clusters.tsv" | tr -d ' ')
n_closest=$(wc -l < "${PREFIX}.l1_closest.tsv" | tr -d ' ')
if [[ "$n_clusters" -ne "$n_closest" ]]; then
    echo "ERROR: cluster count (${n_clusters}) != bedtools closest lines (${n_closest})" >&2
    exit 1
fi

printf '%s\n' "$HEADER" > "$OUT"

exec 3< "${PREFIX}.clusters.tsv"
exec 4< "${PREFIX}.l1_closest.tsv"
while IFS=$'\t' read -r cid chr pos n n_lp n_hg38 n_sec n_sup n_sa mean_mapq max_mapq <&3; do
    if ! IFS=$'\t' read -r _bchr _bstart _bend _cname l1_chr l1_start l1_end l1_name _l1_score l1_strand l1_dist <&4; then
        echo "ERROR: missing bedtools closest line for cluster ${cid}" >&2
        exit 1
    fi

    if [[ -z "${l1_dist:-}" || "$l1_dist" == "." || "$l1_dist" == "-1" || "${l1_dist%%.*}" -gt "${L1_MAX_DIST}" ]]; then
        l1_name="."
        l1_dist="."
        l1_strand_out="."
    else
        l1_strand_out="$l1_strand"
    fi

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

    score=$(awk -v dist="${l1_dist:-999999}" -v n="$n" -v n_sa="$n_sa" -v n_sup="$n_sup" \
        -v mean_mapq="$mean_mapq" -v max_depth="$anchor_max_depth" \
        -v l1_max="${L1_MAX_DIST}" 'BEGIN {
        if (dist == "." || dist == "") {
            l1_pts = 0
        } else if (dist + 0 <= l1_max) {
            l1_pts = int(40 * (1 - (dist + 0) / l1_max))
        } else {
            l1_pts = 0
        }
        read_pts = (n >= 10) ? 30 : int(n * 3)
        split_pts = (n_sa > 0 || n_sup > 0) ? 15 : 0
        mapq_pts = int(mean_mapq * 15 / 60)
        if (mapq_pts > 15) mapq_pts = 15
        depth_pts = int(max_depth)
        if (depth_pts > 15) depth_pts = 15
        printf "%d", l1_pts + read_pts + split_pts + mapq_pts + depth_pts
    }')

    printf '%s\n' \
        "${PREFIX}" "${BAM_LAYER}" "cluster_${cid}" "${chr}" "${pos}" "${n}" "${n_lp}" "${n_hg38}" \
        "${n_sec}" "${n_sup}" "${n_sa}" "${mean_mapq}" "${max_mapq}" \
        "${anchor_mean_depth}" "${anchor_max_depth}" "${l1_name}" "${l1_dist}" "${l1_strand_out}" "${score}" \
        >> "$OUT"
done
exec 3<&-
exec 4<&-

rm -f "${PREFIX}.discordant_reads.tsv" "${PREFIX}.clusters.tsv" "${PREFIX}.clusters.bed" \
    "${PREFIX}.clusters.sorted.bed" "${PREFIX}.line_l1.sorted.bed" "${PREFIX}.l1_closest.tsv"
