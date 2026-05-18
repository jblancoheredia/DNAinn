#!/usr/bin/env bash
# Cluster lp_hg38_discordant anchors (150 bp), annotate nearest L1 (500 bp), score candidates.
# Usage: line_probe_candidates.sh <bam> <line_l1.bed> <prefix> <bam_layer> <suffix>
#
# candidate_score (0-100, higher = stronger candidate):
#   L1 proximity: up to 40 (linear decay from anchor to 500 bp)
#   discordant read support: up to 30 (3 points/read, cap at 10 reads)
#   split/supplementary/SA: 15 if any
#   mean MAPQ: up to 15
#   anchor max depth: up to 15
set -euo pipefail

BAM="$1"
LINE_ANNOT="$2"
PREFIX="$3"
BAM_LAYER="$4"
SUFFIX="$5"

L1_MAX_DIST=500
CLUSTER_DIST=150
OUT="${PREFIX}.line_probe.candidates.${SUFFIX}.tsv"

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
    if (chr == "*" || pos < 1) next
    sec = and($2, 256) ? 1 : 0
    sup = and($2, 2048) ? 1 : 0
    sa = ($0 ~ /SA:Z:/) ? 1 : 0
    print $1, chr, pos, $5, lp_align, hg38_mate, sec, sup, sa
}' > "${PREFIX}.discordant_reads.tsv"

if [[ ! -s "${PREFIX}.discordant_reads.tsv" ]]; then
    printf '%s\n' \
        'sample	bam_layer	cluster_id	chr	pos	n_discordant_reads	n_lp_alignment	n_hg38_mate_on_lp	n_secondary	n_supplementary	n_sa_tag	mean_mapq	max_mapq	anchor_mean_depth	anchor_max_depth	nearest_l1_name	nearest_l1_dist_bp	nearest_l1_strand	candidate_score' \
        > "$OUT"
    exit 0
fi

sort -k2,2 -k3,3n "${PREFIX}.discordant_reads.tsv" | awk -v OFS='\t' -v cluster_dist="${CLUSTER_DIST}" '
function emit(    i) {
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

awk -v OFS='\t' '{
    print "cluster_" $1, $2, $3 - 1, $3
}' "${PREFIX}.clusters.tsv" > "${PREFIX}.clusters.bed"

bedtools closest -a "${PREFIX}.clusters.bed" -b "$LINE_ANNOT" -d -t first -k 1 > "${PREFIX}.l1_closest.tsv"

printf '%s\n' \
    'sample	bam_layer	cluster_id	chr	pos	n_discordant_reads	n_lp_alignment	n_hg38_mate_on_lp	n_secondary	n_supplementary	n_sa_tag	mean_mapq	max_mapq	anchor_mean_depth	anchor_max_depth	nearest_l1_name	nearest_l1_dist_bp	nearest_l1_strand	candidate_score' \
    > "$OUT"

exec 3< "${PREFIX}.clusters.tsv"
exec 4< "${PREFIX}.l1_closest.tsv"
while IFS=$'\t' read -r cid chr pos n n_lp n_hg38 n_sec n_sup n_sa mean_mapq max_mapq <&3; do
    IFS=$'\t' read -r _cname _bchr _bstart _bend l1_chr l1_start l1_end l1_name _l1_score l1_strand l1_dist <&4 || true

    if [[ -z "${l1_dist:-}" || "$l1_dist" == "." || "$l1_dist" == "-1" ]]; then
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
        if (dist == "." || dist == "") l1_pts = 0
        else if (dist + 0 <= l1_max) l1_pts = int(40 * (1 - (dist + 0) / l1_max))
        else l1_pts = 0
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

rm -f "${PREFIX}.discordant_reads.tsv" "${PREFIX}.clusters.tsv" "${PREFIX}.clusters.bed" "${PREFIX}.l1_closest.tsv"
