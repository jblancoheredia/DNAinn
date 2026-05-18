/*
 * LP coverage depth, alignment-class counts, and hg38 anchor depth on realigned BAMs.
 */

process LINE_PROBE_QC_RAW {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::samtools=1.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/umi_aligner:1.0.0' :
        'blancojmskcc/umi_aligner:1.0.0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.line_probe.depth.raw.tsv")            , emit: depth
    tuple val(meta), path("*.line_probe.counts.raw.tsv")           , emit: counts
    tuple val(meta), path("*.line_probe.hg38_anchor.depth.raw.tsv"), emit: hg38_anchor_depth
    path "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    samtools depth -aa -r LP:1-120 ${bam} > ${prefix}.line_probe.depth.raw.tsv

    LP_CONCORDANT_COUNT=\$(samtools view -F 2304 ${bam} \\
        | awk '\$3=="LP" && (\$7=="LP" || \$7=="=")' | wc -l)

    LP_DISCORDANT_COUNT=\$(samtools view -F 4 ${bam} \\
        | awk '(\$3!="LP" && \$7=="LP") || (\$3=="LP" && \$7!="LP" && \$7!="=" && \$7!="*")' | wc -l)

    printf "sample\tbam_layer\talignment_class\tread_count\n" > ${prefix}.line_probe.counts.raw.tsv
    printf "${prefix}\tRaw\tlp_concordant_pair\t\${LP_CONCORDANT_COUNT}\n" >> ${prefix}.line_probe.counts.raw.tsv
    printf "${prefix}\tRaw\tlp_hg38_discordant\t\${LP_DISCORDANT_COUNT}\n" >> ${prefix}.line_probe.counts.raw.tsv

    samtools view -h -F 4 ${bam} \\
        | awk '\$1~/^@/ || (\$3!="LP" && \$7=="LP") || (\$3=="LP" && \$7!="LP" && \$7!="=" && \$7!="*")' \\
        | samtools view -bS - -o ${prefix}.lp_discordant_reads.bam

    samtools index ${prefix}.lp_discordant_reads.bam

    samtools view -F 4 ${prefix}.lp_discordant_reads.bam \\
        | awk 'OFS="\t" {
            if (\$3 == "LP") {
                chr = \$7
                pos = \$8
            } else {
                chr = \$3
                pos = \$4
            }
            start = pos - 1 - 100
            if (start < 0) start = 0
            end = pos - 1 + 100
            print chr, start, end
        }' \\
        | sort -k1,1 -k2,2n \\
        | uniq > ${prefix}.hg38_anchor_sites.bed

    printf "hg38_region\tchr\tpos\tdepth\n" > ${prefix}.line_probe.hg38_anchor.depth.raw.tsv

    while IFS=\$'\t' read -r CHR START END; do
        SITE="\${CHR}:\$((START+1))-\${END}"
        REGION="\${CHR}:\$((START+1))-\${END}"
        samtools depth -aa -r "\$REGION" ${bam} \\
            | awk -v site="\$SITE" 'OFS="\t" {print site, \$1, \$2, \$3}'
    done < ${prefix}.hg38_anchor_sites.bed >> ${prefix}.line_probe.hg38_anchor.depth.raw.tsv

    rm -f ${prefix}.lp_discordant_reads.bam ${prefix}.lp_discordant_reads.bam.bai ${prefix}.hg38_anchor_sites.bed

    cat <<END_VERSIONS > versions.yml
"${task.process}":
    samtools: \$(samtools --version | head -1 | sed 's/samtools //')
END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.line_probe.depth.raw.tsv
    touch ${prefix}.line_probe.counts.raw.tsv
    touch ${prefix}.line_probe.hg38_anchor.depth.raw.tsv

    cat <<END_VERSIONS > versions.yml
"${task.process}":
    samtools: \$(samtools --version | head -1 | sed 's/samtools //')
END_VERSIONS
    """
}

process LINE_PROBE_QC_CON {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::samtools=1.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/umi_aligner:1.0.0' :
        'blancojmskcc/umi_aligner:1.0.0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.line_probe.depth.con.tsv")            , emit: depth
    tuple val(meta), path("*.line_probe.counts.con.tsv")           , emit: counts
    tuple val(meta), path("*.line_probe.hg38_anchor.depth.con.tsv"), emit: hg38_anchor_depth
    path "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    samtools depth -aa -r LP:1-120 ${bam} > ${prefix}.line_probe.depth.con.tsv

    LP_CONCORDANT_COUNT=\$(samtools view -F 2304 ${bam} \\
        | awk '\$3=="LP" && (\$7=="LP" || \$7=="=")' | wc -l)

    LP_DISCORDANT_COUNT=\$(samtools view -F 4 ${bam} \\
        | awk '(\$3!="LP" && \$7=="LP") || (\$3=="LP" && \$7!="LP" && \$7!="=" && \$7!="*")' | wc -l)

    printf "sample\tbam_layer\talignment_class\tread_count\n" > ${prefix}.line_probe.counts.con.tsv
    printf "${prefix}\tConsensus\tlp_concordant_pair\t\${LP_CONCORDANT_COUNT}\n" >> ${prefix}.line_probe.counts.con.tsv
    printf "${prefix}\tConsensus\tlp_hg38_discordant\t\${LP_DISCORDANT_COUNT}\n" >> ${prefix}.line_probe.counts.con.tsv

    samtools view -h -F 4 ${bam} \\
        | awk '\$1~/^@/ || (\$3!="LP" && \$7=="LP") || (\$3=="LP" && \$7!="LP" && \$7!="=" && \$7!="*")' \\
        | samtools view -bS - -o ${prefix}.lp_discordant_reads.bam

    samtools index ${prefix}.lp_discordant_reads.bam

    samtools view -F 4 ${prefix}.lp_discordant_reads.bam \\
        | awk 'OFS="\t" {
            if (\$3 == "LP") {
                chr = \$7
                pos = \$8
            } else {
                chr = \$3
                pos = \$4
            }
            start = pos - 1 - 100
            if (start < 0) start = 0
            end = pos - 1 + 100
            print chr, start, end
        }' \\
        | sort -k1,1 -k2,2n \\
        | uniq > ${prefix}.hg38_anchor_sites.bed

    printf "hg38_region\tchr\tpos\tdepth\n" > ${prefix}.line_probe.hg38_anchor.depth.con.tsv

    while IFS=\$'\t' read -r CHR START END; do
        SITE="\${CHR}:\$((START+1))-\${END}"
        REGION="\${CHR}:\$((START+1))-\${END}"
        samtools depth -aa -r "\$REGION" ${bam} \\
            | awk -v site="\$SITE" 'OFS="\t" {print site, \$1, \$2, \$3}'
    done < ${prefix}.hg38_anchor_sites.bed >> ${prefix}.line_probe.hg38_anchor.depth.con.tsv

    rm -f ${prefix}.lp_discordant_reads.bam ${prefix}.lp_discordant_reads.bam.bai ${prefix}.hg38_anchor_sites.bed

    cat <<END_VERSIONS > versions.yml
"${task.process}":
    samtools: \$(samtools --version | head -1 | sed 's/samtools //')
END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.line_probe.depth.con.tsv
    touch ${prefix}.line_probe.counts.con.tsv
    touch ${prefix}.line_probe.hg38_anchor.depth.con.tsv

    cat <<END_VERSIONS > versions.yml
"${task.process}":
    samtools: \$(samtools --version | head -1 | sed 's/samtools //')
END_VERSIONS
    """
}

process LINE_PROBE_QC_DUP {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::samtools=1.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/umi_aligner:1.0.0' :
        'blancojmskcc/umi_aligner:1.0.0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.line_probe.depth.dup.tsv")            , emit: depth
    tuple val(meta), path("*.line_probe.counts.dup.tsv")           , emit: counts
    tuple val(meta), path("*.line_probe.hg38_anchor.depth.dup.tsv"), emit: hg38_anchor_depth
    path "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    samtools depth -aa -r LP:1-120 ${bam} > ${prefix}.line_probe.depth.dup.tsv

    LP_CONCORDANT_COUNT=\$(samtools view -F 2304 ${bam} \\
        | awk '\$3=="LP" && (\$7=="LP" || \$7=="=")' | wc -l)

    LP_DISCORDANT_COUNT=\$(samtools view -F 4 ${bam} \\
        | awk '(\$3!="LP" && \$7=="LP") || (\$3=="LP" && \$7!="LP" && \$7!="=" && \$7!="*")' | wc -l)

    printf "sample\tbam_layer\talignment_class\tread_count\n" > ${prefix}.line_probe.counts.dup.tsv
    printf "${prefix}\tDuplex\tlp_concordant_pair\t\${LP_CONCORDANT_COUNT}\n" >> ${prefix}.line_probe.counts.dup.tsv
    printf "${prefix}\tDuplex\tlp_hg38_discordant\t\${LP_DISCORDANT_COUNT}\n" >> ${prefix}.line_probe.counts.dup.tsv

    samtools view -h -F 4 ${bam} \\
        | awk '\$1~/^@/ || (\$3!="LP" && \$7=="LP") || (\$3=="LP" && \$7!="LP" && \$7!="=" && \$7!="*")' \\
        | samtools view -bS - -o ${prefix}.lp_discordant_reads.bam

    samtools index ${prefix}.lp_discordant_reads.bam

    samtools view -F 4 ${prefix}.lp_discordant_reads.bam \\
        | awk 'OFS="\t" {
            if (\$3 == "LP") {
                chr = \$7
                pos = \$8
            } else {
                chr = \$3
                pos = \$4
            }
            start = pos - 1 - 100
            if (start < 0) start = 0
            end = pos - 1 + 100
            print chr, start, end
        }' \\
        | sort -k1,1 -k2,2n \\
        | uniq > ${prefix}.hg38_anchor_sites.bed

    printf "hg38_region\tchr\tpos\tdepth\n" > ${prefix}.line_probe.hg38_anchor.depth.dup.tsv

    while IFS=\$'\t' read -r CHR START END; do
        SITE="\${CHR}:\$((START+1))-\${END}"
        REGION="\${CHR}:\$((START+1))-\${END}"
        samtools depth -aa -r "\$REGION" ${bam} \\
            | awk -v site="\$SITE" 'OFS="\t" {print site, \$1, \$2, \$3}'
    done < ${prefix}.hg38_anchor_sites.bed >> ${prefix}.line_probe.hg38_anchor.depth.dup.tsv

    rm -f ${prefix}.lp_discordant_reads.bam ${prefix}.lp_discordant_reads.bam.bai ${prefix}.hg38_anchor_sites.bed

    cat <<END_VERSIONS > versions.yml
"${task.process}":
    samtools: \$(samtools --version | head -1 | sed 's/samtools //')
END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.line_probe.depth.dup.tsv
    touch ${prefix}.line_probe.counts.dup.tsv
    touch ${prefix}.line_probe.hg38_anchor.depth.dup.tsv

    cat <<END_VERSIONS > versions.yml
"${task.process}":
    samtools: \$(samtools --version | head -1 | sed 's/samtools //')
END_VERSIONS
    """
}

process LINE_PROBE_QC_SIM {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::samtools=1.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/umi_aligner:1.0.0' :
        'blancojmskcc/umi_aligner:1.0.0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.line_probe.depth.sim.tsv")            , emit: depth
    tuple val(meta), path("*.line_probe.counts.sim.tsv")           , emit: counts
    tuple val(meta), path("*.line_probe.hg38_anchor.depth.sim.tsv"), emit: hg38_anchor_depth
    path "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    samtools depth -aa -r LP:1-120 ${bam} > ${prefix}.line_probe.depth.sim.tsv

    LP_CONCORDANT_COUNT=\$(samtools view -F 2304 ${bam} \\
        | awk '\$3=="LP" && (\$7=="LP" || \$7=="=")' | wc -l)

    LP_DISCORDANT_COUNT=\$(samtools view -F 4 ${bam} \\
        | awk '(\$3!="LP" && \$7=="LP") || (\$3=="LP" && \$7!="LP" && \$7!="=" && \$7!="*")' | wc -l)

    printf "sample\tbam_layer\talignment_class\tread_count\n" > ${prefix}.line_probe.counts.sim.tsv
    printf "${prefix}\tSimplex\tlp_concordant_pair\t\${LP_CONCORDANT_COUNT}\n" >> ${prefix}.line_probe.counts.sim.tsv
    printf "${prefix}\tSimplex\tlp_hg38_discordant\t\${LP_DISCORDANT_COUNT}\n" >> ${prefix}.line_probe.counts.sim.tsv

    samtools view -h -F 4 ${bam} \\
        | awk '\$1~/^@/ || (\$3!="LP" && \$7=="LP") || (\$3=="LP" && \$7!="LP" && \$7!="=" && \$7!="*")' \\
        | samtools view -bS - -o ${prefix}.lp_discordant_reads.bam

    samtools index ${prefix}.lp_discordant_reads.bam

    samtools view -F 4 ${prefix}.lp_discordant_reads.bam \\
        | awk 'OFS="\t" {
            if (\$3 == "LP") {
                chr = \$7
                pos = \$8
            } else {
                chr = \$3
                pos = \$4
            }
            start = pos - 1 - 100
            if (start < 0) start = 0
            end = pos - 1 + 100
            print chr, start, end
        }' \\
        | sort -k1,1 -k2,2n \\
        | uniq > ${prefix}.hg38_anchor_sites.bed

    printf "hg38_region\tchr\tpos\tdepth\n" > ${prefix}.line_probe.hg38_anchor.depth.sim.tsv

    while IFS=\$'\t' read -r CHR START END; do
        SITE="\${CHR}:\$((START+1))-\${END}"
        REGION="\${CHR}:\$((START+1))-\${END}"
        samtools depth -aa -r "\$REGION" ${bam} \\
            | awk -v site="\$SITE" 'OFS="\t" {print site, \$1, \$2, \$3}'
    done < ${prefix}.hg38_anchor_sites.bed >> ${prefix}.line_probe.hg38_anchor.depth.sim.tsv

    rm -f ${prefix}.lp_discordant_reads.bam ${prefix}.lp_discordant_reads.bam.bai ${prefix}.hg38_anchor_sites.bed

    cat <<END_VERSIONS > versions.yml
"${task.process}":
    samtools: \$(samtools --version | head -1 | sed 's/samtools //')
END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.line_probe.depth.sim.tsv
    touch ${prefix}.line_probe.counts.sim.tsv
    touch ${prefix}.line_probe.hg38_anchor.depth.sim.tsv

    cat <<END_VERSIONS > versions.yml
"${task.process}":
    samtools: \$(samtools --version | head -1 | sed 's/samtools //')
END_VERSIONS
    """
}
