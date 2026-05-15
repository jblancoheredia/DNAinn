process LINE_PROBE_RAW {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:03dc1d2818d9de56938078b8b78b82d967c1f820-0' :
        'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:03dc1d2818d9de56938078b8b78b82d967c1f820-0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path probe_fasta
    path bwa_dir

    output:
    tuple val(meta), path("*.line_probe.raw.bam")                , emit: bam
    tuple val(meta), path("*.line_probe.raw.bam.bai")            , emit: bai
    tuple val(meta), path("*.line_probe.depth.raw.tsv")          , emit: depth
    tuple val(meta), path("*.line_probe.counts.raw.tsv")         , emit: counts
    tuple val(meta), path("*.line_probe.insertion.depth.raw.tsv"), emit: summary
    path "versions.yml"                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-Y -M'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fasta_name = probe_fasta.getName()

    """
    BWA_INDEX_PREFIX=\$(find -L ./ -name "*.amb" | sed 's/.amb//')

    samtools collate -u -O ${bam} ${prefix}.collate \\
        | samtools fastq \\
            -n \\
            -F 0x900 \\
            -1 ${prefix}.R1.fastq \\
            -2 ${prefix}.R2.fastq \\
            -0 /dev/null \\
            -s /dev/null \\
            -

    bwa mem \\
        -t ${task.cpus} \\
        ${args} \\
        -R ${meta.read_group} \\
        \$BWA_INDEX_PREFIX \\
        ${prefix}.R1.fastq \\
        ${prefix}.R2.fastq \\
        | samtools sort -@ ${task.cpus} -o ${prefix}.line_probe.raw.bam

    samtools index ${prefix}.line_probe.raw.bam

    samtools depth -aa -H -r LP:1-120 ${prefix}.line_probe.raw.bam > ${prefix}.line_probe.depth.raw.tsv

    PROBE_COUNT=\$(samtools view -c -F 2304 \\
        -e 'rname=="LP" && rnext=="LP"' \\
        ${prefix}.line_probe.raw.bam)

    INSERT_COUNT=\$(samtools view -c -F 2304 \\
        -e 'rname!="LP" && rnext=="LP"' \\
        ${prefix}.line_probe.raw.bam)

    printf "sample\tbam_type\tcategory\tread_count\n" > ${prefix}.line_probe.counts.raw.tsv
    printf "${prefix}\traw\tprobe\t\${PROBE_COUNT}\n"    >> ${prefix}.line_probe.counts.raw.tsv
    printf "${prefix}\traw\tinsertion\t\${INSERT_COUNT}\n" >> ${prefix}.line_probe.counts.raw.tsv

    samtools view -F 2304 \\
        -e 'rname!="LP" && rnext=="LP"' \\
        -b ${prefix}.line_probe.raw.bam \\
        -o ${prefix}.insertion_reads.bam

    samtools index ${prefix}.insertion_reads.bam

    samtools view -F 2304 ${prefix}.insertion_reads.bam \\
        | awk 'OFS="\t" {
            start = (\$4 - 1) - 100
            if (start < 0) start = 0
            end = (\$4 - 1) + 100
            print \$3, start, end
        }' \\
        | sort -k1,1 -k2,2n \\
        | uniq > ${prefix}.insertion_sites.bed

    printf "site\tchr\tpos\tdepth\n" > ${prefix}.line_probe.insertion.depth.raw.tsv

    while IFS=\$'\t' read -r CHR START END; do
        SITE="\${CHR}:\$((START+1))-\${END}"
        REGION="\${CHR}:\$((START+1))-\${END}"
        samtools depth -aa -r "\$REGION" ${prefix}.insertion_reads.bam \\
            | awk -v site="\$SITE" 'OFS="\t" {print site, \$1, \$2, \$3}'
    done < ${prefix}.insertion_sites.bed >> ${prefix}.line_probe.insertion.depth.raw.tsv

    rm -f ${prefix}.insertion_reads.bam ${prefix}.insertion_reads.bam.bai ${prefix}.insertion_sites.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | tail -n 1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.line_probe.raw.bam
    touch ${prefix}.line_probe.raw.bam.bai
    touch ${prefix}.line_probe.depth.raw.tsv
    touch ${prefix}.line_probe.counts.raw.tsv
    touch ${prefix}.line_probe.insertion.depth.raw.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | tail -n 1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}

process LINE_PROBE_CON {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:03dc1d2818d9de56938078b8b78b82d967c1f820-0' :
        'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:03dc1d2818d9de56938078b8b78b82d967c1f820-0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path probe_fasta
    path bwa_dir

    output:
    tuple val(meta), path("*.line_probe.con.bam")                , emit: bam
    tuple val(meta), path("*.line_probe.con.bam.bai")            , emit: bai
    tuple val(meta), path("*.line_probe.depth.con.tsv")          , emit: depth
    tuple val(meta), path("*.line_probe.counts.con.tsv")         , emit: counts
    tuple val(meta), path("*.line_probe.insertion.depth.con.tsv"), emit: summary
    path "versions.yml"                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-Y -M'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fasta_name = probe_fasta.getName()

    """
    BWA_INDEX_PREFIX=\$(find -L ./ -name "*.amb" | sed 's/.amb//')

    samtools collate -u -O ${bam} ${prefix}.collate \\
        | samtools fastq \\
            -n \\
            -F 0x900 \\
            -1 ${prefix}.R1.fastq \\
            -2 ${prefix}.R2.fastq \\
            -0 /dev/null \\
            -s /dev/null \\
            -

    bwa mem \\
        -t ${task.cpus} \\
        ${args} \\
        -R ${meta.read_group} \\
        \$BWA_INDEX_PREFIX \\
        ${prefix}.R1.fastq \\
        ${prefix}.R2.fastq \\
        | samtools sort -@ ${task.cpus} -o ${prefix}.line_probe.con.bam

    samtools index ${prefix}.line_probe.con.bam

    samtools depth -aa -H -r LP:1-120 ${prefix}.line_probe.con.bam > ${prefix}.line_probe.depth.con.tsv

    PROBE_COUNT=\$(samtools view -c -F 2304 \\
        -e 'rname=="LP" && rnext=="LP"' \\
        ${prefix}.line_probe.con.bam)

    INSERT_COUNT=\$(samtools view -c -F 2304 \\
        -e 'rname!="LP" && rnext=="LP"' \\
        ${prefix}.line_probe.con.bam)

    printf "sample\tbam_type\tcategory\tread_count\n" > ${prefix}.line_probe.counts.con.tsv
    printf "${prefix}\tcon\tprobe\t\${PROBE_COUNT}\n"    >> ${prefix}.line_probe.counts.con.tsv
    printf "${prefix}\tcon\tinsertion\t\${INSERT_COUNT}\n" >> ${prefix}.line_probe.counts.con.tsv

    samtools view -F 2304 \\
        -e 'rname!="LP" && rnext=="LP"' \\
        -b ${prefix}.line_probe.con.bam \\
        -o ${prefix}.insertion_reads.bam

    samtools index ${prefix}.insertion_reads.bam

    samtools view -F 2304 ${prefix}.insertion_reads.bam \\
        | awk 'OFS="\t" {
            start = (\$4 - 1) - 100
            if (start < 0) start = 0
            end = (\$4 - 1) + 100
            print \$3, start, end
        }' \\
        | sort -k1,1 -k2,2n \\
        | uniq > ${prefix}.insertion_sites.bed

    printf "site\tchr\tpos\tdepth\n" > ${prefix}.line_probe.insertion.depth.con.tsv

    while IFS=\$'\t' read -r CHR START END; do
        SITE="\${CHR}:\$((START+1))-\${END}"
        REGION="\${CHR}:\$((START+1))-\${END}"
        samtools depth -aa -r "\$REGION" ${prefix}.insertion_reads.bam \\
            | awk -v site="\$SITE" 'OFS="\t" {print site, \$1, \$2, \$3}'
    done < ${prefix}.insertion_sites.bed >> ${prefix}.line_probe.insertion.depth.con.tsv

    rm -f ${prefix}.insertion_reads.bam ${prefix}.insertion_reads.bam.bai ${prefix}.insertion_sites.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | tail -n 1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.line_probe.con.bam
    touch ${prefix}.line_probe.con.bam.bai
    touch ${prefix}.line_probe.depth.con.tsv
    touch ${prefix}.line_probe.counts.con.tsv
    touch ${prefix}.line_probe.insertion.depth.con.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | tail -n 1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}

process LINE_PROBE_DUP {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:03dc1d2818d9de56938078b8b78b82d967c1f820-0' :
        'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:03dc1d2818d9de56938078b8b78b82d967c1f820-0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path probe_fasta
    path bwa_dir

    output:
    tuple val(meta), path("*.line_probe.dup.bam")                , emit: bam
    tuple val(meta), path("*.line_probe.dup.bam.bai")            , emit: bai
    tuple val(meta), path("*.line_probe.depth.dup.tsv")          , emit: depth
    tuple val(meta), path("*.line_probe.counts.dup.tsv")         , emit: counts
    tuple val(meta), path("*.line_probe.insertion.depth.dup.tsv"), emit: summary
    path "versions.yml"                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-Y -M'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fasta_name = probe_fasta.getName()

    """
    BWA_INDEX_PREFIX=\$(find -L ./ -name "*.amb" | sed 's/.amb//')

    samtools collate -u -O ${bam} ${prefix}.collate \\
        | samtools fastq \\
            -n \\
            -F 0x900 \\
            -1 ${prefix}.R1.fastq \\
            -2 ${prefix}.R2.fastq \\
            -0 /dev/null \\
            -s /dev/null \\
            -

    bwa mem \\
        -t ${task.cpus} \\
        ${args} \\
        -R ${meta.read_group} \\
        \$BWA_INDEX_PREFIX \\
        ${prefix}.R1.fastq \\
        ${prefix}.R2.fastq \\
        | samtools sort -@ ${task.cpus} -o ${prefix}.line_probe.dup.bam

    samtools index ${prefix}.line_probe.dup.bam

    samtools depth -aa -H -r LP:1-120 ${prefix}.line_probe.dup.bam > ${prefix}.line_probe.depth.dup.tsv

    PROBE_COUNT=\$(samtools view -c -F 2304 \\
        -e 'rname=="LP" && rnext=="LP"' \\
        ${prefix}.line_probe.dup.bam)

    INSERT_COUNT=\$(samtools view -c -F 2304 \\
        -e 'rname!="LP" && rnext=="LP"' \\
        ${prefix}.line_probe.dup.bam)

    printf "sample\tbam_type\tcategory\tread_count\n" > ${prefix}.line_probe.counts.dup.tsv
    printf "${prefix}\tdup\tprobe\t\${PROBE_COUNT}\n"    >> ${prefix}.line_probe.counts.dup.tsv
    printf "${prefix}\tdup\tinsertion\t\${INSERT_COUNT}\n" >> ${prefix}.line_probe.counts.dup.tsv

    samtools view -F 2304 \\
        -e 'rname!="LP" && rnext=="LP"' \\
        -b ${prefix}.line_probe.dup.bam \\
        -o ${prefix}.insertion_reads.bam

    samtools index ${prefix}.insertion_reads.bam

    samtools view -F 2304 ${prefix}.insertion_reads.bam \\
        | awk 'OFS="\t" {
            start = (\$4 - 1) - 100
            if (start < 0) start = 0
            end = (\$4 - 1) + 100
            print \$3, start, end
        }' \\
        | sort -k1,1 -k2,2n \\
        | uniq > ${prefix}.insertion_sites.bed

    printf "site\tchr\tpos\tdepth\n" > ${prefix}.line_probe.insertion.depth.dup.tsv

    while IFS=\$'\t' read -r CHR START END; do
        SITE="\${CHR}:\$((START+1))-\${END}"
        REGION="\${CHR}:\$((START+1))-\${END}"
        samtools depth -aa -r "\$REGION" ${prefix}.insertion_reads.bam \\
            | awk -v site="\$SITE" 'OFS="\t" {print site, \$1, \$2, \$3}'
    done < ${prefix}.insertion_sites.bed >> ${prefix}.line_probe.insertion.depth.dup.tsv

    rm -f ${prefix}.insertion_reads.bam ${prefix}.insertion_reads.bam.bai ${prefix}.insertion_sites.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | tail -n 1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.line_probe.dup.bam
    touch ${prefix}.line_probe.dup.bam.bai
    touch ${prefix}.line_probe.depth.dup.tsv
    touch ${prefix}.line_probe.counts.dup.tsv
    touch ${prefix}.line_probe.insertion.depth.dup.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | tail -n 1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}

process LINE_PROBE_SIM {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:03dc1d2818d9de56938078b8b78b82d967c1f820-0' :
        'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:03dc1d2818d9de56938078b8b78b82d967c1f820-0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path probe_fasta
    path bwa_dir

    output:
    tuple val(meta), path("*.line_probe.sim.bam")                , emit: bam
    tuple val(meta), path("*.line_probe.sim.bam.bai")            , emit: bai
    tuple val(meta), path("*.line_probe.depth.sim.tsv")          , emit: depth
    tuple val(meta), path("*.line_probe.counts.sim.tsv")         , emit: counts
    tuple val(meta), path("*.line_probe.insertion.depth.sim.tsv"), emit: summary
    path "versions.yml"                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-Y -M'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fasta_name = probe_fasta.getName()

    """
    BWA_INDEX_PREFIX=\$(find -L ./ -name "*.amb" | sed 's/.amb//')

    samtools collate -u -O ${bam} ${prefix}.collate \\
        | samtools fastq \\
            -n \\
            -F 0x900 \\
            -1 ${prefix}.R1.fastq \\
            -2 ${prefix}.R2.fastq \\
            -0 /dev/null \\
            -s /dev/null \\
            -

    bwa mem \\
        -t ${task.cpus} \\
        ${args} \\
        -R ${meta.read_group} \\
        \$BWA_INDEX_PREFIX \\
        ${prefix}.R1.fastq \\
        ${prefix}.R2.fastq \\
        | samtools sort -@ ${task.cpus} -o ${prefix}.line_probe.sim.bam

    samtools index ${prefix}.line_probe.sim.bam

    samtools depth -aa -H -r LP:1-120 ${prefix}.line_probe.sim.bam > ${prefix}.line_probe.depth.sim.tsv

    PROBE_COUNT=\$(samtools view -c -F 2304 \\
        -e 'rname=="LP" && rnext=="LP"' \\
        ${prefix}.line_probe.sim.bam)

    INSERT_COUNT=\$(samtools view -c -F 2304 \\
        -e 'rname!="LP" && rnext=="LP"' \\
        ${prefix}.line_probe.sim.bam)

    printf "sample\tbam_type\tcategory\tread_count\n" > ${prefix}.line_probe.counts.sim.tsv
    printf "${prefix}\tsim\tprobe\t\${PROBE_COUNT}\n"    >> ${prefix}.line_probe.counts.sim.tsv
    printf "${prefix}\tsim\tinsertion\t\${INSERT_COUNT}\n" >> ${prefix}.line_probe.counts.sim.tsv

    samtools view -F 2304 \\
        -e 'rname!="LP" && rnext=="LP"' \\
        -b ${prefix}.line_probe.sim.bam \\
        -o ${prefix}.insertion_reads.bam

    samtools index ${prefix}.insertion_reads.bam

    samtools view -F 2304 ${prefix}.insertion_reads.bam \\
        | awk 'OFS="\t" {
            start = (\$4 - 1) - 100
            if (start < 0) start = 0
            end = (\$4 - 1) + 100
            print \$3, start, end
        }' \\
        | sort -k1,1 -k2,2n \\
        | uniq > ${prefix}.insertion_sites.bed

    printf "site\tchr\tpos\tdepth\n" > ${prefix}.line_probe.insertion.depth.sim.tsv

    while IFS=\$'\t' read -r CHR START END; do
        SITE="\${CHR}:\$((START+1))-\${END}"
        REGION="\${CHR}:\$((START+1))-\${END}"
        samtools depth -aa -r "\$REGION" ${prefix}.insertion_reads.bam \\
            | awk -v site="\$SITE" 'OFS="\t" {print site, \$1, \$2, \$3}'
    done < ${prefix}.insertion_sites.bed >> ${prefix}.line_probe.insertion.depth.sim.tsv

    rm -f ${prefix}.insertion_reads.bam ${prefix}.insertion_reads.bam.bai ${prefix}.insertion_sites.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | tail -n 1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.line_probe.sim.bam
    touch ${prefix}.line_probe.sim.bam.bai
    touch ${prefix}.line_probe.depth.sim.tsv
    touch ${prefix}.line_probe.counts.sim.tsv
    touch ${prefix}.line_probe.insertion.depth.sim.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | tail -n 1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
