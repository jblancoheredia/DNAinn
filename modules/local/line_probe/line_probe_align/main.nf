/*
 * Realign reads to hg38+LP (BWA-MEM2) and index.
 */

process LINE_PROBE_ALIGN_RAW {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::fgbio=2.4.0 bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/umi_aligner:1.0.0' :
        'blancojmskcc/umi_aligner:1.0.0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path probe_fasta
    path bwa_dir

    output:
    tuple val(meta), path("*.line_probe.raw.bam")    , emit: bam
    tuple val(meta), path("*.line_probe.raw.bam.bai"), emit: bai
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-B 3 -K 10000000 -a -Y -M'
    def prefix = task.ext.prefix ?: "${meta.id}"

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

    bwa-mem2 mem \\
        -t ${task.cpus} \\
        ${args} \\
        -R ${meta.read_group} \\
        \$BWA_INDEX_PREFIX \\
        ${prefix}.R1.fastq \\
        ${prefix}.R2.fastq \\
        | samtools sort -@ ${task.cpus} -o ${prefix}.line_probe.raw.bam

    samtools index ${prefix}.line_probe.raw.bam

    cat <<END_VERSIONS > versions.yml
"${task.process}":
    samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    bwa: \$(bwa 2>&1 | grep Version | sed 's/Version: //')
END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.line_probe.raw.bam
    touch ${prefix}.line_probe.raw.bam.bai

    cat <<END_VERSIONS > versions.yml
"${task.process}":
    samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    bwa: \$(bwa 2>&1 | grep Version | sed 's/Version: //')
END_VERSIONS
    """
}

process LINE_PROBE_ALIGN_CON {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::fgbio=2.4.0 bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/umi_aligner:1.0.0' :
        'blancojmskcc/umi_aligner:1.0.0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path probe_fasta
    path bwa_dir

    output:
    tuple val(meta), path("*.line_probe.con.bam")    , emit: bam
    tuple val(meta), path("*.line_probe.con.bam.bai"), emit: bai
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-B 3 -K 10000000 -a -Y -M'
    def prefix = task.ext.prefix ?: "${meta.id}"

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

    bwa-mem2 mem \\
        -t ${task.cpus} \\
        ${args} \\
        -R ${meta.read_group} \\
        \$BWA_INDEX_PREFIX \\
        ${prefix}.R1.fastq \\
        ${prefix}.R2.fastq \\
        | samtools sort -@ ${task.cpus} -o ${prefix}.line_probe.con.bam

    samtools index ${prefix}.line_probe.con.bam

    cat <<END_VERSIONS > versions.yml
"${task.process}":
    samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    bwa: \$(bwa 2>&1 | grep Version | sed 's/Version: //')
END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.line_probe.con.bam
    touch ${prefix}.line_probe.con.bam.bai

    cat <<END_VERSIONS > versions.yml
"${task.process}":
    samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    bwa: \$(bwa 2>&1 | grep Version | sed 's/Version: //')
END_VERSIONS
    """
}

process LINE_PROBE_ALIGN_DUP {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::fgbio=2.4.0 bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/umi_aligner:1.0.0' :
        'blancojmskcc/umi_aligner:1.0.0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path probe_fasta
    path bwa_dir

    output:
    tuple val(meta), path("*.line_probe.dup.bam")    , emit: bam
    tuple val(meta), path("*.line_probe.dup.bam.bai"), emit: bai
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-B 3 -K 10000000 -a -Y -M'
    def prefix = task.ext.prefix ?: "${meta.id}"

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

    bwa-mem2 mem \\
        -t ${task.cpus} \\
        ${args} \\
        -R ${meta.read_group} \\
        \$BWA_INDEX_PREFIX \\
        ${prefix}.R1.fastq \\
        ${prefix}.R2.fastq \\
        | samtools sort -@ ${task.cpus} -o ${prefix}.line_probe.dup.bam

    samtools index ${prefix}.line_probe.dup.bam

    cat <<END_VERSIONS > versions.yml
"${task.process}":
    samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    bwa: \$(bwa 2>&1 | grep Version | sed 's/Version: //')
END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.line_probe.dup.bam
    touch ${prefix}.line_probe.dup.bam.bai

    cat <<END_VERSIONS > versions.yml
"${task.process}":
    samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    bwa: \$(bwa 2>&1 | grep Version | sed 's/Version: //')
END_VERSIONS
    """
}

process LINE_PROBE_ALIGN_SIM {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::fgbio=2.4.0 bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/umi_aligner:1.0.0' :
        'blancojmskcc/umi_aligner:1.0.0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path probe_fasta
    path bwa_dir

    output:
    tuple val(meta), path("*.line_probe.sim.bam")    , emit: bam
    tuple val(meta), path("*.line_probe.sim.bam.bai"), emit: bai
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-B 3 -K 10000000 -a -Y -M'
    def prefix = task.ext.prefix ?: "${meta.id}"

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

    bwa-mem2 mem \\
        -t ${task.cpus} \\
        ${args} \\
        -R ${meta.read_group} \\
        \$BWA_INDEX_PREFIX \\
        ${prefix}.R1.fastq \\
        ${prefix}.R2.fastq \\
        | samtools sort -@ ${task.cpus} -o ${prefix}.line_probe.sim.bam

    samtools index ${prefix}.line_probe.sim.bam

    cat <<END_VERSIONS > versions.yml
"${task.process}":
    samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    bwa: \$(bwa 2>&1 | grep Version | sed 's/Version: //')
END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.line_probe.sim.bam
    touch ${prefix}.line_probe.sim.bam.bai

    cat <<END_VERSIONS > versions.yml
"${task.process}":
    samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    bwa: \$(bwa 2>&1 | grep Version | sed 's/Version: //')
END_VERSIONS
    """
}
