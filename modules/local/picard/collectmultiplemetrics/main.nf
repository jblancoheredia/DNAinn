process PICARD_COLLECTMULTIPLEMETRICS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.3.0--hdfd78af_0' :
        'quay.io/biocontainers/picard:3.3.0--hdfd78af_0' }"

    input:
    tuple val(meta) , path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*_metrics"), emit: metrics
    tuple val(meta), path("*.pdf")    , emit: pdf, optional: true
    path  "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--REFERENCE_SEQUENCE ${fasta}" : ""
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard CollectMultipleMetrics] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    picard \\
        -Xmx${avail_mem}M \\
        CollectMultipleMetrics \\
        $args \\
        --INPUT $bam \\
        --OUTPUT ${prefix}.CollectMultipleMetrics \\
        $reference

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard CollectMultipleMetrics --version 2>&1 | grep -o 'Version.*' | cut -f2- -d:)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.CollectMultipleMetrics.alignment_summary_metrics
    touch ${prefix}.CollectMultipleMetrics.insert_size_metrics
    touch ${prefix}.CollectMultipleMetrics.quality_distribution.pdf
    touch ${prefix}.CollectMultipleMetrics.base_distribution_by_cycle_metrics
    touch ${prefix}.CollectMultipleMetrics.quality_by_cycle_metrics
    touch ${prefix}.CollectMultipleMetrics.read_length_histogram.pdf
    touch ${prefix}.CollectMultipleMetrics.base_distribution_by_cycle.pdf
    touch ${prefix}.CollectMultipleMetrics.quality_by_cycle.pdf
    touch ${prefix}.CollectMultipleMetrics.insert_size_histogram.pdf
    touch ${prefix}.CollectMultipleMetrics.quality_distribution_metrics

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard CollectMultipleMetrics --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}

process PICARD_COLLECTMULTIPLEMETRICS_RAW {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.3.0--hdfd78af_0' :
        'quay.io/biocontainers/picard:3.3.0--hdfd78af_0' }"

    input:
    tuple val(meta) , path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*_metrics"), emit: metrics
    tuple val(meta), path("*.pdf")    , emit: pdf, optional: true
    path  "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--REFERENCE_SEQUENCE ${fasta}" : ""
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard CollectMultipleMetrics] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    picard \\
        -Xmx${avail_mem}M \\
        CollectMultipleMetrics \\
        $args \\
        --INPUT $bam \\
        --OUTPUT ${prefix}.raw.CollectMultipleMetrics \\
        $reference

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard CollectMultipleMetrics --version 2>&1 | grep -o 'Version.*' | cut -f2- -d:)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.raw.CollectMultipleMetrics.alignment_summary_metrics
    touch ${prefix}.raw.CollectMultipleMetrics.insert_size_metrics
    touch ${prefix}.raw.CollectMultipleMetrics.quality_distribution.pdf
    touch ${prefix}.raw.CollectMultipleMetrics.base_distribution_by_cycle_metrics
    touch ${prefix}.raw.CollectMultipleMetrics.quality_by_cycle_metrics
    touch ${prefix}.raw.CollectMultipleMetrics.read_length_histogram.pdf
    touch ${prefix}.raw.CollectMultipleMetrics.base_distribution_by_cycle.pdf
    touch ${prefix}.raw.CollectMultipleMetrics.quality_by_cycle.pdf
    touch ${prefix}.raw.CollectMultipleMetrics.insert_size_histogram.pdf
    touch ${prefix}.raw.CollectMultipleMetrics.quality_distribution_metrics

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard CollectMultipleMetrics --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}

process PICARD_COLLECTMULTIPLEMETRICS_FIN {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.3.0--hdfd78af_0' :
        'quay.io/biocontainers/picard:3.3.0--hdfd78af_0' }"

    input:
    tuple val(meta) , path(bam)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*_metrics"), emit: metrics
    tuple val(meta), path("*.pdf")    , emit: pdf, optional: true
    path  "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--REFERENCE_SEQUENCE ${fasta}" : ""
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard CollectMultipleMetrics] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    picard \\
        -Xmx${avail_mem}M \\
        SortSam \\
        --INPUT ${bam} \\
        --OUTPUT ${prefix}.sorted.bam \\
        --SORT_ORDER coordinate \\
        --CREATE_INDEX \\
        --VALIDATION_STRINGENCY SILENT \\
        ${reference}

    picard \\
        -Xmx${avail_mem}M \\
        CollectMultipleMetrics \\
        --INPUT ${prefix}.sorted.bam \\
        --OUTPUT ${prefix}.fin.CollectMultipleMetrics \\
        ${reference}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard CollectMultipleMetrics --version 2>&1 | grep -o 'Version.*' | cut -f2- -d:)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fin.CollectMultipleMetrics.alignment_summary_metrics
    touch ${prefix}.fin.CollectMultipleMetrics.insert_size_metrics
    touch ${prefix}.fin.CollectMultipleMetrics.quality_distribution.pdf
    touch ${prefix}.fin.CollectMultipleMetrics.base_distribution_by_cycle_metrics
    touch ${prefix}.fin.CollectMultipleMetrics.quality_by_cycle_metrics
    touch ${prefix}.fin.CollectMultipleMetrics.read_length_histogram.pdf
    touch ${prefix}.fin.CollectMultipleMetrics.base_distribution_by_cycle.pdf
    touch ${prefix}.fin.CollectMultipleMetrics.quality_by_cycle.pdf
    touch ${prefix}.fin.CollectMultipleMetrics.insert_size_histogram.pdf
    touch ${prefix}.fin.CollectMultipleMetrics.quality_distribution_metrics

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard CollectMultipleMetrics --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
