process PRESEQ_GCEXTRAP {
    tag "$meta.id"
    label 'process_high'
    label 'error_retry'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/preseq:3.2.0--hdcf5f25_6':
        'biocontainers/preseq:3.2.0--hdcf5f25_6' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.gc_extrap.txt"), emit: gc_extrap
    tuple val(meta), path("*.log")          , emit: log
    path  "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    args = task.attempt > 1 ? args.join(' -defects') : args  // Disable testing for defects
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paired_end = meta.single_end ? '' : '-pe'
    """
    preseq \\
        gc_extrap \\
        -B \\
        -output ${prefix}.gc_extrap.txt \\
        ${bam}
    cp .command.err ${prefix}.command.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        preseq: \$(echo \$(preseq 2>&1) | sed 's/^.*Version: //; s/Usage:.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gc_extrap.txt
    touch ${prefix}.command.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        preseq: \$(echo \$(preseq 2>&1) | sed 's/^.*Version: //; s/Usage:.*\$//')
    END_VERSIONS
    """
}
