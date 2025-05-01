process GETBASECOUNTS_MULTISAMPLE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/getbasecounts:1.4.0':
        'blancojmskcc/getbasecounts:1.4.0' }"

    input:
    tuple val(meta) , path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(vcf)

    output:
    tuple val(meta), path("*.tvs"), emit: tvs
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    GetBaseCountsMultiSample \\
        --fasta ${fasta} \\
        --bam ${bam} \\
        --vcf ${vcf} \\
        --thread ${task.cpus} \\
        --output ${prefix}_getbasecounts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getbasecounts: \$(GetBaseCountsMultiSample |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_getbasecounts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getbasecounts: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
