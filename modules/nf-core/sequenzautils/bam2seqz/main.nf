process SEQUENZAUTILS_BAM2SEQZ {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sequenza-utils:3.0.0--py38h6ed170a_2' :
        'quay.io/biocontainers/sequenza-utils:3.0.0--py38h6ed170a_2' }"

    input:
    tuple val(meta),  path(tbam), path(tbai)
    tuple val(meta1), path(fasta)
    tuple val(meta2), path(fastai)
    path (nbam)
    path (nbai)
    path (wigfile)

    output:
    tuple val(meta), path("*.gz"), emit: seqz
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sequenza-utils \\
        bam2seqz \\
        --normal ${tbam} \\
        --tumor  ${tbam} \\
        --normal2 ${nbam} \\
        --fasta ${fasta} \\
        -gc ${wigfile} \\
        -o ${prefix}.gz \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: \$(echo \$(sequenza-utils 2>&1) | sed 's/^.*is version //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: \$(echo \$(sequenza-utils 2>&1) | sed 's/^.*is version //; s/ .*\$//')
    END_VERSIONS
    """
}
