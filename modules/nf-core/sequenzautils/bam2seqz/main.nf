process SEQUENZAUTILS_BAM2SEQZ {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/sequenza:3.0.0' :
        'blancojmskcc/sequenza:3.0.0' }"

    input:
    tuple val(meta) , path(tumourbam), path(tumourbai)
    tuple val(meta3), path(fasta)
    tuple val(meta0), path(fastai)
    tuple val(meta1), path(normalbam)
    tuple val(meta2), path(normalbai)
    path(wigfile)

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
        $args \\
        -n $normalbam \\
        -t $tumourbam \\
        --fasta $fasta \\
        -gc $wigfile \\
        -o ${prefix}.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: \$(echo \$(sequenza-utils 2>&1) | sed 's/^.*is version //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: \$(echo \$(sequenza-utils 2>&1) | sed 's/^.*is version //; s/ .*\$//')
    END_VERSIONS
    """
}
