process SEQUENZA_FITS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/sequenza:3.0.0' :
        'blancojmskcc/sequenza:3.0.0' }"

    input:
    tuple val(meta),  path(seqz)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val(meta), path("*.pdf"), emit: pdf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """ 
    cat << EOF > run_sequenza.${meta.id}.R
    library(sequenza)

    Sys.setenv("VROOM_CONNECTION_SIZE" = 1000000000)

    ${meta.id} <- sequenza.extract("${seqz}", verbose = FALSE)

    CP <- sequenza.fit(${meta.id})

    sequenza.results(sequenza.extract = ${meta.id}, 
                    cp.table = CP, 
                    sample.id = "${meta.id}", 
                    out.dir = "${meta.id}.output")

    EOF

    Rscript run_sequenza.${meta.id}.R

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: \$(echo \$(sequenza-utils 2>&1) | sed 's/^.*is version //; s/ .*\$//')
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    touch ${prefix}.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: \$(echo \$(sequenza-utils 2>&1) | sed 's/^.*is version //; s/ .*\$//')
    END_VERSIONS
    """
}
