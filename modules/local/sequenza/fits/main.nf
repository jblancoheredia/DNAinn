process SEQUENZA_FITS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/sequenza:3.0.1' :
        'blancojmskcc/sequenza:3.0.1' }"

    input:
    tuple val(meta),  path(seqz)

    output:
    tuple val(meta), path("*.tsv")  , emit: tsv
    tuple val(meta), path("*.pdf")  , emit: pdf
    tuple val(meta), path("*.RData"), emit: rdata
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """ 
    cat << EOF > run_sequenza.${prefix}.R
    Sys.setenv("VROOM_CONNECTION_SIZE" = 1000000000)

    library(sequenza)

    data.file <- "${seqz}"

    output.dir <- "sequenza_results"

    sample.id <- "${prefix}"

    seqz.data <- sequenza.extract(data.file)

    CP <- sequenza.fit(seqz.data)

    sequenza.results(sequenza.extract = seqz.data,
                 cp.table = CP,
                 sample.id = sample.id,
                 out.dir = output.dir)

    EOF

    Rscript run_sequenza.${prefix}.R

    mv sequenza_results/${prefix}* .

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
    touch ${prefix}.RData

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: \$(echo \$(sequenza-utils 2>&1) | sed 's/^.*is version //; s/ .*\$//')
    END_VERSIONS
    """
}
