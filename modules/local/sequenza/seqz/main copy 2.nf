process SEQUENZA_SEQZ {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/sequenza:3.0.0' :
        'blancojmskcc/sequenza:3.0.0' }"

    input:
    tuple val(meta),  path(tbam), path(tbai)
    tuple val(meta1), path(fasta)
    tuple val(meta2), path(fastai)
    path (nbam)
    path (nbai)
    path (wigfile)

    output:
    tuple val(meta), path("*_bin50.seqz.gz"), emit: seqz
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
        -f "illumina" \\
        -gc ${wigfile} \\
        -o ${prefix}.seqz.gz \\
        $args

    sequenza-utils \\
        seqz_binning \\
        --seqz ${prefix}.seqz.gz \\
        --window 50 \\
        -o ${prefix}_bin50.seqz.gz

    cat << EOF > run_sequenza.${meta.id}.R
    library(sequenza)

    ${meta.id} <- sequenza.extract("${meta.id}_bin50.seqz.gz", verbose = FALSE)

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
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_bin50.seqz.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: \$(echo \$(sequenza-utils 2>&1) | sed 's/^.*is version //; s/ .*\$//')
    END_VERSIONS
    """
}
