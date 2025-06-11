process POLYSOLVER_CALLHLATYPE {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::gatk=3.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/polysolver:4.0.0' :
        'blancojmskcc/polysolver:4.0.0' }"

    input:
    tuple val(meta) , path(bam)
    tuple val(meta2), path(bai)
    val ethnicity

    output:
    tuple val(meta), path("*winners.hla.txt"), emit: hla
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    shell_call_hla_type \\
        ${bam} \\
        ${ethnicity} \\
        1 \\
        hg19 \\
        STDFQ \\
        0 \\
        ${prefix} \\
        $args 

    cp ${prefix}/winners.hla.txt ${prefix}_polysolver.winners.hla.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        polysolver: 4.0.0
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    
    touch ${prefix}_polysolver.winners.hla.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        polysolver: 4.0.0
    END_VERSIONS
    """
}
