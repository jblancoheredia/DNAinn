process INSURVEYOR {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/insurveyor:1.1.2' :
        'blancojmskcc/insurveyor:1.1.2' }"

    input:
    tuple val(meta),  path(bam), path(bai)
    tuple val(meta1), path(fasta)
    tuple val(meta2), path(fastai)

    output:
    tuple val(meta), path("*_insurveyor.vcf")   , emit: vcf
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    insurveyor.py \\
        --threads ${task.cpus} \\
        ${bam} \\
        . \\
        ${fasta}

    mv out.pass.vcf.gz ${prefix}_insurveyor.vcf.gz

    gunzip  ${prefix}_insurveyor.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        insurveyor: \$(echo \$(insurveyor.py --version 2>&1) | sed 's/^.*INSurVeyor v//; s/ .*\$//')
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_insurveyor.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        insurveyor: \$(echo \$(insurveyor.py --version 2>&1) | sed 's/^.*INSurVeyor v//; s/ .*\$//')
    END_VERSIONS
    """
}
