process BCFTOOLS_MPILEUP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'quay.io/biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta),  path(bam), path(bai)
    tuple val(meta1), path(fai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(intervals)

    output:
    tuple val(meta), path("*vcf.gz")     , emit: vcf
    tuple val(meta), path("*vcf.gz.tbi") , emit: tbi
    tuple val(meta), path("*stats.txt")  , emit: stats
    path  "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def intervals = intervals ? "--regions-file ${intervals}" : ""
    """
    echo "${meta.id}" > sample_name.list

    bcftools \\
        mpileup \\
        ${intervals} \\
        ${args} \\
        ${fasta} \\
        ${bam} | \\
        bcftools \\
        call \\
        ${args2} \\
        ${prefix}.pileup.vcf.tmp

    bcftools sort -o ${prefix}.pileup.vcf ${prefix}.pileup.vcf.tmp

    rm ${prefix}.pileup.vcf.tmp

    bgzip ${prefix}.pileup.vcf

    tabix ${prefix}.pileup.vcf.gz

    bcftools stats ${prefix}.pileup.vcf.gz > ${prefix}.pileup_bcftools_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pileup.vcf.gz
    touch ${prefix}.pileup.vcf.gz.tbi
    touch ${prefix}.pileup_bcftools_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
