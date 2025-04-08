process BCFTOOLS_MERGE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'quay.io/biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta) , path(vcfs), path(tbis)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    path(bed)

    output:
    tuple val(meta), path("*.merged.vcf.gz"),       emit: vcf
    tuple val(meta), path("*.merged.vcf.gz.tbi"),   emit: index, optional: true
    path "versions.yml",                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools norm \\
        -f ${fasta} \\
        -m \\
        -any \\
        ${prefix}.vardict.vcf.gz \\
        -o ${prefix}.vardict.norm.vcf.gz

    tabix ${prefix}.vardict.norm.vcf.gz

    bcftools norm \\
        -f ${fasta} \\
        -m \\
        -any \\
        ${prefix}.lofreq.vcf.gz \\ \\
        -o ${prefix}.lofreq.norm.vcf.gz

    tabix ${prefix}.lofreq.norm.vcf.gz

    bcftools norm \\
        -f ${fasta} \\
        -m \\
        -any \\
        ${prefix}.mutect2.filtered.vcf.gz \\
        -o ${prefix}.mutect2.norm.vcf.gz

    tabix ${prefix}.mutect2.norm.vcf.gz

    bcftools norm \\
        -f ${fasta} \\
        -m \\
        -any \\
        ${prefix}.freebayes.vcf.gz \\
        -o ${prefix}.freebayes.norm.vcf.gz

    tabix ${prefix}.freebayes.norm.vcf.gz

    bcftools merge \\
        --regions-file ${bed} \\
        --threads ${task.cpus} \\
        --output ${prefix}.merged.vcf \\
        --force-samples \\
        --info-rules TYPE:join,AF:max \\
        ${prefix}.lofreq.norm.vcf.gz \\
        ${prefix}.mutect2.norm.vcf.gz \\
        ${prefix}.vardict.norm.vcf.gz  \\
        ${prefix}.freebayes.norm.vcf.gz \\

    bgzip ${prefix}.merged.vcf

    tabix ${prefix}.merged.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
    stub:
    def args = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.merged.vcf.gz
    touch ${prefix}.merged.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
