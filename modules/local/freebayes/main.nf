process FREEBAYES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/freebayes:1.3.6--hbfe0e7f_2' :
        'quay.io/biocontainers/freebayes:1.3.6--hbfe0e7f_2' }"

    input:
    tuple val(meta),  path(bam), path(bai), path(fasta), path(fasta_fai)
    tuple val(meta3), path(fasta)
    tuple val(meta4), path(fasta_fai)
    path(target_bed)

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"),    emit: vcf
    path  "versions.yml",                                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat ${target_bed} | awk '{ print \$1 \"\t\" \$2 \"\t\" \$3 }' > tmp.bed

    freebayes \\
        -f ${fasta} \\
        --target tmp.bed \\
        ${args} \\
        ${bam} > ${prefix}.vcf

    cat ${prefix}.vcf | awk '\$1 ~ /^#/ {print \$0;next} {print \$0 | \"sort -k1,1 -k2,2n\"}' > ${prefix}.freebayes.vcf

    bgzip ${prefix}.freebayes.vcf

    tabix ${prefix}.freebayes.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freebayes: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g' )
    END_VERSIONS
    """
}
