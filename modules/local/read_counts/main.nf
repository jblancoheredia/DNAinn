process UMI_READ_COUNTS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'quay.io/biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta) , path(unmap_bam),
                      path(umifx_bam),
                      path(group_bam),
                      path(calls_bam),
                      path(filtr_bam), path(filtr_bai),
                      path(cnsns_bam), path(cnsns_bai),
                      path(duplx_bam), path(duplx_bai),
                      path(smplx_bam), path(smplx_bai)

    output:
    tuple val(meta), path("*.total_umi_counts.tsv"), emit: tsv
    path  "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    total_raw=\$(samtools view -c "${unmap_bam}")
    total_fix=\$(samtools view -c "${umifx_bam}")
    total_group=\$(samtools view -c "${group_bam}")
    total_calls=\$(samtools view -c "${calls_bam}")
    total_filtr=\$(samtools view -c "${filtr_bam}")
    total_cnsns=\$(samtools view -c "${cnsns_bam}")
    total_duplex=\$(samtools view -c -f 0x2 "${duplx_bam}")
    total_simplex=\$(samtools view -c -f 0x2 "${smplx_bam}")
    echo -e "Sample\tTotal_Raw_Reads\tTotal_UMI_Fixed_Reads\tTotal_UMI_Grouped_Reads\tTotal_UMI_Called\tTotal_UMI_Filtered\tTotal_UMI_Consensus\tTotal_UMI_Duplex\tTotal_UMI_Simplex" > ${prefix}.total_umi_counts.tsv
    echo -e "${prefix}\t\${total_raw}\t\${total_fix}\t\${total_group}\t\${total_calls}\t\${total_filtr}\t\${total_cnsns}\t\${total_duplex}\t\${total_simplex}" >> ${prefix}.total_umi_counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.total_umi_counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}

process DR_READ_COUNTS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'quay.io/biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta) , path(trimm_bam),
                      path(dedup_bam),
                      path(recal_bam)

    output:
    tuple val(meta), path("*.total_dr_counts.tsv"), emit: tsv
    path  "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    total_trimm=\$(samtools view -c "${trimm_bam}")
    total_dedup=\$(samtools view -c "${dedup_bam}")
    total_recal=\$(samtools view -c "${recal_bam}")
    echo -e "Sample\tTotal_Trimmed_Reads\tTotal_DeDuplicated_Reads\tTotal_ReCallibrated_Reads" > ${prefix}.total_dr_counts.tsv
    echo -e "${prefix}\t\${total_trimm}\t\${total_dedup}\t\${total_recal}" >> ${prefix}.total_dr_counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.total_dr_counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
