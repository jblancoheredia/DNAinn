process COPYNCAT {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/copyncat_dnainn:1.0.0':
        'blancojmskcc/copyncat_dnainn:1.0.0' }"

    input:
    tuple val(meta), \
        path(cnvkit_call), \
        path(cnvkit_cns), \
        path(cnvkit_vcf), \
        path(controlfreec_cnv), \
        path(controlfreec_config), \
        path(sequenza_segments), \
        path(sequenza_confints), \
        path(sequenza_alternative), \
        path(oncocnv_profile), \
        path(facets_vcf)

    output:
    tuple val(meta), path("*_CNV_MERGED.tsv")       , emit: tsv
    tuple val(meta), path("*_CNV_PLOIDY_PURITY.tsv"), emit: summary
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def facets_vcf_decompressed = facets_vcf.name.endsWith('.gz') ? facets_vcf.name.replaceAll(/[.]gz$/, '') : facets_vcf.name
    """
    if [[ "${facets_vcf}" == *.gz ]]; then
        bgzip -d -f ${facets_vcf}
    fi

    CopyNcat \\
        --sample ${prefix} \\
        --cnvkit-call ${cnvkit_call} \\
        --cnvkit-cns ${cnvkit_cns} \\
        --cnvkit-vcf ${cnvkit_vcf} \\
        --controlfreec-cnv ${controlfreec_cnv} \\
        --controlfreec-config ${controlfreec_config} \\
        --sequenza-segments ${sequenza_segments} \\
        --sequenza-confints ${sequenza_confints} \\
        --sequenza-alternative ${sequenza_alternative} \\
        --oncocnv-profile ${oncocnv_profile} \\
        --facets-vcf ${facets_vcf_decompressed} \\
        --out-tsv ${prefix}_CNV_MERGED.tsv \\
        --out-summary ${prefix}_CNV_PLOIDY_PURITY.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        copyncat: \$(CopyNcat --version 2>&1 | sed -e "s/CopyNcat v//g")
    END_VERSIONS
    """
    stub:
    def args = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_CNV_MERGED.tsv
    touch ${prefix}_CNV_PLOIDY_PURITY.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        copyncat: \$(CopyNcat --version 2>&1 | sed -e "s/CopyNcat v//g")
    END_VERSIONS
    """
}
