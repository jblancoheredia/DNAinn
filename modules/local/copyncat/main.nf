process COPYNCAT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/copyncat_dnainn:1.0.0':
        'blancojmskcc/copyncat_dnainn:1.0.0' }"

    input:
    tuple val(meta), path(vcfs), path(tbis)

    output:
    tuple val(meta), path("*_CNV_MERGED.tsv")       , emit: tsv
    tuple val(meta), path("*_CNV_PLOIDY_PURITY.tsv"), emit: tsv
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    for f in *.vcf.gz; do
        bgzip -d \"\$f\"
    done
    
    rm *.vcf.gz.tbi

    CopyNcat \\
        --sample ${prefix} \\
        --cnvkit-call ${prefix}_sort.call.cns \\
        --cnvkit-cns ${prefix}_sort.cns \\
        --cnvkit-vcf ${prefix}.cnvcall.vcf \\
        --controlfreec-cnv ${prefix}.pileup.gz_CNVs \\
        --controlfreec-config config.txt \\
        --sequenza-segments ${prefix}_segments.txt \\
        --sequenza-confints ${prefix}_confints_CP.txt \\
        --sequenza-alternative ${prefix}_alternative_solutions.txt \\
        --oncocnv-profile ${prefix}_sort.profile.txt \\
        --facets-vcf ${prefix}.vcf.gz \\
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