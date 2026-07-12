process VARIANTMERGE_MAF {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/vcfcalls2tsv:2.2.0':
        'blancojmskcc/vcfcalls2tsv:2.2.0' }"

    input:
    tuple val(meta), path(tier_tsv)

    output:
    tuple val(meta), path("*.clinical_maf.tsv"), emit: clinical_maf_tsv
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def maf_label = task.ext.maf_label ?: 'high_confidence'
    """
    HC2MAF \\
        ${tier_tsv} \\
        --output ${prefix}.${maf_label}.clinical_maf.tsv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        high_confidence_to_clinical_maf: \$(echo \$(HC2MAF -v 2>&1 | grep "high_confidence_to_clinical_maf.py" | sed 's/^high_confidence_to_clinical_maf.py //'))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def maf_label = task.ext.maf_label ?: 'high_confidence'
    """
    touch ${prefix}.${maf_label}.clinical_maf.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        high_confidence_to_clinical_maf: \$(echo \$(HC2MAF -v 2>&1 | grep "high_confidence_to_clinical_maf.py" | sed 's/^high_confidence_to_clinical_maf.py //'))
    END_VERSIONS
    """
}
