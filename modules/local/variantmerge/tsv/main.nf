process VARIANTMERGE_TSV {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/vcfcalls2tsv:2.1.0':
        'blancojmskcc/vcfcalls2tsv:2.1.0' }"

    input:
    tuple val(meta), path(vcfs), path(tbis)

    output:
    tuple val(meta), path("*.merged_variants.tsv"), emit: tsv
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    VCF2TSV \\
        ${vcfs} \\
        --output ${prefix}.merged_variants.tsv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcfcalls2tsv: \$(echo \$(VCF2TSV -v 2>&1 | grep "vcfcalls2tsv.py" | sed 's/^vcfcalls2tsv.py //'))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.merged_variants.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcfcalls2tsv: \$(echo \$(VCF2TSV -v 2>&1 | grep "vcfcalls2tsv.py" | sed 's/^vcfcalls2tsv.py //'))
    END_VERSIONS
    """
}
