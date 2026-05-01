process MAPK_SNVINDEL {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/vcfcalls2tsv:2.0.0':
        'blancojmskcc/vcfcalls2tsv:2.0.0' }"

    input:
    tuple val(meta), path(variants_tsv)
    path mapk_genes_tsv

    output:
    tuple val(meta), path("${prefix}.mapk_snvindel.alterations.tsv")   , emit: alterations
    tuple val(meta), path("${prefix}.mapk_snvindel.gene_summary.tsv")  , emit: gene_summary
    tuple val(meta), path("${prefix}.mapk_snvindel.sample_summary.tsv"), emit: sample_summary
    path "versions.yml"                                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    python3 ${moduleDir}/resources/usr/bin/mapk_snvindel.py \\
        --sample-id "${meta.id}" \\
        --variants-tsv ${variants_tsv} \\
        --mapk-genes ${mapk_genes_tsv} \\
        --out-prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mapk_snvindel: \$(mapk_snvindel.py --version | sed 's/^mapk_snvindel //')
    END_VERSIONS
    """
}
