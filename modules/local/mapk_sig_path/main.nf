process MAPK_SIG_PATH {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/vcfcalls2tsv:2.0.0':
        'blancojmskcc/vcfcalls2tsv:2.0.0' }"

    input:
    tuple val(meta), path(variants_tsv), path(cnv_tsv), path(sv_tsv)
    path mapk_genes_tsv

    output:
    tuple val(meta), path("${prefix}.mapk_sig_path.alterations.tsv")   , emit: alterations
    tuple val(meta), path("${prefix}.mapk_sig_path.gene_summary.tsv")  , emit: gene_summary
    tuple val(meta), path("${prefix}.mapk_sig_path.sample_summary.tsv"), emit: sample_summary
    path "versions.yml"                                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    python3 ${moduleDir}/resources/usr/bin/mapk_sig_path.py \\
        --sv-tsv ${sv_tsv} \\
        --cnv-tsv ${cnv_tsv} \\
        --out-prefix ${prefix} \\
        --sample-id "${meta.id}" \\
        --mapk-genes ${mapk_genes_tsv} \\
        --variants-tsv ${variants_tsv}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mapk_sig_path: \$(python3 ${moduleDir}/resources/usr/bin/mapk_sig_path.py --version | sed 's/^mapk_sig_path //')
    END_VERSIONS
    """
}
