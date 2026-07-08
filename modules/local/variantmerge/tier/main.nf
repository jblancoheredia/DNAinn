process VARIANTMERGE_TIER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/vcfcalls2tsv:2.1.0':
        'blancojmskcc/vcfcalls2tsv:2.1.0' }"

    input:
    tuple val(meta), path(merged_tsv)

    output:
    tuple val(meta), path("*.merged_variants.tiered.tsv")         , emit: tiered_tsv
    tuple val(meta), path("*.merged_variants.high_confidence.tsv"), emit: high_confidence_tsv
    tuple val(meta), path("*.merged_variants.review.tsv")         , emit: review_tsv
    tuple val(meta), path("*.merged_variants.likely_noise.tsv")   , emit: likely_noise_tsv
    path "versions.yml"                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    TIERSnMERGE \\
        ${merged_tsv} \\
        --output ${prefix}.merged_variants.tiered.tsv \\
        --high-confidence-output ${prefix}.merged_variants.high_confidence.tsv \\
        --review-output ${prefix}.merged_variants.review.tsv \\
        --noise-output ${prefix}.merged_variants.likely_noise.tsv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tier_merged_calls: \$(echo \$(TIERSnMERGE -v 2>&1 | grep "tier_merged_calls.py" | sed 's/^tier_merged_calls.py //'))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.merged_variants.tiered.tsv
    touch ${prefix}.merged_variants.high_confidence.tsv
    touch ${prefix}.merged_variants.review.tsv
    touch ${prefix}.merged_variants.likely_noise.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tier_merged_calls: \$(echo \$(TIERSnMERGE -v 2>&1 | grep "tier_merged_calls.py" | sed 's/^tier_merged_calls.py //'))
    END_VERSIONS
    """
}
