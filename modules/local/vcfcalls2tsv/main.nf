process VCFCALLS2TSV {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/vcfcalls2tsv:2.1.0':
        'blancojmskcc/vcfcalls2tsv:2.1.0' }"

    input:
    tuple val(meta), path(vcfs), path(tbis)

    output:
    tuple val(meta), path("*.merged_variants.tsv")                , emit: tsv
    tuple val(meta), path("*.merged_variants.review.tsv")         , emit: review_tsv
    tuple val(meta), path("*.merged_variants.tiered.tsv")         , emit: tiered_tsv
    tuple val(meta), path("*.high_confidence.clinical_maf.tsv")   , emit: clinical_maf_tsv
    tuple val(meta), path("*.merged_variants.likely_noise.tsv")   , emit: likely_noise_tsv
    tuple val(meta), path("*.merged_variants.high_confidence.tsv"), emit: high_confidence_tsv
    path "versions.yml"                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    VCF2TSV \\
        ${prefix}.freebayes.ann.vcf  \\
        ${prefix}.lofreq.somatic_final.snvs.ann.vcf \\
        ${prefix}.lofreq.somatic_final.indels.ann.vcf \\
        ${prefix}.mutect2.filtered.ann.vcf \\
        ${prefix}.vardict.ann.vcf \\
        --output ${prefix}.merged_variants.tsv

    TIERSnMERGE \\
        ${prefix}.merged_variants.tsv \\
        --output ${prefix}.merged_variants.tiered.tsv \\
        --high-confidence-output ${prefix}.merged_variants.high_confidence.tsv \\
        --review-output ${prefix}.merged_variants.review.tsv \\
        --noise-output ${prefix}.merged_variants.likely_noise.tsv

    HC2MAF \\
        ${prefix}.merged_variants.high_confidence.tsv \\
        --output ${prefix}.high_confidence.clinical_maf.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcfcalls2tsv: \$(echo \$(VCF2TSV -v 2>&1 | grep "vcfcalls2tsv.py" | sed 's/^vcfcalls2tsv.py //'))
        tier_merged_calls: \$(echo \$(TIERSnMERGE -v 2>&1 | grep "tier_merged_calls.py" | sed 's/^tier_merged_calls.py //'))
        high_confidence_to_clinical_maf: \$(echo \$(HC2MAF -v 2>&1 | grep "high_confidence_to_clinical_maf.py" | sed 's/^high_confidence_to_clinical_maf.py //'))
    END_VERSIONS
    """
    stub:
    def args = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.merged_variants.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcfcalls2tsv: \$(echo \$(VCF2TSV -v 2>&1 | grep "vcfcalls2tsv.py" | sed 's/^vcfcalls2tsv.py //'))
        tier_merged_calls: \$(echo \$(TIERSnMERGE -v 2>&1 | grep "tier_merged_calls.py" | sed 's/^tier_merged_calls.py //'))
        high_confidence_to_clinical_maf: \$(echo \$(HC2MAF -v 2>&1 | grep "high_confidence_to_clinical_maf.py" | sed 's/^high_confidence_to_clinical_maf.py //'))    
    END_VERSIONS
    """
}