process VCFCALLS2TSV {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/vcfcalls2tsv:2.0.0':
        'blancojmskcc/vcfcalls2tsv:2.0.0' }"

    input:
    tuple val(meta), path(vcfs), path(tbis)

    output:
    tuple val(meta), path("*.merged_variants.tsv")  , emit: tsv
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    vcfcalls2tsv.py \\
        ${prefix}.freebayes.vcf.gz  \\
        ${prefix}.lofreq.somatic_final.snvs.vcf.gz \\
        ${prefix}.lofreq.somatic_final.indels.vcf.gz \\
        ${prefix}.mutect2.filtered.vcf.gz \\
        ${prefix}.vardict.vcf.gz \\
        --output ${prefix}.merged_variants.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcfcalls2tsv: \$(echo \$(vcfcalls2tsv.py -v 2>&1 | grep "vcfcalls2tsv.py" | sed 's/^vcfcalls2tsv.py //'))
    END_VERSIONS
    """
    stub:
    def args = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.merged_variants.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcfcalls2tsv: \$(echo \$(vcfcalls2tsv.py -v 2>&1 | grep "vcfcalls2tsv.py" | sed 's/^vcfcalls2tsv.py //'))
    END_VERSIONS
    """
}