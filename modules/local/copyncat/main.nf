process VCFCALLS2TSV {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/vcfcalls2tsv:1.0.0':
        'blancojmskcc/vcfcalls2tsv:1.0.0' }"

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
    for f in *.vcf.gz; do
        bgzip -d \"\$f\"
    done
    
    rm *.vcf.gz.tbi

    vcfcalls2tsv.py \\
        ${prefix}.freebayes.vcf  \\
        ${prefix}.lofreq.somatic_final.snvs.vcf \\
        ${prefix}.lofreq.somatic_final.indels.vcf \\
        ${prefix}.mutect2.filtered.vcf \\
        ${prefix}.vardict.vcf \\
        --output ${prefix}.merged_variants.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcfcalls2tsv: "1.0.0"
    END_VERSIONS
    """
    stub:
    def args = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.merged_variants.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcfcalls2tsv: "1.0.0"
    END_VERSIONS
    """
}