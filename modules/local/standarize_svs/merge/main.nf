process STANDARDIZE_SVS_MERGE {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge/label/cf202003::python"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/standardize_svs:1.4.0' :
        'blancojmskcc/standardize_svs:1.4.0' }"

    input:
    tuple val(meta) , path(delly_vcf)
    tuple val(meta1), path(svaba_vcf)
    tuple val(meta2), path(manta_vcf)
    tuple val(meta3), path(tiddit_vcf)
    tuple val(meta4), path(gridss_vcf)

    output:
    tuple val(meta), path("*.delly.standard.vcf")   , emit: delly_vcf
    tuple val(meta), path("*.delly.summary.tsv")    , emit: delly_tsv
    tuple val(meta), path("*.svaba.standard.vcf")   , emit: svaba_vcf
    tuple val(meta), path("*.svaba.summary.tsv")    , emit: svaba_tsv
    tuple val(meta), path("*.manta.standard.vcf")   , emit: manta_vcf
    tuple val(meta), path("*.manta.summary.tsv")    , emit: manta_tsv
    tuple val(meta), path("*.tiddit.standard.vcf")  , emit: tiddit_vcf
    tuple val(meta), path("*.tiddit.summary.tsv")   , emit: tiddit_tsv
    tuple val(meta), path("*.gridss.standard.vcf")  , emit: gridss_vcf
    tuple val(meta), path("*.gridss.summary.tsv")   , emit: gridss_tsv
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python /usr/local/bin/standardize_svs_merge.py \$(ls *.vcf)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        standardize_svs: "v1.4.0"
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.delly.standard.vcf
    touch ${prefix}.delly.summary.tsv 
    touch ${prefix}.svaba.standard.vcf
    touch ${prefix}.svaba.summary.tsv 
    touch ${prefix}.manta.standard.vcf
    touch ${prefix}.manta.summary.tsv 
    touch ${prefix}.tiddit.standard.vcf
    touch ${prefix}.tiddit.summary.tsv
    touch ${prefix}.gridss.standard.vcf
    touch ${prefix}.gridss.summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        standardize_svs: "v1.4.0"
    END_VERSIONS
    """
}

