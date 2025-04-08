process STANDARDIZE_SVS_ANNOTE {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge/label/cf202003::python"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/standardize_svs:1.4.0' :
        'blancojmskcc/standardize_svs:1.4.0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.summary.tsv")  , emit: tsv
    tuple val(meta), path("*.standard.vcf") , emit: vcf
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python /usr/local/bin/standardize_svs_annote.py ${vcf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        standardize_svs: "v1.4.0"
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.survivor.standard.vcf
    touch ${prefix}.survivor.summary.tsv 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        standardize_svs: "v1.4.0"
    END_VERSIONS
    """
}

