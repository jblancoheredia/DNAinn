process STANDARDIZE_SVS_FILTER {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge/label/cf202003::python"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/standardize_svs:1.4.0' :
        'blancojmskcc/standardize_svs:1.4.0' }"

    input:
    tuple val(meta) , path(recall_vcf)

    output:
    tuple val(meta), path("*.recall.standard.vcf")  , emit: recall_vcf
    tuple val(meta), path("*.recall.summary.tsv")   , emit: recall_tsv
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python /usr/local/bin/standardize_svs_filter.py ${recall_vcf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        standardize_svs: "v1.4.0"
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.recall.standard.vcf
    touch ${prefix}.recall.summary.tsv 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        standardize_svs: "v1.4.0"
    END_VERSIONS
    """
}

