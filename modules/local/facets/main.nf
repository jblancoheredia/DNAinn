process FACETS_CNV {
    tag "$meta.patient"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/cnv_facets:0.16.1' :
        'blancojmskcc/cnv_facets:0.16.1' }"

    input:
    tuple val(meta),  path(tbam), path(tbai), path(nbam), path(nbai)
    path(common_vcf)
    path(common_vcf_tbi)

    output:
    tuple val(meta), path("*.csv.gz"),      emit: pileup
    tuple val(meta), path("*.vcf.gz"),      emit: vcf
    tuple val(meta), path("*.cnv.png"),     emit: png
    tuple val(meta), path("*.cov.pdf"),     emit: coverage
    tuple val(meta), path("*.spider.pdf"),  emit: spider
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.patient}"
    """
    cnv_facets.R \\
        -t ${tbam} \\
        -n ${nbam} \\
        -o ${prefix} \\
        -vcf ${common_vcf} \\
        --snp-nprocs ${task.cpus} \\
        $args 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        facets: \$(cnv_facets.R -v |& grep -o 'facets=[0-9.]*' | cut -d= -f2)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.patient}"
    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.cnv.png
    touch ${prefix}.cov.pdf
    touch ${prefix}.spider.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        facets: \$(cnv_facets.R -v |& grep -o 'facets=[0-9.]*' | cut -d= -f2)
    END_VERSIONS
    """
}
