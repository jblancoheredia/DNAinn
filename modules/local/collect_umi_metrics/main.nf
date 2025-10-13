process COLLECT_UMI_METRICS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/collect_umi_consensus_metrics:1.0.3':
        'blancojmskcc/collect_umi_consensus_metrics:1.0.3' }"

    input:
    tuple val(meta) , path(cnsns_bam), path{cnsns_bai},
                      path(duplx_bam), path{duplx_bai},
                      path(smplx_bam), path{smplx_bai}

    output:
    tuple val(meta), path("*.dup_umi_counts.txt")         , emit: dup_umi_counts
    tuple val(meta), path("*.dup_family_sizes.txt")       , emit: dup_family_sizes
    tuple val(meta), path("*.duplex_umi_counts.txt")      , emit: duplex_umi_counts
    tuple val(meta), path("*.duplex_family_sizes.txt")    , emit: duplex_family_sizes
    tuple val(meta), path("*.duplex_yield_metrics.txt")   , emit: duplex_yield_metrics
    tuple val(meta), path("*.sim_umi_counts.txt")         , emit: sim_umi_counts
    tuple val(meta), path("*.sim_family_sizes.txt")       , emit: sim_family_sizes
    tuple val(meta), path("*.simplex_umi_counts.txt")     , emit: simplex_umi_counts
    tuple val(meta), path("*.simplex_family_sizes.txt")   , emit: simplex_family_sizes
    tuple val(meta), path("*.simplex_yield_metrics.txt")  , emit: simplex_yield_metrics
    tuple val(meta), path("*.cons_umi_counts.txt")        , emit: cons_umi_counts
    tuple val(meta), path("*.cons_family_sizes.txt")      , emit: cons_family_sizes
    tuple val(meta), path("*.consensus_umi_counts.txt")   , emit: consensus_umi_counts
    tuple val(meta), path("*.consensus_family_sizes.txt") , emit: consensus_family_sizes
    tuple val(meta), path("*.consensus_yield_metrics.txt"), emit: consensus_yield_metrics
    path  "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    export MPLCONFIGDIR=.

    mkdir duplex/
    mkdir simplex/
    mkdir consensus/

    collect_umi_metrics -i ${duplx_bam} -o duplex/${prefix}
    collect_umi_metrics -i ${smplx_bam} -o simplex/${prefix}
    collect_umi_metrics -i ${cnsns_bam} -o consensus/${prefix}

    mv duplex/${prefix}.umi_counts.txt ${prefix}.dup_umi_counts.txt
    mv duplex/${prefix}.family_sizes.txt ${prefix}.dup_family_sizes.txt
    mv duplex/${prefix}.duplex_umi_counts.txt ${prefix}.duplex_umi_counts.txt
    mv duplex/${prefix}.duplex_family_sizes.txt ${prefix}.duplex_family_sizes.txt
    mv duplex/${prefix}.duplex_yield_metrics.txt ${prefix}.duplex_yield_metrics.txt

    mv simplex/${prefix}.umi_counts.txt ${prefix}.sim_umi_counts.txt
    mv simplex/${prefix}.family_sizes.txt ${prefix}.sim_family_sizes.txt
    mv simplex/${prefix}.duplex_umi_counts.txt ${prefix}.simplex_umi_counts.txt
    mv simplex/${prefix}.duplex_family_sizes.txt ${prefix}.simplex_family_sizes.txt
    mv simplex/${prefix}.duplex_yield_metrics.txt ${prefix}.simplex_yield_metrics.txt

    mv consensus/${prefix}.umi_counts.txt ${prefix}.cons_umi_counts.txt
    mv consensus/${prefix}.family_sizes.txt ${prefix}.cons_family_sizes.txt
    mv consensus/${prefix}.duplex_umi_counts.txt ${prefix}.consensus_umi_counts.txt
    mv consensus/${prefix}.duplex_family_sizes.txt ${prefix}.consensus_family_sizes.txt
    mv consensus/${prefix}.duplex_yield_metrics.txt ${prefix}.consensus_yield_metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        collect_umi_metrics: 1.0.3
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.dup_umi_counts.txt
    touch ${prefix}.dup_family_sizes.txt
    touch ${prefix}.duplex_umi_counts.txt
    touch ${prefix}.duplex_family_sizes.txt
    touch ${prefix}.duplex_yield_metrics.txt
    touch ${prefix}.sim_umi_counts.txt
    touch ${prefix}.sim_family_sizes.txt
    touch ${prefix}.simplex_umi_counts.txt
    touch ${prefix}.simplex_family_sizes.txt
    touch ${prefix}.simplex_yield_metrics.txt
    touch ${prefix}.cons_umi_counts.txt
    touch ${prefix}.cons_family_sizes.txt
    touch ${prefix}.consensus_umi_counts.txt
    touch ${prefix}.consensus_family_sizes.txt
    touch ${prefix}.consensus_yield_metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        collect_umi_metrics: 1.0.3
    END_VERSIONS
    """
}
