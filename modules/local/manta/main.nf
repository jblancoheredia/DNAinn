process MANTA {
    tag "$meta.patient"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/manta:1.6.0--h9ee0642_1' :
        'quay.io/biocontainers/manta:1.6.0--h9ee0642_1' }"

    input:
    tuple val(meta) , path(tumour_bam), path(tumour_bai), path(normal_bam), path(normal_bai)
    tuple val(meta2), path(target_bed) 
    tuple val(meta3), path(bed_index)
    tuple val(meta4), path(fasta)
    tuple val(meta5), path(fai)
    path(config)

    output:
    tuple val(meta), path("*.candidate_small_indels.vcf.gz.tbi", optional: true) , emit: candidate_small_indels_vcf_tbi
    tuple val(meta), path("*.candidate_small_indels.vcf.gz", optional: true)     , emit: candidate_small_indels_vcf
    tuple val(meta), path("*.candidate_sv.vcf.gz.tbi", optional: true)           , emit: candidate_sv_vcf_tbi
    tuple val(meta), path("*.diploid_sv.vcf.gz.tbi", optional: true)             , emit: diploid_sv_vcf_tbi
    tuple val(meta), path("*.somatic_sv.vcf.gz.tbi", optional: true)             , emit: somatic_sv_vcf_tbi
    tuple val(meta), path("*.candidate_sv.vcf.gz", optional: true)               , emit: candidate_sv_vcf
    tuple val(meta), path("*.diploid_sv.vcf.gz", optional: true)                 , emit: diploid_sv_vcf
    tuple val(meta), path("*.somatic_sv.vcf.gz", optional: true)                 , emit: somatic_sv_vcf
    tuple val(meta), path("*.tsv", optional: true)                               , emit: metrics_tsv
    tuple val(meta), path("*.txt", optional: true)                               , emit: metrics_txt
    path "versions.yml"                                                          , emit: versions
    tuple val(meta), path("*.manta.unfiltered.vcf")                              , emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.patient}"
    def options_manta = target_bed ? "--callRegions $target_bed" : ""
    """
    configManta.py \\
        --tumorBam ${tumour_bam} \\
        --reference ${fasta} \\
        --runDir manta_tumour \\
        ${options_manta}    
    python manta_tumour/runWorkflow.py -m local -j ${task.cpus}

    mv manta_tumour/results/stats/svCandidateGenerationStats.tsv ${prefix}.manta_tumour.svCandidateGenerationStats.tsv
    mv manta_tumour/results/variants/candidateSmallIndels.vcf.gz.tbi ${prefix}.candidate_small_indels.vcf.gz.tbi
    mv manta_tumour/results/stats/alignmentStatsSummary.txt ${prefix}.manta_tumour.alignmentStatsSummary.txt
    mv manta_tumour/results/variants/candidateSmallIndels.vcf.gz ${prefix}.candidate_small_indels.vcf.gz
    mv manta_tumour/results/stats/svLocusGraphStats.tsv ${prefix}.manta_tumour.svLocusGraphStats.tsv
    mv manta_tumour/results/variants/candidateSV.vcf.gz.tbi ${prefix}.manta.unfiltered.vcf.gz.tbi
    mv manta_tumour/results/variants/tumorSV.vcf.gz.tbi ${prefix}.tumour_only_sv.vcf.gz.tbi
    mv manta_tumour/results/variants/candidateSV.vcf.gz ${prefix}.manta.unfiltered.vcf.gz
    mv manta_tumour/results/variants/tumorSV.vcf.gz ${prefix}.tumour_only.sv.vcf.gz

    zcat ${prefix}.tumour_only.sv.vcf.gz | grep -v "#" | gzip >> ${prefix}.manta.unfiltered.vcf.gz

    gunzip ${prefix}.manta.unfiltered.vcf.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        manta: \$( configManta.py --version )
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.patient}"
    """
    touch *.manta.unfiltered.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        manta: \$( configManta.py --version )
    END_VERSIONS
    """
}
