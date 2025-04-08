process MANTA_TUMORONLY {
    tag "$meta.id"
    label 'process_medium'
    label 'error_retry'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/manta:1.6.0--h9ee0642_1' :
        'quay.io/biocontainers/manta:1.6.0--h9ee0642_1' }"

    input:
    tuple val(meta),  path(input), path(input_index)
    tuple val(meta1), path(target_bed)
    tuple val(meta2), path(target_bed_tbi)
    tuple val(meta3), path(fasta)
    tuple val(meta4), path(fai)
    path(config)

    output:
    tuple val(meta), path("*candidate_small_indels.vcf.gz.tbi"), emit: candidate_small_indels_vcf_tbi
    tuple val(meta), path("*candidate_small_indels.vcf.gz")    , emit: candidate_small_indels_vcf
    tuple val(meta), path("*candidate_sv.vcf.gz.tbi")          , emit: candidate_sv_vcf_tbi
    tuple val(meta), path("*candidate_sv.vcf.gz")              , emit: candidate_sv_vcf
    tuple val(meta), path("*tumor_sv.vcf.gz.tbi")              , emit: tumor_sv_vcf_tbi
    tuple val(meta), path("*tumor_sv.vcf.gz")                  , emit: tumor_sv_vcf
    tuple val(meta), path("*.unfiltered.vcf")                  , emit: vcf
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def options_manta = target_bed ? "--callRegions $target_bed" : ""
    def config_option = config ? "--config ${config}" : ""
    """
    configManta.py \\
        --tumorBam $input \\
        --reference $fasta \\
        ${config_option} \\
        --runDir manta \\
        $options_manta \\
        $args

    python manta/runWorkflow.py -m local -j $task.cpus

    mv manta/results/variants/candidateSmallIndels.vcf.gz \\
        ${prefix}.candidate_small_indels.vcf.gz
    mv manta/results/variants/candidateSmallIndels.vcf.gz.tbi \\
        ${prefix}.candidate_small_indels.vcf.gz.tbi
    mv manta/results/variants/candidateSV.vcf.gz \\
        ${prefix}.candidate_sv.vcf.gz
    mv manta/results/variants/candidateSV.vcf.gz.tbi \\
        ${prefix}.candidate_sv.vcf.gz.tbi
    mv manta/results/variants/tumorSV.vcf.gz \\
        ${prefix}.tumor_sv.vcf.gz
    mv manta/results/variants/tumorSV.vcf.gz.tbi \\
        ${prefix}.tumor_sv.vcf.gz.tbi

    cp ${prefix}.tumor_sv.vcf.gz ${prefix}.unfiltered.vcf.gz
    gunzip ${prefix}.unfiltered.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        manta: \$( configManta.py --version )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.candidate_small_indels.vcf.gz
    touch ${prefix}.candidate_small_indels.vcf.gz.tbi
    echo "" | gzip > ${prefix}.candidate_sv.vcf.gz
    touch ${prefix}.candidate_sv.vcf.gz.tbi
    echo "" | gzip > ${prefix}.tumor_sv.vcf.gz
    touch ${prefix}.tumor_sv.vcf.gz.tbi
    touch ${prefix}.unfiltered.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        manta: \$( configManta.py --version )
    END_VERSIONS
    """
}
