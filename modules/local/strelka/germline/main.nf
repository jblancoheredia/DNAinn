process STRELKA_GERMLINE {
    tag "$meta.id"
    label 'process_medium'
    label 'error_retry'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/strelka:2.9.10--h9ee0642_1' :
        'quay.io/biocontainers/strelka:2.9.10--h9ee0642_1' }"

    input:
    tuple val(meta),  path(input), path(input_index)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    path(target_bed)
    path(target_bed_idx)

    output:
    tuple val(meta), path("*variants.vcf.gz")    , emit: vcf
    tuple val(meta), path("*variants.vcf.gz.tbi"), emit: vcf_tbi
    tuple val(meta), path("*genome.vcf.gz")      , emit: genome_vcf
    tuple val(meta), path("*genome.vcf.gz.tbi")  , emit: genome_vcf_tbi
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    configureStrelkaGermlineWorkflow.py \\
        --bam $input \\
        --referenceFasta $fasta \\
        --callRegions ${target_bed} \\
        $args \\
        --runDir strelka

    sed -i s/"isEmail = isLocalSmtp()"/"isEmail = False"/g strelka/runWorkflow.py

    python strelka/runWorkflow.py -m local -j $task.cpus
    mv strelka/results/variants/genome.*.vcf.gz     ${prefix}.genome.vcf.gz
    mv strelka/results/variants/genome.*.vcf.gz.tbi ${prefix}.genome.vcf.gz.tbi
    mv strelka/results/variants/variants.vcf.gz     ${prefix}.variants.vcf.gz
    mv strelka/results/variants/variants.vcf.gz.tbi ${prefix}.variants.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        strelka: \$( configureStrelkaGermlineWorkflow.py --version )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.genome.vcf.gz
    touch ${prefix}.genome.vcf.gz.tbi
    echo "" | gzip > ${prefix}.variants.vcf.gz
    touch ${prefix}.variants.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        strelka: \$( configureStrelkaSomaticWorkflow.py --version )
    END_VERSIONS
    """
}
