process SCALPEL_EXPORT {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://drgiovianco/scalpel:0.5.3' :
        'drgiovianco/scalpel:0.5.3' }"

    input:
    tuple val(meta),  path(bam), path(bai)
    tuple val(meta1), path(fasta)
    tuple val(meta1), path(fai)
    tuple val(meta1), path(bed)

    output:
    tuple val(meta),  path("*.db")  , emit: db
    tuple val(meta2), path("*.vcf") , emit: vcf
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.5.3' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    scalpel-discovery \\
        --single \\
        --numprocs ${task.cpus} \\
        --bam $bam \\
        --bed $bed \\
        --dir . \\
        --ref $fasta

    scalpel-export \\
        --single \\
        --db ${prefix}_database.db \\
        --bed $bed \\
        --dir . \\
        --ref genome.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scalpel: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.5.3' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}_database.db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scalpel: $VERSION
    END_VERSIONS
    """
}    
