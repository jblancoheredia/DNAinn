process SAMTOOLS_MPILEUP {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'quay.io/biocontainers/samtools:1.21--h50ea8bc_0' }"
    input:
    tuple val(meta),  path(bam), path(bai)
    tuple val(meta1), path(intervals)
    path  fasta

    output:
    tuple val(meta), path("*.pileup.gz")    , emit: pup
    tuple val(meta), path("*.pileup.gz.tbi"), emit: tbi
    path  "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def intervals = intervals ? "-l ${intervals}" : ""
    """
    samtools mpileup \\
        --fasta-ref $fasta \\
        --output ${prefix}.pileup \\
        $args \\
        $intervals \\
        $bam

    bgzip ${prefix}.pileup

    tabix ${prefix}.pileup.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pileup.gz
    touch ${prefix}.pileup.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
