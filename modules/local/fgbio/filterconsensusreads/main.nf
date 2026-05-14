process FGBIO_FILTERCONSENSUSREADS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::fgbio=2.0.2 bioconda::samtools=1.16.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-69f5207f538e4de9ef3bae6f9a95c5af56a88ab8:82d3ec41f9f1227f7183d344be46f73365efa704-0' :
        'quay.io/biocontainers/mulled-v2-69f5207f538e4de9ef3bae6f9a95c5af56a88ab8:82d3ec41f9f1227f7183d344be46f73365efa704-0' }"

    input:
    tuple val(meta), path(consensus_bam)
    path(fasta)
    path(fai)
    val(min_reads)
    val(min_baseq)
    val(max_base_error_rate)
    val(max_read_error_rate)
    val(max_no_call_fraction)

    output:
    tuple val(meta), path("*.con.filtered.bam"),         path("*.con.filtered.bam.bai")           , emit: suplex_bam_bai
    tuple val(meta), path("*.con.duplex.filtered.bam"),  path("*.con.duplex.filtered.bam.bai")    , emit: duplex_bam_bai
    tuple val(meta), path("*.con.simplex.filtered.bam"), path("*.con.simplex.filtered.bam.bai")   , emit: simplex_bam_bai
    path "versions.yml"                                                                           , emit: versions

    script:
    def fgbio_args = task.ext.fgbio_args ?: ''
    def samtools_args = task.ext.samtools_args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem_gb = 8
    if (!task.memory) {
        log.info '[fgbio FilterConsensusReads] Available memory not known - defaulting to 8GB. Specify process memory requirements to change this.'
    } else {
        mem_gb = task.memory.giga
    }

    """
    fgbio \\
        -Xmx${mem_gb}g \\
        --tmp-dir=. \\
        --compression=0 \\
        FilterConsensusReads \\
        --input ${consensus_bam} \\
        --output /dev/stdout \\
        --ref ${fasta} \\
        --min-reads ${min_reads} \\
        --min-base-quality ${min_baseq} \\
        --max-base-error-rate ${max_base_error_rate} \\
        --max-read-error-rate ${max_read_error_rate} \\
        --max-no-call-fraction ${max_no_call_fraction} \\
        --reverse-per-base-tags true \\
        | samtools sort \\
        -@ ${task.cpus} \\
        -o ${prefix}.con.duplex.filtered.bam

    samtools index ${prefix}.con.duplex.filtered.bam

    samtools view -h ${consensus_bam} | egrep \"(^@|bD:i:0)\" | samtools sort -n > ${prefix}.tmp.bam

    samtools fixmate -u ${prefix}.tmp.bam - | samtools view -e 'flag.paired' -Sb - > ${prefix}.con.simplex.unfiltered.bam

    rm ${prefix}.tmp.bam

    fgbio \\
        -Xmx${mem_gb}g \\
        --tmp-dir=. \\
        --compression=0 \\
        FilterConsensusReads \\
        --input ${prefix}.con.simplex.unfiltered.bam \\
        --output /dev/stdout \\
        --ref ${fasta} \\
        --min-reads 3 3 0 \\
        --min-base-quality ${min_baseq} \\
        --max-base-error-rate ${max_base_error_rate} \\
        --max-read-error-rate ${max_read_error_rate} \\
        --max-no-call-fraction ${max_no_call_fraction} \\
        --reverse-per-base-tags true \\
        | samtools sort \\
        -@ ${task.cpus} \\
        -o ${prefix}.con.simplex.filtered.bam    

    samtools index ${prefix}.con.simplex.filtered.bam

    samtools merge -o ${prefix}.con.filtered.bam ${prefix}.con.duplex.filtered.bam ${prefix}.con.simplex.filtered.bam

    samtools index ${prefix}.con.filtered.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.con.filtered.bam
    touch ${prefix}.con.filtered.bam.bai
    touch ${prefix}.con.duplex.filtered.bam
    touch ${prefix}.con.duplex.filtered.bam.bai
    touch ${prefix}.con.simplex.filtered.bam
    touch ${prefix}.con.simplex.filtered.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """
}
