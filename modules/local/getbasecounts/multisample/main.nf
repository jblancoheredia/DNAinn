process GETBASECOUNTS_MULTISAMPLE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/getbasecounts:1.2.2':
        'blancojmskcc/getbasecounts:1.2.2' }"

    input:
    tuple val(meta),  path(vcfs), path(tbis)
    tuple val(meta1), path(bam),  path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*.tvs"), emit: tvs
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    GetBaseCountsMultiSample \\
        --fasta ${fasta} \\
        --bam ${bam} \\
        --vcf ${prefix}.freebayes.vcf  \\
        --thread ${task.cpus} \\
        --output ${prefix}_freebayes_getbasecounts.tsv

    GetBaseCountsMultiSample \\
        --fasta ${fasta} \\
        --bam ${bam} \\
        --vcf ${prefix}.lofreq.somatic_final.snvs.vcf  \\
        --thread ${task.cpus} \\
        --output ${prefix}_lofreq_snvs_getbasecounts.tsv

    GetBaseCountsMultiSample \\
        --fasta ${fasta} \\
        --bam ${bam} \\
        --vcf ${prefix}.lofreq.somatic_final.indels.vcf  \\
        --thread ${task.cpus} \\
        --output ${prefix}_lofreq_indels_getbasecounts.tsv

    GetBaseCountsMultiSample \\
        --fasta ${fasta} \\
        --bam ${bam} \\
        --vcf ${prefix}.mutect2.filtered.vcf  \\
        --thread ${task.cpus} \\
        --output ${prefix}_mutect2_getbasecounts.tsv

    GetBaseCountsMultiSample \\
        --fasta ${fasta} \\
        --bam ${bam} \\
        --vcf ${prefix}.vardict.vcf  \\
        --thread ${task.cpus} \\
        --output ${prefix}_vardict_getbasecounts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getbasecounts: \$(GetBaseCountsMultiSample |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_getbasecounts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getbasecounts: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
