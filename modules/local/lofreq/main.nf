process LOFREQ {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/lofreq:2.1.5' :
        'blancojmskcc/lofreq:2.1.5' }"

    input:
    tuple val(meta) , path(tbam), path(tbai), path(nbam), path(nbai)
    tuple val(meta1), path(intervals)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    path(known_sites_tbi)
    path(known_sites)
    
    output:
    tuple val(meta), path("*.somatic_final.snvs.vcf.gz")  , path("*.somatic_final.snvs.vcf.gz.tbi")  , emit: vcf_snvs
    tuple val(meta), path("*.somatic_final.indels.vcf.gz"), path("*.somatic_final.indels.vcf.gz.tbi"), emit: vcf_indels
    path "versions.yml"                                                                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    lofreq indelqual \\
        --dindel \\
        -f ${fasta} \\
        -o Tumour_dindel.bam \\
        ${tbam}

    samtools index -@ ${task.cpus} Tumour_dindel.bam

    lofreq indelqual \\
        --dindel \\
        -f ${fasta} \\
        -o Normal_dindel.bam \\
        ${nbam}

    samtools index -@ ${task.cpus} Normal_dindel.bam

    lofreq \\
        somatic \\
        ${args} \\
        -f ${fasta} \\
        -l ${intervals} \\
        -d ${known_sites} \\
        -n Normal_dindel.bam \\
        -t Tumour_dindel.bam \\
        --threads ${task.cpus} \\
        -o ${prefix}.lofreq.

    rm Normal_dindel.bam* Tumour_dindel.bam*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lofreq: \$(echo \$(lofreq version 2>&1) | sed 's/^version: //; s/ *commit.*\$//')
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
	touch ${prefix}.lofreq.somatic_final.indels.vcf.gz
	touch ${prefix}.lofreq.somatic_final.indels.vcf.gz.tbi
	touch ${prefix}.lofreq.somatic_final.snvs.vcf.gz
	touch ${prefix}.lofreq.somatic_final.snvs.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lofreq: \$(echo \$(lofreq version 2>&1) | sed 's/^version: //; s/ *commit.*\$//')
    END_VERSIONS
    """
}
