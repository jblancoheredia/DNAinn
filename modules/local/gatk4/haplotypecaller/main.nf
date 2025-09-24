process GATK4_HAPLOTYPECALLER {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ced519873646379e287bc28738bdf88e975edd39a92e7bc6a34bccd37153d9d0/data'
        : 'community.wave.seqera.io/library/gatk4_gcnvkernel:edb12e4f0bf02cd3'}"

    input:
    tuple val(meta),  path(t_bam), path(t_bai), path(n_bam), path(n_bai)
    tuple val(meta2), path(dbsnp), path(dbsnp_tbi)
    tuple val(meta3), path(intervals)
    tuple val(meta4), path(fasta)
    tuple val(meta5), path(dict)
    tuple val(meta6), path(fai)
    
    output:
    tuple val(meta), path("*.vcf.gz"),        emit: vcf
    tuple val(meta), path("*.tbi"),           emit: tbi, optional: true
    tuple val(meta), path("*.realigned.bam"), emit: bam, optional: true
    path "versions.yml",                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK HaplotypeCaller] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        HaplotypeCaller \\
        --native-pair-hmm-threads ${task.cpus} \\
        --bam-output ${prefix}.realigned.bam \\
        --emit-ref-confidence BP_RESOLUTION \\
        --output ${prefix}.vcf.gz \\
        --intervals ${intervals} \\
        --pcr-indel-model NONE \\
        --reference ${fasta} \\
        --dbsnp ${dbsnp} \\
        --input ${t_bam} \\
        --tmp-dir . \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bamout_command = args.contains("--bam-writer-type") ? "--bam-output ${prefix.replaceAll('.g\\s*$', '')}.realigned.bam" : ""

    def stub_realigned_bam = bamout_command ? "touch ${prefix.replaceAll('.g\\s*$', '')}.realigned.bam" : ""
    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    ${stub_realigned_bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
