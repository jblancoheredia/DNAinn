process GATK4_APPLYBQSR {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.5.0.0--py36hdfd78af_0':
        'quay.io/biocontainers/gatk4:4.5.0.0--py36hdfd78af_0' }"

    input:
    tuple val(meta),  path(input), path(input_index)
    tuple val(meta2), path(bqsr_table)
    path(fasta)
    path(fai)
    path(dict)

    output:
    tuple val(meta), path("*.bam") , emit: bam,  optional: true
    tuple val(meta), path("*.cram"), emit: cram, optional: true
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK ApplyBQSR] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        ApplyBQSR \\
        --input $input \\
        --output ${prefix}.${input.getExtension()} \\
        --reference $fasta \\
        --bqsr-recal-file $bqsr_table \\
        --tmp-dir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_extension = input.getExtension()
    def output_extension = input_extension == 'bam' ? 'bam' : 'cram'
    """
    touch ${prefix}.${output_extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
