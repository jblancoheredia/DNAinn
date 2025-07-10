process GATK4_MUTECT2 {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.5.0.0--py36hdfd78af_0':
        'quay.io/biocontainers/gatk4:4.5.0.0--py36hdfd78af_0' }"

    input:
    tuple val(meta),  path(tbam), path(tbai), path(nbam), path(nbai)
    tuple val(meta1), path(intervals)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(dict)
    tuple val(meta5), path(txt)
    path(germline_resource)
    path(germline_resource_tbi)
    path(panel_of_normals)
    path(panel_of_normals_tbi)

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.tbi"),       emit: vcf
    tuple val(meta), path("*.stats"),                       emit: stats
    tuple val(meta), path("*.f1r2.tar.gz"), optional:true,  emit: f1r2
    path "versions.yml",                                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK Mutect2] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    normal_sample_name=\$(cat ${txt})

    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        Mutect2 \\
        -I ${tbam} \\
        -I ${nbam} \\
        -R ${fasta} \\
        -normal \$normal_sample_name \\
        -O ${prefix}.vcf.gz \\
        --intervals ${intervals} \\
        --panel-of-normals ${panel_of_normals} \\
        --germline-resource ${germline_resource} \\
        --tmp-dir . \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    touch ${prefix}.vcf.gz.stats
    touch ${prefix}.f1r2.tar.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
