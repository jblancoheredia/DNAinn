process MSISENSORPRO_RAW {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/msisensor-pro:1.2.0--hfc31af2_0' :
        'quay.io/biocontainers/msisensor-pro:1.2.0--hfc31af2_0' }"

    input:
    tuple val(meta),   path(bam)
    tuple val(meta1) , path(bai)
    tuple val(meta2),  path(hom_ms)

    output:
    tuple val(meta), path("${meta.id}.summary.raw"), emit: summary
    tuple val(meta), path("${meta.id}.msi_uns.raw"), emit: msi_uns
    tuple val(meta), path("${meta.id}.msi_dis.raw"), emit: msi_dis
    tuple val(meta), path("${meta.id}.msi_all.raw"), emit: msi_all
    path "versions.yml",                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    msisensor-pro \\
        pro \\
        -d ${hom_ms} \\
        -t ${bam} \\
        -o ${prefix} \\
        $args

    mv ${prefix}            ${prefix}.summary.raw
    mv ${prefix}_unstable   ${prefix}.msi_uns.raw
    mv ${prefix}_dis        ${prefix}.msi_dis.raw
    mv ${prefix}_all        ${prefix}.msi_all.raw

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        msisensor-pro: \$(msisensor-pro 2>&1 | sed -nE 's/Version:\\sv([0-9]\\.[0-9])/\\1/ p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.summary.raw
    touch ${prefix}.msi_uns.raw
    touch ${prefix}.msi_dis.raw
    touch ${prefix}.msi_all.raw

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        msisensor-pro: \$(msisensor-pro 2>&1 | sed -nE 's/Version:\\sv([0-9]\\.[0-9])/\\1/ p')
    END_VERSIONS
    """
}

process MSISENSORPRO_FIN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/msisensor-pro:1.2.0--hfc31af2_0' :
        'quay.io/biocontainers/msisensor-pro:1.2.0--hfc31af2_0' }"

    input:
    tuple val(meta),   path(bam), path(bai)
    tuple val(meta2),  path(hom_ms)

    output:
    tuple val(meta), path("${meta.id}.summary.fin"), emit: summary
    tuple val(meta), path("${meta.id}.msi_uns.fin"), emit: msi_uns
    tuple val(meta), path("${meta.id}.msi_dis.fin"), emit: msi_dis
    tuple val(meta), path("${meta.id}.msi_all.fin"), emit: msi_all
    path "versions.yml",                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    msisensor-pro \\
        pro \\
        -d ${hom_ms} \\
        -t ${bam} \\
        -o ${prefix} \\
        $args

    mv ${prefix}            ${prefix}.summary.fin
    mv ${prefix}_unstable   ${prefix}.msi_uns.fin
    mv ${prefix}_dis        ${prefix}.msi_dis.fin
    mv ${prefix}_all        ${prefix}.msi_all.fin

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        msisensor-pro: \$(msisensor-pro 2>&1 | sed -nE 's/Version:\\sv([0-9]\\.[0-9])/\\1/ p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.summary.fin
    touch ${prefix}.msi_uns.fin
    touch ${prefix}.msi_dis.fin
    touch ${prefix}.msi_all.fin

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        msisensor-pro: \$(msisensor-pro 2>&1 | sed -nE 's/Version:\\sv([0-9]\\.[0-9])/\\1/ p')
    END_VERSIONS
    """
}
