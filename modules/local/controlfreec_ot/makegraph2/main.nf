process CONTROLFREEC_OT_MAKEGRAPH2 {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/control-freec:11.6--h1b792b2_1' :
        'quay.io/biocontainers/control-freec:11.6--h1b792b2_1' }"

    input:
    tuple val(meta),  path(ratio)
    tuple val(meta1), path(baf)
    path(makeGraph2)

    output:
    tuple val(meta), path("*_BAF.png")       , emit: png_baf
    tuple val(meta), path("*_ratio.log2.png"), emit: png_ratio_log2
    tuple val(meta), path("*_ratio.png")     , emit: png_ratio

    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def baf = baf ?: ""
    def VERSION = '11.6' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    cat ${makeGraph2} | R --slave --args ${ratio} ${baf}

    mv ${prefix}_BAF.txt.png ${prefix}_BAF.png
    mv ${prefix}_ratio.txt.log2.png ${prefix}_ratio.log2.png
    mv ${prefix}_ratio.txt.png ${prefix}_ratio.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        controlfreec: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '11.6' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}_BAF.png
    touch ${prefix}_ratio.log2.png
    touch ${prefix}_ratio.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        controlfreec: $VERSION
    END_VERSIONS
    """
}
