process CNVKIT_REFERENCE {
    tag "$fasta"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cnvkit:0.9.11--pyhdfd78af_0':
        'quay.io/biocontainers/cnvkit:0.9.11--pyhdfd78af_0' }"

    input:
    tuple val(meta),  path(fasta)
    tuple val(meta1), path(targets)
    tuple val(meta2), path(antitargets)

    output:
    path("versions.yml"), emit: versions
    tuple val(meta), path("*.cnn"), emit: cnn

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: targets.BaseName

    """
    cnvkit.py \\
        reference \\
        --fasta $fasta \\
        --targets $targets \\
        --antitargets $antitargets \\
        --output ${prefix}.reference.cnn \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvkit: \$(cnvkit.py version | sed -e "s/cnvkit v//g")
    END_VERSIONS
    """
}
