process FASTP {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::fastp=0.23.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastp:0.23.4--h5f740d0_0' :
        'quay.io/biocontainers/fastp:0.23.4--h5f740d0_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    path "versions.yml"                       , emit: versions
    tuple val(meta), path('*.html')           , emit: html
    tuple val(meta), path('*.json')           , emit: json
    tuple val(meta), path('*.fail.fastq.gz')  , optional:true, emit: reads_fail
    tuple val(meta), path('*.fastp.fastq.gz') , optional:true, emit: reads

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fastp \\
        --in1 ${reads[0]} \\
        --out1 ${prefix}_R1.fastp.fastq.gz \\
        ${meta.single_end ? "" : "--in2 ${reads[1]}"} \\
        ${meta.single_end ? "" : "--out2 ${prefix}_R2.fastp.fastq.gz"} \\
        --json ${prefix}.fastp.json \\
        --html ${prefix}.fastp.html \\
        --failed_out ${prefix}.paired.fail.fastq.gz \\
        --unpaired1 ${prefix}_R1.fail.fastq.gz \\
        --unpaired2 ${prefix}_R2.fail.fastq.gz \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_1.fastp.fastq.gz
    ${meta.single_end ? "" : "touch ${prefix}_2.fastp.fastq.gz"}
    touch ${prefix}.fastp.json
    touch ${prefix}.fastp.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
    END_VERSIONS
    """
}
