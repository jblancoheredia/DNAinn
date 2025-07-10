process VARDICTJAVA {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-731b8c4cf44d76e9aa181af565b9eee448d82a8c:edd70e76f3529411a748168f6eb1a61f29702123-0' :
        'biocontainers/mulled-v2-731b8c4cf44d76e9aa181af565b9eee448d82a8c:edd70e76f3529411a748168f6eb1a61f29702123-0' }"

    input:
    tuple val(meta),  path(tbam), path(tbai), path(nbam), path(nbai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)
    path(bed)

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.tbi"), emit: vcf
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export JAVA_OPTS='"-Xms${task.memory.toMega()/4}m" "-Xmx${task.memory.toGiga()}g" "-Dsamjdk.reference_fasta=${fasta}"'

    vardict-java \\
        ${args} \\
        -b \"${tbam}|${nbam}\" \\
        -th ${task.cpus} \\
        -G ${fasta} \\
        ${bed} \\
        | var2vcf_paired.pl \\
        ${args2} \\
        | bgzip --threads ${task.cpus} > ${prefix}.vcf.gz

    tabix -p vcf ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vardict-java: \$( realpath \$( command -v vardict-java ) | sed 's/.*java-//;s/-.*//' )
        var2vcf_paired.pl: \$( var2vcf_paired.pl -h | sed '2!d;s/.* //' )
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vardict-java: \$( realpath \$( command -v vardict-java ) | sed 's/.*java-//;s/-.*//' )
        var2vcf_paired.pl: \$( var2vcf_paired.pl -h | sed '2!d;s/.* //' )
    END_VERSIONS
    """
}