process VARDICTJAVA {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-731b8c4cf44d76e9aa181af565b9eee448d82a8c:edd70e76f3529411a748168f6eb1a61f29702123-0' :
        'biocontainers/mulled-v2-731b8c4cf44d76e9aa181af565b9eee448d82a8c:edd70e76f3529411a748168f6eb1a61f29702123-0' }"

    input:
    tuple val(meta), path(bams), path(bais)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)
    path(bed)

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.tbi"), emit: vcf
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-c 1 -S 2 -E 3'
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Don't run test scripts when -fisher has been used by vardictjava
    def run_test = !args.contains("-fisher")

    def somatic = bams instanceof List && bams.size() == 2 ? true : false
    def input = somatic ? "-b \"${bams[0]}|${bams[1]}\"" : "-b ${bams}"
    def test = run_test ? somatic ? "| testsomatic.R" : "| teststrandbias.R" : ""
    def convert_to_vcf = somatic ? "var2vcf_paired.pl" : "var2vcf_valid.pl"
    """
    export JAVA_OPTS='"-Xms${task.memory.toMega()/4}m" "-Xmx${task.memory.toGiga()}g" "-Dsamjdk.reference_fasta=${fasta}"'
    vardict-java \\
        ${args} \\
        ${input} \\
        -th ${task.cpus} \\
        -G ${fasta} \\
        ${bed} \\
    ${test} \\
    | ${convert_to_vcf} \\
        ${args2} \\
    | bgzip ${args3} --threads ${task.cpus} > ${prefix}.vcf.gz

    tabix -p vcf ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vardict-java: \$( realpath \$( command -v vardict-java ) | sed 's/.*java-//;s/-.*//' )
        var2vcf_valid.pl: \$( var2vcf_valid.pl -h | sed '2!d;s/.* //' )
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
        var2vcf_valid.pl: \$( var2vcf_valid.pl -h | sed '2!d;s/.* //' )
    END_VERSIONS
    """
}
