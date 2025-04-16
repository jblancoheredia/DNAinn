process HLAHD {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/hlahd:1.7.1' :
        'blancojmskcc/hlahd:1.7.1' }"

    input:

    tuple val(meta) , path(fastqs)
    val(read_length)

    output:
    tuple val(meta), path("*.txt"), emit: hla_type
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def length = read_length == 100 ? 95 :
                 read_length == 150 ? 145 :
                 { error("Unsupported read length: ${read_length}") }()
    def freq_data  = "/opt/hla-hd/hlahd.1.7.1/freq_data/"
    def gene_split = "/opt/hla-hd/hlahd.1.7.1/HLA_gene.split.3.50.0.txt"
    def dictionary = "/opt/hla-hd/hlahd.1.7.1/dictionary/"
    """
    hlahd.sh \\
        -t ${task.cpus} \\
        -m ${length} \\
        -c 0.95 \\
        -f ${freq_data} \\
        ${fastqs} \\
        ${gene_split} \\
        ${dictionary} \\
        ${prefix} \\
        .

    cp ${prefix}/result/${prefix}_final.result.txt ${prefix}_hla_hd.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hlahd: "1.7.1"
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hlahd: "1.7.1"
    END_VERSIONS
    """
}
