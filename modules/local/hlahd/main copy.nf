process HLAHD {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/hlahd:1.7.1' :
        'blancojmskcc/hlahd:1.7.1' }"

    input:

    tuple val(meta), path(bam), path(bai)
    val(genome)
    val(read_length)

    output:
    tuple val(meta), path("*.txt"), emit: hla_type
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def region = genome == 'GRCh38' ? '6:28510120-33480577' :
                 genome == 'GRCh37' ? '6:28477797-33448354' :
                 { error("Unsupported genome: ${genome}") }()
    def length = read_length == 100 ? 95 :
                 read_length == 150 ? 145 :
                 { error("Unsupported read length: ${read_length}") }()
    def freq_data  = "/opt/hla-hd/hlahd.1.7.1/freq_data/"
    def gene_split = "/opt/hla-hd/hlahd.1.7.1/HLA_gene.split.3.50.0.txt"
    def dictionary = "/opt/hla-hd/hlahd.1.7.1/dictionary/"
    """
    samtools \\
        view \\
        -@ ${task.cpus} \\
        -h \\
        -b ${bam} \\
        ${region} \\
        > ${prefix}.mhc.bam

    samtools \\
        view \\
        -@ ${task.cpus} \\
        -b \\
        -f 4 \\
        ${bam} \\
        > ${prefix}.unmap.bam

    samtools \\
        merge \\
        -@ ${task.cpus} \\
        ${prefix}.merged.bam \\
        ${prefix}.unmap.bam \\
        ${prefix}.mhc.bam

    picard \\
        SamToFastq \\
        -I ${prefix}.merged.bam \\
        -F ${prefix}.hlatmp.1.fastq \\
        -F2 ${prefix}.hlatmp.2.fastq \\
        -FU ${prefix}.hlatmp.0.fastq \\
        --VALIDATION_STRINGENCY SILENT


    cat ${prefix}.hlatmp.1.fastq | awk '{if(NR%4 == 1){O=\$0;gsub(\"/1\",\" 1\",O);print O}else{print \$0}}' > ${prefix}.hla.1.fastq
    cat ${prefix}.hlatmp.2.fastq | awk '{if(NR%4 == 1){O=\$0;gsub(\"/2\",\" 2\",O);print O}else{print \$0}}' > ${prefix}.hla.2.fastq

    hlahd.sh \\
        -t ${task.cpus} \\
        -m ${length} \\
        -c 0.95 \\
        -f ${freq_data} \\
        ${prefix}.hla.1.fastq \\
        ${prefix}.hla.2.fastq \\
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
