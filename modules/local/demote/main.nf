process DEMOTE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/demote:1.0':
        'blancojmskcc/demote:1.0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.demoted.bam"), path("*.demoted.bam.bai"), emit: bam
    path "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def prefix   = "${meta.id}"
    def fixbam   = "${prefix}.demoted.fixmate.bam"
    def namebam  = "${prefix}.demoted.name.bam"
    def unsorted = "${prefix}.demoted.unsorted.bam"
    def coordbam = "${prefix}.demoted.fixmate.sorted.bam"
    def cleanbam = "${prefix}.demoted.bam"
    """
    Demote \\
        --in  ${bam} \\
        --out ${unsorted} \\
        --log ${prefix}.demote.log \\
        ${args}

    samtools sort -@ ${task.cpus} -n -o ${namebam} ${unsorted}

    samtools fixmate -@ ${task.cpus} -m ${namebam} ${fixbam}

    samtools sort -@ ${task.cpus} -o ${coordbam} ${fixbam}

    samtools markdup -@ ${task.cpus} -r ${coordbam} ${cleanbam}

    samtools view -@ ${task.cpus} -h ${cleanbam} | awk 'BEGIN{OFS="\\t"} /^@/{print; next} {for(i=1;i<=NF;i++) if(i==11) {gsub(/[K-~]/,\"F\",\$i)} print}' | samtools view -@ ${task.cpus} -b -o ${cleanbam}.tmp

    samtools index -@ ${task.cpus} ${cleanbam}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        demoted: 1.0
    END_VERSIONS
    """
    stub:
    def prefix = "${meta.id}"
    """
    touch ${prefix}.demoted.bam
    touch ${prefix}.demoted.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        demoted: 1.0
    END_VERSIONS
    """
}

