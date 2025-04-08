process TELOMEREHUNTER_FULL {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/telomerehunter:1.1.0' :
        'blancojmskcc/telomerehunter:1.1.0' }"

    input:
    tuple val(meta) , path(tumour_bam)
    tuple val(meta1), path(tumour_bai)
    path(normal_bam)
    path(normal_bai)
    path(banding)

    output:
    tuple val(meta), path("*_fusions.bam"), path("*_fusions.bai"), emit: bam
    tuple val(meta), path("*.pdf")                               , emit: pdf
    tuple val(meta), path("*.tsv")                               , emit: tsv
    path "versions.yml"                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def outputDir = "${prefix}_telomerehunter"
    """
    mkdir -p ${outputDir}

    telomerehunter --outPath ${outputDir}/ \\
    --pid ${prefix} \\
    --inputBamTumor ${tumour_bam} \\
    --inputBamControl ${normal_bam} \\
    --bandingFile ${banding} \\
    --repeatThreshold 4 \\
    --mappingQualityThreshold 6 \\
    --repeats TTAGGG TGAGGG TCAGGG TTCGGG TTGGGG TTTGGG ATAGGG CATGGG CTAGGG GTAGGG TAAGGG \\
    --repeatsContext TTAGGG TGAGGG TCAGGG TTCGGG TTGGGG TTTGGG ATAGGG CATGGG CTAGGG GTAGGG TAAGGG \\
    --parallel \\
    -p1 -p2 -p3 -p4 -p5 -p6 -p7

    cp -r ${outputDir}/${prefix}/* .

    samtools view -H tumor_TelomerCnt_${prefix}/${prefix}_filtered.bam > ${prefix}_fusions
    samtools view tumor_TelomerCnt_${prefix}/${prefix}_filtered.bam >> ${prefix}_fusions
    samtools view tumor_TelomerCnt_${prefix}/${prefix}_filtered_intrachromosomal.bam >> ${prefix}_fusions
    samtools view tumor_TelomerCnt_${prefix}/${prefix}_filtered_intratelomeric.bam >> ${prefix}_fusions
    samtools view tumor_TelomerCnt_${prefix}/${prefix}_filtered_junctionspanning.bam >> ${prefix}_fusions
    samtools view tumor_TelomerCnt_${prefix}/${prefix}_filtered_subtelomeric.bam >> ${prefix}_fusions

    cp ${projectDir}/bin/sed.sh .
    bash sed.sh ${prefix}

    samtools \\
    sort \\
    -@ ${task.cpus} \\
    ${prefix}_fusions \\
    > ${prefix}_fusions.bam

    samtools \\
    index \\
    -@ ${task.cpus} \\
    ${prefix}_fusions.bam \\
    ${prefix}_fusions.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        telomerehunter: "1.1.0"
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    touch ${prefix}.pdf
    touch ${prefix}_fusions.bam
    touch ${prefix}_fusions.bai
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        telomerehunter: "1.1.0"
    END_VERSIONS
    """
}
