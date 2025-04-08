process REPEATSEQ {
    tag "$meta.id"
    label 'process_high'

    memory { 6.GB * task.attempt }
    maxRetries 3
    errorStrategy 'retry'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/repeatseq:0.8.2' :
        'blancojmskcc/repeatseq:0.8.2' }"

    input:
    tuple val(meta),  path(bams) 
    tuple val(meta0), path(bais)
    tuple val(meta1), path(fasta)
    tuple val(meta2), path(fai)
    path(rep_regions)

    output:
    tuple val(meta),  path("*.vcf"),        emit: vcf
    tuple val(meta),  path("*.calls"),      emit: calls
    tuple val(meta2), path("*.repeatseq"),  emit: repeatseq
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat <<EOF > ${prefix}_repeatseq.vcf.tmp
    ##fileformat=VCFv4.1
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype likelihood">
    ##INFO=<ID=AL,Number=A,Type=Integer,Description="Allele Length Offset(s)">
    ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
    ##INFO=<ID=RU,Number=1,Type=String,Description="Repeat Unit">
    ##INFO=<ID=RL,Number=1,Type=Integer,Description="Reference Length of Repeat">
    EOF

    repeatseq \\
        -repeatseq \\
        -calls \\
        ${prefix}_1.bam \\
        ${fasta} \\
        ${rep_regions}

    cat *.calls >> ${prefix}_repeatseq.calls.tmp
    cat *.vcf | sed 8d >> ${prefix}_repeatseq.vcf.tmp
    cat *.repeatseq >> ${prefix}_repeatseq.repeatseq.tmp
    rm *.vcf
    rm *.calls
    rm *.repeatseq

    repeatseq \\
        -repeatseq \\
        -calls \\
        ${prefix}_10.bam \\
        ${fasta} \\
        ${rep_regions}

    cat *.calls >> ${prefix}_repeatseq.calls.tmp
    cat *.vcf | sed 8d >> ${prefix}_repeatseq.vcf.tmp
    cat *.repeatseq >> ${prefix}_repeatseq.repeatseq.tmp
    rm *.vcf
    rm *.calls
    rm *.repeatseq

    repeatseq \\
        -repeatseq \\
        -calls \\
        ${prefix}_11.bam \\
        ${fasta} \\
        ${rep_regions}

    cat *.calls >> ${prefix}_repeatseq.calls.tmp
    cat *.vcf | sed 8d >> ${prefix}_repeatseq.vcf.tmp
    cat *.repeatseq >> ${prefix}_repeatseq.repeatseq.tmp
    rm *.vcf
    rm *.calls
    rm *.repeatseq

    repeatseq \\
        -repeatseq \\
        -calls \\
        ${prefix}_12.bam \\
        ${fasta} \\
        ${rep_regions}

    cat *.calls >> ${prefix}_repeatseq.calls.tmp
    cat *.vcf | sed 8d >> ${prefix}_repeatseq.vcf.tmp
    cat *.repeatseq >> ${prefix}_repeatseq.repeatseq.tmp
    rm *.vcf
    rm *.calls
    rm *.repeatseq

    repeatseq \\
        -repeatseq \\
        -calls \\
        ${prefix}_13.bam \\
        ${fasta} \\
        ${rep_regions}

    cat *.calls >> ${prefix}_repeatseq.calls.tmp
    cat *.vcf | sed 8d >> ${prefix}_repeatseq.vcf.tmp
    cat *.repeatseq >> ${prefix}_repeatseq.repeatseq.tmp
    rm *.vcf
    rm *.calls
    rm *.repeatseq

    repeatseq \\
        -repeatseq \\
        -calls \\
        ${prefix}_14.bam \\
        ${fasta} \\
        ${rep_regions}

    cat *.calls >> ${prefix}_repeatseq.calls.tmp
    cat *.vcf | sed 8d >> ${prefix}_repeatseq.vcf.tmp
    cat *.repeatseq >> ${prefix}_repeatseq.repeatseq.tmp
    rm *.vcf
    rm *.calls
    rm *.repeatseq

    repeatseq \\
        -repeatseq \\
        -calls \\
        ${prefix}_15.bam \\
        ${fasta} \\
        ${rep_regions}

    cat *.calls >> ${prefix}_repeatseq.calls.tmp
    cat *.vcf | sed 8d >> ${prefix}_repeatseq.vcf.tmp
    cat *.repeatseq >> ${prefix}_repeatseq.repeatseq.tmp
    rm *.vcf
    rm *.calls
    rm *.repeatseq

    repeatseq \\
        -repeatseq \\
        -calls \\
        ${prefix}_16.bam \\
        ${fasta} \\
        ${rep_regions}

    cat *.calls >> ${prefix}_repeatseq.calls.tmp
    cat *.vcf | sed 8d >> ${prefix}_repeatseq.vcf.tmp
    cat *.repeatseq >> ${prefix}_repeatseq.repeatseq.tmp
    rm *.vcf
    rm *.calls
    rm *.repeatseq

    repeatseq \\
        -repeatseq \\
        -calls \\
        ${prefix}_17.bam \\
        ${fasta} \\
        ${rep_regions}

    cat *.calls >> ${prefix}_repeatseq.calls.tmp
    cat *.vcf | sed 8d >> ${prefix}_repeatseq.vcf.tmp
    cat *.repeatseq >> ${prefix}_repeatseq.repeatseq.tmp
    rm *.vcf
    rm *.calls
    rm *.repeatseq

    repeatseq \\
        -repeatseq \\
        -calls \\
        ${prefix}_18.bam \\
        ${fasta} \\
        ${rep_regions}

    cat *.calls >> ${prefix}_repeatseq.calls.tmp
    cat *.vcf | sed 8d >> ${prefix}_repeatseq.vcf.tmp
    cat *.repeatseq >> ${prefix}_repeatseq.repeatseq.tmp
    rm *.vcf
    rm *.calls
    rm *.repeatseq

    repeatseq \\
        -repeatseq \\
        -calls \\
        ${prefix}_19.bam \\
        ${fasta} \\
        ${rep_regions}

    cat *.calls >> ${prefix}_repeatseq.calls.tmp
    cat *.vcf | sed 8d >> ${prefix}_repeatseq.vcf.tmp
    cat *.repeatseq >> ${prefix}_repeatseq.repeatseq.tmp
    rm *.vcf
    rm *.calls
    rm *.repeatseq

    repeatseq \\
        -repeatseq \\
        -calls \\
        ${prefix}_1.bam \\
        ${fasta} \\
        ${rep_regions}

    cat *.calls >> ${prefix}_repeatseq.calls.tmp
    cat *.vcf | sed 8d >> ${prefix}_repeatseq.vcf.tmp
    cat *.repeatseq >> ${prefix}_repeatseq.repeatseq.tmp
    rm *.vcf
    rm *.calls
    rm *.repeatseq

    repeatseq \\
        -repeatseq \\
        -calls \\
        ${prefix}_2.bam \\
        ${fasta} \\
        ${rep_regions}

    cat *.calls >> ${prefix}_repeatseq.calls.tmp
    cat *.vcf | sed 8d >> ${prefix}_repeatseq.vcf.tmp
    cat *.repeatseq >> ${prefix}_repeatseq.repeatseq.tmp
    rm *.vcf
    rm *.calls
    rm *.repeatseq

    repeatseq \\
        -repeatseq \\
        -calls \\
        ${prefix}_20.bam \\
        ${fasta} \\
        ${rep_regions}

    cat *.calls >> ${prefix}_repeatseq.calls.tmp
    cat *.vcf | sed 8d >> ${prefix}_repeatseq.vcf.tmp
    cat *.repeatseq >> ${prefix}_repeatseq.repeatseq.tmp
    rm *.vcf
    rm *.calls
    rm *.repeatseq

    repeatseq \\
        -repeatseq \\
        -calls \\
        ${prefix}_21.bam \\
        ${fasta} \\
        ${rep_regions}

    cat *.calls >> ${prefix}_repeatseq.calls.tmp
    cat *.vcf | sed 8d >> ${prefix}_repeatseq.vcf.tmp
    cat *.repeatseq >> ${prefix}_repeatseq.repeatseq.tmp
    rm *.vcf
    rm *.calls
    rm *.repeatseq

    repeatseq \\
        -repeatseq \\
        -calls \\
        ${prefix}_22.bam \\
        ${fasta} \\
        ${rep_regions}

    cat *.calls >> ${prefix}_repeatseq.calls.tmp
    cat *.vcf | sed 8d >> ${prefix}_repeatseq.vcf.tmp
    cat *.repeatseq >> ${prefix}_repeatseq.repeatseq.tmp
    rm *.vcf
    rm *.calls
    rm *.repeatseq

    repeatseq \\
        -repeatseq \\
        -calls \\
        ${prefix}_3.bam \\
        ${fasta} \\
        ${rep_regions}

    cat *.calls >> ${prefix}_repeatseq.calls.tmp
    cat *.vcf | sed 8d >> ${prefix}_repeatseq.vcf.tmp
    cat *.repeatseq >> ${prefix}_repeatseq.repeatseq.tmp
    rm *.vcf
    rm *.calls
    rm *.repeatseq

    repeatseq \\
        -repeatseq \\
        -calls \\
        ${prefix}_4.bam \\
        ${fasta} \\
        ${rep_regions}

    cat *.calls >> ${prefix}_repeatseq.calls.tmp
    cat *.vcf | sed 8d >> ${prefix}_repeatseq.vcf.tmp
    cat *.repeatseq >> ${prefix}_repeatseq.repeatseq.tmp
    rm *.vcf
    rm *.calls
    rm *.repeatseq

    repeatseq \\
        -repeatseq \\
        -calls \\
        ${prefix}_5.bam \\
        ${fasta} \\
        ${rep_regions}

    cat *.calls >> ${prefix}_repeatseq.calls.tmp
    cat *.vcf | sed 8d >> ${prefix}_repeatseq.vcf.tmp
    cat *.repeatseq >> ${prefix}_repeatseq.repeatseq.tmp
    rm *.vcf
    rm *.calls
    rm *.repeatseq

    repeatseq \\
        -repeatseq \\
        -calls \\
        ${prefix}_6.bam \\
        ${fasta} \\
        ${rep_regions}

    cat *.calls >> ${prefix}_repeatseq.calls.tmp
    cat *.vcf | sed 8d >> ${prefix}_repeatseq.vcf.tmp
    cat *.repeatseq >> ${prefix}_repeatseq.repeatseq.tmp
    rm *.vcf
    rm *.calls
    rm *.repeatseq

    repeatseq \\
        -repeatseq \\
        -calls \\
        ${prefix}_7.bam \\
        ${fasta} \\
        ${rep_regions}

    cat *.calls >> ${prefix}_repeatseq.calls.tmp
    cat *.vcf | sed 8d >> ${prefix}_repeatseq.vcf.tmp
    cat *.repeatseq >> ${prefix}_repeatseq.repeatseq.tmp
    rm *.vcf
    rm *.calls
    rm *.repeatseq

    repeatseq \\
        -repeatseq \\
        -calls \\
        ${prefix}_8.bam \\
        ${fasta} \\
        ${rep_regions}

    cat *.calls >> ${prefix}_repeatseq.calls.tmp
    cat *.vcf | sed 8d >> ${prefix}_repeatseq.vcf.tmp
    cat *.repeatseq >> ${prefix}_repeatseq.repeatseq.tmp
    rm *.vcf
    rm *.calls
    rm *.repeatseq

    repeatseq \\
        -repeatseq \\
        -calls \\
        ${prefix}_9.bam \\
        ${fasta} \\
        ${rep_regions}

    cat *.calls >> ${prefix}_repeatseq.calls.tmp
    cat *.vcf | sed 8d >> ${prefix}_repeatseq.vcf.tmp
    cat *.repeatseq >> ${prefix}_repeatseq.repeatseq.tmp
    rm *.vcf
    rm *.calls
    rm *.repeatseq

    repeatseq \\
        -repeatseq \\
        -calls \\
        ${prefix}_MT.bam \\
        ${fasta} \\
        ${rep_regions}

    cat *.calls >> ${prefix}_repeatseq.calls.tmp
    cat *.vcf | sed 8d >> ${prefix}_repeatseq.vcf.tmp
    cat *.repeatseq >> ${prefix}_repeatseq.repeatseq.tmp
    rm *.vcf
    rm *.calls
    rm *.repeatseq

    repeatseq \\
        -repeatseq \\
        -calls \\
        ${prefix}_XY.bam \\
        ${fasta} \\
        ${rep_regions}

    cat *.calls >> ${prefix}_repeatseq.calls.tmp
    cat *.vcf | sed 8d >> ${prefix}_repeatseq.vcf.tmp
    cat *.repeatseq >> ${prefix}_repeatseq.repeatseq.tmp
    rm *.vcf
    rm *.calls
    rm *.repeatseq

    mv ${prefix}_repeatseq.vcf.tmp ${prefix}_repeatseq.vcf
    mv ${prefix}_repeatseq.calls.tmp ${prefix}_repeatseq.calls
    mv ${prefix}_repeatseq.repeatseq.tmp ${prefix}_repeatseq.repeatseq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatseq: \$(repeatseq -v 2>&1 | grep "RepeatSeq v" | sed 's/.*RepeatSeq v//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_repeatseq.vcf
    touch ${prefix}_repeatseq.calls
    touch ${prefix}_repeatseq.repeatseq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatseq: \$(repeatseq -v 2>&1 | grep "RepeatSeq v" | sed 's/.*RepeatSeq v//')
    END_VERSIONS
    """
}    
