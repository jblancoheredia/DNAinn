process SAMBLASTER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-19fa9f1a5c3966b63a24166365e81da35738c5ab:60ebac4ad9c6530c0d7bf6844f52ec6916e1e0b1-0' :
        'quay.io/biocontainers/mulled-v2-19fa9f1a5c3966b63a24166365e81da35738c5ab:60ebac4ad9c6530c0d7bf6844f52ec6916e1e0b1-0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path('*splits.fastq.gz'), emit: split_reads
    tuple val(meta), path("*splits.bam")     , emit: split_bam
    tuple val(meta), path("*tagged.bam")     , emit: bam
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools view -h ${bam} | \\
    samblaster ${args} -s ${prefix}_splits.sam | \\
    samtools view -Sb - > ${prefix}_tagged.bam

    samtools view ${prefix}_splits.sam | \\
        awk '{
            name=\$1;
            sub(/_[12]\$/, "", name);
            sequence=\$10;
            quality=\$11;

            if (NR % 2 == 1) {
                print "@" name >> "${prefix}_R1.splits.fastq";
                print sequence >> "${prefix}_R1.splits.fastq";
                print "+" >> "${prefix}_R1.splits.fastq";
                print quality >> "${prefix}_R1.splits.fastq"
            }
            else {
                print "@" name >> "${prefix}_R2.splits.fastq";
                print sequence >> "${prefix}_R2.splits.fastq";
                print "+" >> "${prefix}_R2.splits.fastq";
                print quality >> "${prefix}_R2.splits.fastq"
            }
        }'

    gzip ${prefix}_R1.splits.fastq
    gzip ${prefix}_R2.splits.fastq

    samtools view -b ${prefix}_splits.sam -o ${prefix}_splits.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samblaster: \$( samblaster -h 2>&1 | head -n 1 | sed 's/^samblaster: Version //' )
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_R1.splits.fastq.gz
    touch ${prefix}_R2.splits.fastq.gz
    touch ${prefix}_splits.bam
    touch ${prefix}_tagged.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samblaster: \$( samblaster -h 2>&1 | head -n 1 | sed 's/^samblaster: Version //' )
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
