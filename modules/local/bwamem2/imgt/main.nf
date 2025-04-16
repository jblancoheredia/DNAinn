process ALIGN_HLA_IMGT {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::fgbio=2.4.0 bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/umi_aligner:1.0.0' :
        'blancojmskcc/umi_aligner:1.0.0' }"

    input:
    tuple val(meta),  path(reads)
    path(bwa2)
    path(fasta)
    path(fasta_fai)

    output:
    tuple val(meta), path("*.hla.imgt.mapped.bam"), emit: bam
    tuple val(meta), path("*.hla.imgt.mapped.bai"), emit: bai
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.bwa_args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # The real path to the BWAMEM2 index prefix
    BWA_INDEX_PREFIX=\$(find -L ./ -name "*.amb" | sed 's/.amb//')

    bwa-mem2 mem \\
        ${args} \\
        -t ${task.cpus} \\
        \$BWA_INDEX_PREFIX \\
        ${reads[0]} \\
        ${reads[1]} | \\
        samtools view -bS - | \\
        samtools sort  -@ ${task.cpus} -o ${prefix}_sorted.bam - && \\
        samtools index -@ ${task.cpus} ${prefix}_sorted.bam && \\
        samtools calmd -@ ${task.cpus} -b ${prefix}_sorted.bam ${fasta} > ${prefix}.hla.imgt.mapped.bam && \\
        samtools index -@ ${task.cpus} ${prefix}.hla.imgt.mapped.bam ${prefix}.hla.imgt.mapped.bai
        
    rm ${prefix}_sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>&1 | tail -n 1)
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.hla.imgt.mapped.bam
    touch ${prefix}.hla.imgt.mapped.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>&1 | tail -n 1)
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

}

process ALIGN_BAM_FIN {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::fgbio=2.4.0 bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/umi_aligner:1.0.0' :
        'blancojmskcc/umi_aligner:1.0.0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta0), path(duplex_bam), path(duplex_bai)
    tuple val(meta1), path(simplex_bam), path(simplex_bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)
    tuple val(meta4), path(dict)
    tuple val(meta5), path(bwa_dir)

    output:
    tuple val(meta), path("*.mapped.bam"),          emit: bam
    tuple val(meta), path("*.mapped.duplex.bam"),   emit: duplex_bam
    tuple val(meta), path("*.mapped.simplex.bam"),  emit: simplex_bam
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def samtools_fastq_args = task.ext.samtools_fastq_args ?: ''
    def samtools_sort_args = task.ext.samtools_sort_args ?: ''
    def bwa_args = task.ext.bwa_args ?: ''
    def fgbio_args = task.ext.fgbio_args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fgbio_mem_gb = 4
    def extra_command = ""

    if (!task.memory) {
        log.info '[fgbio ZipperBams] Available memory not known - defaulting to 4GB. Specify process memory requirements to change this.'
    } else if (fgbio_mem_gb > task.memory.giga) {
        if (task.memory.giga < 2) {
            fgbio_mem_gb = 1
        } else {
            fgbio_mem_gb = task.memory.giga - 1
        }
    }
    """
    # The real path to the BWA index prefix`
    BWA_INDEX_PREFIX=\$(find -L ./ -name "*.amb" | sed 's/.amb//')

    samtools fastq ${samtools_fastq_args} ${bam} \\
        | bwa-mem2 mem ${bwa_args} -t ${task.cpus} -p -K 150000000 -Y \$BWA_INDEX_PREFIX - \\
        | fgbio -Xmx${fgbio_mem_gb}g \\
            --compression 1 \\
            --async-io=true \\
            ZipperBams \\
            --unmapped ${bam} \\
            --ref ${fasta} \\
            --output ${prefix}.mapped.bam \\
            --tags-to-reverse Consensus \\
            --tags-to-revcomp Consensus \\
            --sort TemplateCoordinate

    samtools fastq ${samtools_fastq_args} ${duplex_bam} \\
        | bwa-mem2 mem ${bwa_args} -t ${task.cpus} -p -K 150000000 -Y \$BWA_INDEX_PREFIX - \\
        | fgbio -Xmx${fgbio_mem_gb}g \\
            --compression 1 \\
            --async-io=true \\
            ZipperBams \\
            --unmapped ${duplex_bam} \\
            --ref ${fasta} \\
            --output ${prefix}.mapped.duplex.bam \\
            --tags-to-reverse Consensus \\
            --tags-to-revcomp Consensus \\
            --sort TemplateCoordinate

    samtools fastq ${samtools_fastq_args} ${simplex_bam} \\
        | bwa-mem2 mem ${bwa_args} -t ${task.cpus} -p -K 150000000 -Y \$BWA_INDEX_PREFIX - \\
        | fgbio -Xmx${fgbio_mem_gb}g \\
            --compression 1 \\
            --async-io=true \\
            ZipperBams \\
            --unmapped ${simplex_bam} \\
            --ref ${fasta} \\
            --output ${prefix}.mapped.simplex.bam \\
            --tags-to-reverse Consensus \\
            --tags-to-revcomp Consensus \\
            --sort TemplateCoordinate

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>&1 | tail -n 1)
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mapped.bam
    touch ${prefix}.mapped.duplex.bam
    touch ${prefix}.mapped.simplex.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>&1 | tail -n 1)
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

}