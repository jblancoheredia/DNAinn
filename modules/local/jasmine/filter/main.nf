process JASMINE_FILTER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/jasminesv:1.1.5--hdfd78af_0':
        'quay.io/biocontainers/jasminesv:1.1.5--hdfd78af_0' }"

    input:
    tuple val(meta) , path(delly_vcf)
    tuple val(meta1), path(svaba_vcf)
    tuple val(meta2), path(manta_vcf)
    tuple val(meta3), path(tiddit_vcf)
    tuple val(meta4), path(gridss_vcf)
    tuple val(meta5), path(recall_vcf)
    tuple val(meta6), path(fasta)
    tuple val(meta7), path(fasta_fai)
    tuple val(meta8), path(bam) , path(bai)

    output:
    tuple val(meta), path("*.filtered.bed"), emit: bed
    tuple val(meta), path("*.filtered.vcf"), emit: vcf
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def iris_argument = args != '' ? "--run_iris iris_args=${args}" : ""
    """
    ls *.vcf > vcfs.txt

    jasmine \\
        file_list=vcfs.txt \\
        out_file=${prefix}.filtered.vcf \\
        threads=${task.cpus} \\
        genome_file=${fasta} \\
        --dup_to_ins \\
        --normalize_type \\
        ${iris_argument}

    python ${baseDir}/bin/vcf2bed.py -i ${prefix}.filtered.vcf -o ${prefix}.filtered.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jasminesv: \$(echo \$(jasmine 2>&1 | grep "version" | sed 's/Jasmine version //'))
        bgzip: \$(echo \$(bgzip -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.filtered.vcf
    touch ${prefix}.filtered.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jasminesv: \$(echo \$(jasmine 2>&1 | grep "version" | sed 's/Jasmine version //'))
        bgzip: \$(echo \$(bgzip -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
