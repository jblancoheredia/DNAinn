process ANNOTSV_ANNOTSV {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b2/b202e030802ec909556961b542f15e0b37583755cebf08e899b3042a44f93ddb/data' :
        'community.wave.seqera.io/library/annotsv:3.4.2--6e6cee83703bd24c' }"

    input:
    tuple val(meta) ,
          val(meta1), path(sv_vcf), path(sv_vcf_tbi),
          val(meta2), path(snv_indel_vcf)
    path(gene_transcripts)
    path(candidate_genes)
    val(genome_version)
    path(annotations)

    output:
    tuple val(meta), path("*.annotated.tsv"), emit: tsv
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    AnnotSV \\
        -SVinputFile ${sv_vcf}       \\
        -txFile ${gene_transcripts}    \\
        -genomeBuild ${genome_version}   \\
        -annotationsDir ${annotations}     \\
        -outputFile ${prefix}.annotated.tsv  \\
        -candidateGenesFile ${candidate_genes} \\
        -candidateSnvIndelFiles ${snv_indel_vcf} \\
        -SVminSize 50  \\
        -outputDir ./    \\
        -SVinputInfo 1     \\
        -annotationMode full \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        annotsv: \$(echo \$(AnnotSV -help 2>&1 | head -n1 | sed 's/^AnnotSV //'))
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.annotated.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        annotsv: \$(echo \$(AnnotSV -help 2>&1 | head -n1 | sed 's/^AnnotSV //'))
    END_VERSIONS
    """
}
