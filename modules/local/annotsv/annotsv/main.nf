process ANNOTSV_ANNOTSV {
    tag "$meta.patient_id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/dnainn_annotsv:3.5':
        'blancojmskcc/dnainn_annotsv:3.5' }"

    input:
    tuple val(meta) ,
          val(meta1), path(sv_vcf), path(sv_vcf_tbi),
          val(meta2), path(snv_indel_vcf)
    path(gene_transcripts)
    path(candidate_genes)
    val(genome_version)
    path(annotations)

    output:
    tuple val(meta), path("*.log")          , emit: log
    tuple val(meta), path("*.annotated.tsv"), emit: tsv
    tuple val(meta), path("*.annotated.vcf"), emit: vcf
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.patient_id}"
    """
    \$ANNOTSV/bin/AnnotSV \\
        -SVinputFile ${sv_vcf} \\
        -txFile ${gene_transcripts} \\
        -genomeBuild ${genome_version} \\
        -outputFile ${prefix}.annotated.tsv \\
        -candidateGenesFile ${candidate_genes} \\
        -candidateSnvIndelFiles ${snv_indel_vcf} \\
        -vcf 1
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        annotsv: \$(echo \$(AnnotSV -help 2>&1 | head -n1 | sed 's/^AnnotSV //'))
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.patient_id}"
    """
    touch ${prefix}.annotated.vcf
    touch ${prefix}.annotated.tsv
    touch ${prefix}.annotated.variantconvert.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        annotsv: \$(echo \$(AnnotSV -help 2>&1 | head -n1 | sed 's/^AnnotSV //'))
    END_VERSIONS
    """
}
