process MANTA {
    tag "$meta.patient"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/manta:1.6.0--h9ee0642_1' :
        'quay.io/biocontainers/manta:1.6.0--h9ee0642_1' }"

    input:
    tuple val(meta) , path(tumour_bam), path(tumour_bai), path(normal_bam), path(normal_bai)
    tuple val(meta2), path(target_bed) 
    tuple val(meta3), path(bed_index)
    tuple val(meta4), path(fasta)
    tuple val(meta5), path(fai)
    path(config)

    output:
    tuple val(meta), path("*.tsv",  optional: true), emit: metrics_tsv
    tuple val(meta), path("*.txt",  optional: true), emit: metrics_txt
    tuple val(meta), path("*.manta.unfiltered.vcf"), emit: vcf
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.patient}"
    def options_manta = target_bed ? "--callRegions ${target_bed}" : ""
    """
    configManta.py \\
        --tumorBam ${tumour_bam} \\
        --reference ${fasta} \\
        --runDir manta_tumour \\
        ${options_manta}    
    python manta_tumour/runWorkflow.py -m local -j ${task.cpus}

    mv manta_tumour/results/stats/svCandidateGenerationStats.tsv ${prefix}.manta_tumour.svCandidateGenerationStats.tsv
    mv manta_tumour/results/stats/alignmentStatsSummary.txt ${prefix}.manta_tumour.alignmentStatsSummary.txt
    mv manta_tumour/results/stats/svLocusGraphStats.tsv ${prefix}.manta_tumour.svLocusGraphStats.tsv
    mv manta_tumour/results/variants/candidateSV.vcf.gz.tbi ${prefix}.manta.unfiltered.vcf.gz.tbi
    mv manta_tumour/results/variants/candidateSV.vcf.gz ${prefix}.manta.unfiltered.vcf.gz

    gunzip ${prefix}.manta.unfiltered.vcf.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        manta: \$( configManta.py --version )
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.patient}"
    """
    touch *.manta.unfiltered.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        manta: \$( configManta.py --version )
    END_VERSIONS
    """
}
