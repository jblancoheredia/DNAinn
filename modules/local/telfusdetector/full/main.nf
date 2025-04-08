process TELFUSDETECTOR_FULL {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/telfusdetector:1.0.0' :
        'blancojmskcc/telfusdetector:1.0.0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    val(genome)
    val(purity)
    val(panel)

    output:
    tuple val(meta), path("*_rates.tsv"), emit: tsv
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def hg = (genome == 'GRCh38') ? 'Hg38' : (genome == 'GRCh37') ? 'Hg19' : null
    if (hg == null) {
        error "Unsupported genome: ${genome}. Only GRCh38 and GRCh37 are supported."
    }
    """
    TMP=\${TMPDIR:-/tmp}

    cp "${workflow.projectDir}/bin/TelFusCall.py" .

    cp "${workflow.projectDir}/bin/TelFusRate.py" .

    python TelFusCall.py \\
        --bam ${bam} \\
        --genome ${hg} \\
        --panel ${panel} \\
        --outfolder . \\
        --sample ${prefix} \\
        --tmpfolder \${TMP} \\
        --threads ${task.cpus} \\

    python TelFusRate.py \\
        --fusion_file ${prefix}.summary_fusions.pass.tsv \\
        --variables Filter Type Orientation_full_pair \\
        --purity ${purity} \\
        --outfile ${prefix}.telomere_fusion_rates.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        telfusdetector: "1.0.0"
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.telomere_fusion_rates.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        telfusdetector: "1.0.0"
    END_VERSIONS
    """
}
