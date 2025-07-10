process SVTK_STANDARDIZE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/svtk:0.0.20190615--py39hbcbf7aa_6':
        'quay.io/biocontainers/svtk:0.0.20190615-6' }"

    input:
    tuple val(meta) , path(vcf_delly)
    tuple val(meta1), path(vcf_svaba)
    tuple val(meta2), path(vcf_manta)
    tuple val(meta3), path(vcf_tiddit)
    tuple val(meta4), path(vcf_gridss)
    path fasta_fai

    output:
    tuple val(meta), path("*_delly.std.vcf.gz") , emit: delly_std_vcf
    tuple val(meta), path("*_svaba.std.vcf.gz") , emit: svaba_std_vcf
    tuple val(meta), path("*_manta.std.vcf.gz") , emit: manta_std_vcf
    tuple val(meta), path("*_tiddit.std.vcf.gz"), emit: tiddit_std_vcf
    tuple val(meta), path("*_gridss.std.vcf.gz"), emit: gridss_std_vcf
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix      = task.ext.prefix           ?: "${meta.id}"
    def caller      = args.caller               ?: 'delly'
    def VERSION     = '0.0.20190615'            // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def contigs     = '--contigs ${fasta_fai}'
    """
    svtk standardize \\
        ${contigs} \\
        ${vcf_delly} \\
        ${prefix}_delly.std.vcf.gz

    svtk standardize \\
        ${contigs} \\
        ${vcf_svaba} \\
        ${prefix}_svaba.std.vcf.gz

    svtk standardize \\
        ${contigs} \\
        ${vcf_manta} \\
        ${prefix}_manta.std.vcf.gz

    svtk standardize \\
        ${contigs} \\
        ${vcf_tiddit} \\
        ${prefix}_tiddit.std.vcf.gz

    svtk standardize \\
        ${contigs} \\
        ${vcf_gridss} \\
        ${prefix}_gridss.std.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svtk: ${VERSION}
    END_VERSIONS
    """
}
