process GRAPHTYPER_GENOTYPE_SV {
    tag "$meta.patient_id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/graphtyper:2.7.7':
        'blancojmskcc/graphtyper:2.7.7' }"

    input:
    tuple val(meta) , 
          val(meta2), path(bam), path(bai), 
          val(meta3), path(vcf), path(tbi)
    tuple val(meta4), path(ref)
    tuple val(meta5), path(ref_fai)
    path region_file

    output:
    tuple val(meta), path("*sv_results_large.vcf.gz")    , emit: vcf
    tuple val(meta), path("*sv_results_large.vcf.gz.tbi"), emit: tbi
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def bam_path_text = bam.sort().join('\\n')
    def prefix = task.ext.prefix ?: "${meta.patient_id}"
    def region_text   = region_file.size() > 0 ? "--region_file ${region_file}" : ""
    if (region_file.size() == 0 && ! args.contains("region")) {
        error "GRAPHTYPER_GENOTYPE_SV requires either a region file or a region specified using '--region' in ext.args"
    }
    """
    printf "${bam_path_text}" > bam_list.txt
    
    graphtyper \\
        genotype_sv \\
        ${ref} \\
        ${vcf} \\
        $args \\
        --sams bam_list.txt \\
        --threads ${task.cpus} \\
        ${region_text}

    echo chr{1..22} chrX | tr ' ' '\\n' | while read chrom; do if [[ ! -d sv_results/\${chrom} ]]; then continue; fi; find sv_results/\${chrom} -name \"*.vcf.gz\" | sort; done > vcf_file_list

    bcftools concat --naive --file-list vcf_file_list -Oz -o ${prefix}_sv_results_large.vcf.gz

    tabix -p vcf ${prefix}_sv_results_large.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        graphtyper: \$(graphtyper --help | tail -n 1 | sed 's/^   //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p results/test
    echo | gzip > results/test/region1.vcf.gz
    echo | gzip > results/test/region2.vcf.gz
    touch results/test/region1.vcf.gz.tbi
    touch results/test/region2.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        graphtyper: \$(graphtyper --help | tail -n 1 | sed 's/^   //')
    END_VERSIONS
    """

}
