process SAMTOOLS_SORT_INDEX {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'quay.io/biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta) , path(bam)
    tuple val(meta2), path(fasta)
    path(fai)

    output:
    tuple val(meta), path("*_sort.bam"), path("*.bai"), emit: bam_bai
    tuple val(meta), path("*_sort.bam"),                emit: bam
    tuple val(meta), path("*.bai"),                     emit: bai
    path  "versions.yml",                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_sort"
    def extension = args.contains("--output-fmt sam") ? "sam" :
                    args.contains("--output-fmt cram") ? "cram" :
                    "bam"
    def reference = fasta ? "--reference ${fasta}" : ""
    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"

    """
    samtools cat \\
        ${bam} \\
    | \\
    samtools sort \\
        $args \\
        -T ${prefix} \\
        --threads $task.cpus \\
        ${reference} \\
        -o ${prefix}.${extension} \\
        -
    
    samtools index \\
        -@ ${task.cpus-1} \\
        $args \\
        ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_sort"
    def extension = args.contains("--output-fmt sam") ? "sam" :
                    args.contains("--output-fmt cram") ? "cram" :
                    "bam"
    """
    touch ${prefix}.${extension}
    if [ "${extension}" == "bam" ];
    then
        touch ${prefix}.${extension}.csi
    elif [ "${extension}" == "cram" ];
    then
        touch ${prefix}.${extension}.crai
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}

process SAMTOOLS_SORT_INDEX_RAW {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'quay.io/biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta) , path(bam)
    tuple val(meta2), path(fasta)
    path(fai)

    output:
    tuple val(meta), path("*_sort.bam"), path("*.bai"), emit: bam_bai
    tuple val(meta), path("*_sort.bam"),                emit: bam
    tuple val(meta), path("*.bai"),                     emit: bai
    path  "versions.yml",                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_sort"
    def extension = args.contains("--output-fmt sam") ? "sam" :
                    args.contains("--output-fmt cram") ? "cram" :
                    "bam"
    def reference = fasta ? "--reference ${fasta}" : ""
    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"

    """
    samtools cat \\
        ${bam} \\
    | \\
    samtools sort \\
        $args \\
        -T ${prefix} \\
        --threads $task.cpus \\
        ${reference} \\
        -o ${prefix}.${extension} \\
        -
    
    samtools index \\
        -@ ${task.cpus-1} \\
        $args \\
        ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_sort"
    def extension = args.contains("--output-fmt sam") ? "sam" :
                    args.contains("--output-fmt cram") ? "cram" :
                    "bam"
    """
    touch ${prefix}.${extension}
    if [ "${extension}" == "bam" ];
    then
        touch ${prefix}.${extension}.csi
    elif [ "${extension}" == "cram" ];
    then
        touch ${prefix}.${extension}.crai
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}

process SAMTOOLS_SORT_INDEX_FIN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'quay.io/biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta) , path(bam)
    tuple val(meta3), path(fai)
    tuple val(meta2), path(fasta)
    tuple val(meta0), path(duplex_bam)
    tuple val(meta1), path(simplex_bam)

    output:
    tuple val(meta), path("*_sort.bam")         , path("*_sort.bam.bai")        , emit: bam_bai
    tuple val(meta), path("*_sort.bam")                                         , emit: bam
    tuple val(meta), path("*_sort.bam.bai")                                     , emit: bai
    tuple val(meta), path("*_sort.duplex.bam")  , path("*_sort.duplex.bam.bai") , emit: bam_duplex
    tuple val(meta), path("*_sort.simplex.bam") , path("*_sort.simplex.bam.bai"), emit: bam_simplex
    path  "versions.yml"                                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_sort"
    def extension = args.contains("--output-fmt sam") ? "sam" :
                    args.contains("--output-fmt cram") ? "cram" :
                    "bam"
    def reference = fasta ? "--reference ${fasta}" : ""
    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"

    """
    samtools cat \\
        ${bam} \\
    | \\
    samtools sort \\
        $args \\
        -T ${prefix} \\
        --threads $task.cpus \\
        ${reference} \\
        -o ${prefix}.${extension} \\
        -
    
    samtools index \\
        -@ ${task.cpus-1} \\
        $args \\
        ${prefix}.${extension}

    samtools cat \\
        ${duplex_bam} \\
    | \\
    samtools sort \\
        $args \\
        -T ${prefix} \\
        --threads $task.cpus \\
        ${reference} \\
        -o ${prefix}.duplex.${extension} \\
        -

    samtools index \\
        -@ ${task.cpus-1} \\
        $args \\
        ${prefix}.duplex.${extension}

    samtools cat \\
        ${simplex_bam} \\
    | \\
    samtools sort \\
        $args \\
        -T ${prefix} \\
        --threads $task.cpus \\
        ${reference} \\
        -o ${prefix}.simplex.${extension} \\
        -

    samtools index \\
        -@ ${task.cpus-1} \\
        $args \\
        ${prefix}.simplex.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_sort"
    def extension = args.contains("--output-fmt sam") ? "sam" :
                    args.contains("--output-fmt cram") ? "cram" :
                    "bam"
    """
    touch ${prefix}.${extension}
    touch ${prefix}.duplex.${extension}
    touch ${prefix}.simplex.${extension}
    if [ "${extension}" == "bam" ];
    then
        touch ${prefix}.${extension}.csi
    elif [ "${extension}" == "cram" ];
    then
        touch ${prefix}.${extension}.crai
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
