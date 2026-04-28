process LINE_PROBE_RAW {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:03dc1d2818d9de56938078b8b78b82d967c1f820-0' :
        'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:03dc1d2818d9de56938078b8b78b82d967c1f820-0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path probe_fasta
    path bwa_dir

    output:
    tuple val(meta), path("*.line_probe.raw.bam")        , emit: bam
    tuple val(meta), path("*.line_probe.raw.bam.bai")    , emit: bai
    tuple val(meta), path("*.line_probe.depth.raw.tsv")  , emit: depth
    tuple val(meta), path("*.line_probe.summary.raw.tsv"), emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fasta_name = probe_fasta.getName()

    """
    BWA_INDEX_PREFIX=\$(find -L ./ -name "*.amb" | sed 's/.amb//')

    samtools fastq -@ ${task.cpus} ${prefix}.sorted.bam \\
        | bwa mem -t ${task.cpus} -B 3 -K 100000000 -Y -M -R ${meta.read_group} \$BWA_INDEX_PREFIX - \\
        | samtools sort -@ ${task.cpus} -o ${prefix}.line_probe.raw.bam

    samtools index ${prefix}.line_probe.raw.bam

    samtools depth \\
        -aa \\
        ${prefix}.line_probe.raw.bam \\
        > ${prefix}.line_probe.depth.raw.tsv

    total_records=\$(samtools view -c ${prefix}.line_probe.raw.bam)
    mapped_records=\$(samtools view -c -F 4 ${prefix}.line_probe.raw.bam)
    mapped_primary_records=\$(samtools view -c -F 260 ${prefix}.line_probe.raw.bam)
    proper_pairs=\$(samtools view -c -f 2 ${prefix}.line_probe.raw.bam)
    mean_depth=\$(awk '{sum+=\$3; n++} END{if (n > 0) print sum/n; else print 0}' ${prefix}.line_probe.depth.raw.tsv)
    bases_covered_1x=\$(awk '\$3 >= 1 {n++} END{print n+0}' ${prefix}.line_probe.depth.raw.tsv)
    bases_covered_10x=\$(awk '\$3 >= 10 {n++} END{print n+0}' ${prefix}.line_probe.depth.raw.tsv)
    probe_length=\$(awk 'BEGIN{n=0} !/^>/ {n+=length(\$0)} END{print n}' ${probe_fasta})

    {
        echo -e "sample\\tprobe_fasta\\ttotal_records\\tmapped_records\\tmapped_primary_records\\tproper_pairs\\tmean_depth\\tprobe_length\\tbases_covered_1x\\tbases_covered_10x"
        echo -e "${prefix}\\t${fasta_name}\\t\${total_records}\\t\${mapped_records}\\t\${mapped_primary_records}\\t\${proper_pairs}\\t\${mean_depth}\\t\${probe_length}\\t\${bases_covered_1x}\\t\${bases_covered_10x}"
    } > ${prefix}.line_probe.summary.raw.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | tail -n 1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
	touch ${prefix}..line_probe.raw.bam     
	touch ${prefix}..line_probe.raw.bam.bai   
	touch ${prefix}..line_probe.depth.raw.tsv 
	touch ${prefix}..line_probe.summary.raw.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | tail -n 1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

}

process LINE_PROBE_CON {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:03dc1d2818d9de56938078b8b78b82d967c1f820-0' :
        'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:03dc1d2818d9de56938078b8b78b82d967c1f820-0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path probe_fasta
    path bwa_dir

    output:
    tuple val(meta), path("*.line_probe.con.bam")        , emit: bam
    tuple val(meta), path("*.line_probe.con.bam.bai")    , emit: bai
    tuple val(meta), path("*.line_probe.depth.con.tsv")  , emit: depth
    tuple val(meta), path("*.line_probe.summary.con.tsv"), emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fasta_name = probe_fasta.getName()

    """
    BWA_INDEX_PREFIX=\$(find -L ./ -name "*.amb" | sed 's/.amb//')

    samtools fastq -@ ${task.cpus} ${prefix}.sorted.bam \\
        | bwa mem -t ${task.cpus} -B 3 -K 100000000 -Y -M -R ${meta.read_group} \$BWA_INDEX_PREFIX - \\
        | samtools sort -@ ${task.cpus} -o ${prefix}.line_probe.con.bam

    samtools index ${prefix}.line_probe.con.bam

    samtools depth \\
        -aa \\
        ${prefix}.line_probe.con.bam \\
        > ${prefix}.line_probe.depth.con.tsv

    total_records=\$(samtools view -c ${prefix}.line_probe.con.bam)
    mapped_records=\$(samtools view -c -F 4 ${prefix}.line_probe.con.bam)
    mapped_primary_records=\$(samtools view -c -F 260 ${prefix}.line_probe.con.bam)
    proper_pairs=\$(samtools view -c -f 2 ${prefix}.line_probe.con.bam)
    mean_depth=\$(awk '{sum+=\$3; n++} END{if (n > 0) print sum/n; else print 0}' ${prefix}.line_probe.depth.con.tsv)
    bases_covered_1x=\$(awk '\$3 >= 1 {n++} END{print n+0}' ${prefix}.line_probe.depth.con.tsv)
    bases_covered_10x=\$(awk '\$3 >= 10 {n++} END{print n+0}' ${prefix}.line_probe.depth.con.tsv)
    probe_length=\$(awk 'BEGIN{n=0} !/^>/ {n+=length(\$0)} END{print n}' ${probe_fasta})

    {
        echo -e "sample\\tprobe_fasta\\ttotal_records\\tmapped_records\\tmapped_primary_records\\tproper_pairs\\tmean_depth\\tprobe_length\\tbases_covered_1x\\tbases_covered_10x"
        echo -e "${prefix}\\t${fasta_name}\\t\${total_records}\\t\${mapped_records}\\t\${mapped_primary_records}\\t\${proper_pairs}\\t\${mean_depth}\\t\${probe_length}\\t\${bases_covered_1x}\\t\${bases_covered_10x}"
    } > ${prefix}.line_probe.summary.con.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | tail -n 1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
	touch ${prefix}..line_probe.con.bam     
	touch ${prefix}..line_probe.con.bam.bai   
	touch ${prefix}..line_probe.depth.con.tsv 
	touch ${prefix}..line_probe.summary.con.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | tail -n 1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

}

process LINE_PROBE_DUP {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:03dc1d2818d9de56938078b8b78b82d967c1f820-0' :
        'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:03dc1d2818d9de56938078b8b78b82d967c1f820-0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path probe_fasta
    path bwa_dir

    output:
    tuple val(meta), path("*.line_probe.dup.bam")        , emit: bam
    tuple val(meta), path("*.line_probe.dup.bam.bai")    , emit: bai
    tuple val(meta), path("*.line_probe.depth.dup.tsv")  , emit: depth
    tuple val(meta), path("*.line_probe.summary.dup.tsv"), emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fasta_name = probe_fasta.getName()

    """
    BWA_INDEX_PREFIX=\$(find -L ./ -name "*.amb" | sed 's/.amb//')

    samtools fastq -@ ${task.cpus} ${prefix}.sorted.bam \\
        | bwa mem -t ${task.cpus} -B 3 -K 100000000 -Y -M -R ${meta.read_group} \$BWA_INDEX_PREFIX - \\
        | samtools sort -@ ${task.cpus} -o ${prefix}.line_probe.dup.bam

    samtools index ${prefix}.line_probe.dup.bam

    samtools depth \\
        -aa \\
        ${prefix}.line_probe.dup.bam \\
        > ${prefix}.line_probe.depth.dup.tsv

    total_records=\$(samtools view -c ${prefix}.line_probe.dup.bam)
    mapped_records=\$(samtools view -c -F 4 ${prefix}.line_probe.dup.bam)
    mapped_primary_records=\$(samtools view -c -F 260 ${prefix}.line_probe.dup.bam)
    proper_pairs=\$(samtools view -c -f 2 ${prefix}.line_probe.dup.bam)
    mean_depth=\$(awk '{sum+=\$3; n++} END{if (n > 0) print sum/n; else print 0}' ${prefix}.line_probe.depth.dup.tsv)
    bases_covered_1x=\$(awk '\$3 >= 1 {n++} END{print n+0}' ${prefix}.line_probe.depth.dup.tsv)
    bases_covered_10x=\$(awk '\$3 >= 10 {n++} END{print n+0}' ${prefix}.line_probe.depth.dup.tsv)
    probe_length=\$(awk 'BEGIN{n=0} !/^>/ {n+=length(\$0)} END{print n}' ${probe_fasta})

    {
        echo -e "sample\\tprobe_fasta\\ttotal_records\\tmapped_records\\tmapped_primary_records\\tproper_pairs\\tmean_depth\\tprobe_length\\tbases_covered_1x\\tbases_covered_10x"
        echo -e "${prefix}\\t${fasta_name}\\t\${total_records}\\t\${mapped_records}\\t\${mapped_primary_records}\\t\${proper_pairs}\\t\${mean_depth}\\t\${probe_length}\\t\${bases_covered_1x}\\t\${bases_covered_10x}"
    } > ${prefix}.line_probe.summary.dup.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | tail -n 1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
	touch ${prefix}..line_probe.dup.bam     
	touch ${prefix}..line_probe.dup.bam.bai   
	touch ${prefix}..line_probe.depth.dup.tsv 
	touch ${prefix}..line_probe.summary.dup.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | tail -n 1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

}

process LINE_PROBE_SIM {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:03dc1d2818d9de56938078b8b78b82d967c1f820-0' :
        'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:03dc1d2818d9de56938078b8b78b82d967c1f820-0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path probe_fasta
    path bwa_dir

    output:
    tuple val(meta), path("*.line_probe.sim.bam")        , emit: bam
    tuple val(meta), path("*.line_probe.sim.bam.bai")    , emit: bai
    tuple val(meta), path("*.line_probe.depth.sim.tsv")  , emit: depth
    tuple val(meta), path("*.line_probe.summary.sim.tsv"), emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fasta_name = probe_fasta.getName()

    """
    BWA_INDEX_PREFIX=\$(find -L ./ -name "*.amb" | sed 's/.amb//')

    samtools fastq -@ ${task.cpus} ${prefix}.sorted.bam \\
        | bwa mem -t ${task.cpus} -B 3 -K 100000000 -Y -M -R ${meta.read_group} \$BWA_INDEX_PREFIX - \\
        | samtools sort -@ ${task.cpus} -o ${prefix}.line_probe.sim.bam

    samtools index ${prefix}.line_probe.sim.bam

    samtools depth \\
        -aa \\
        ${prefix}.line_probe.sim.bam \\
        > ${prefix}.line_probe.depth.sim.tsv

    total_records=\$(samtools view -c ${prefix}.line_probe.sim.bam)
    mapped_records=\$(samtools view -c -F 4 ${prefix}.line_probe.sim.bam)
    mapped_primary_records=\$(samtools view -c -F 260 ${prefix}.line_probe.sim.bam)
    proper_pairs=\$(samtools view -c -f 2 ${prefix}.line_probe.sim.bam)
    mean_depth=\$(awk '{sum+=\$3; n++} END{if (n > 0) print sum/n; else print 0}' ${prefix}.line_probe.depth.sim.tsv)
    bases_covered_1x=\$(awk '\$3 >= 1 {n++} END{print n+0}' ${prefix}.line_probe.depth.sim.tsv)
    bases_covered_10x=\$(awk '\$3 >= 10 {n++} END{print n+0}' ${prefix}.line_probe.depth.sim.tsv)
    probe_length=\$(awk 'BEGIN{n=0} !/^>/ {n+=length(\$0)} END{print n}' ${probe_fasta})

    {
        echo -e "sample\\tprobe_fasta\\ttotal_records\\tmapped_records\\tmapped_primary_records\\tproper_pairs\\tmean_depth\\tprobe_length\\tbases_covered_1x\\tbases_covered_10x"
        echo -e "${prefix}\\t${fasta_name}\\t\${total_records}\\t\${mapped_records}\\t\${mapped_primary_records}\\t\${proper_pairs}\\t\${mean_depth}\\t\${probe_length}\\t\${bases_covered_1x}\\t\${bases_covered_10x}"
    } > ${prefix}.line_probe.summary.sim.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | tail -n 1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
	touch ${prefix}..line_probe.sim.bam     
	touch ${prefix}..line_probe.sim.bam.bai   
	touch ${prefix}..line_probe.depth.sim.tsv 
	touch ${prefix}..line_probe.summary.sim.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | tail -n 1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

}
