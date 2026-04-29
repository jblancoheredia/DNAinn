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
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-B 3 -K 100000000 -Y -M'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fasta_name = probe_fasta.getName()

    """
    BWA_INDEX_PREFIX=\$(find -L ./ -name "*.amb" | sed 's/.amb//')

    samtools collate -u -O ${bam} ${prefix}.collate \\
        | samtools fastq \\
            -n \\
            -F 0x900 \\
            -1 /dev/stdout \\
            -2 /dev/stdout \\
            -0 /dev/null \\
            -s /dev/null \\
            - \\
        | bwa mem -p -t ${task.cpus} ${args} -R ${meta.read_group} \$BWA_INDEX_PREFIX - \\
        | samtools sort -@ ${task.cpus} -o ${prefix}.line_probe.raw.bam

    samtools index ${prefix}.line_probe.raw.bam

    samtools depth -aa ${prefix}.line_probe.raw.bam > ${prefix}.line_probe.depth.raw.tsv

    total_records=\$(samtools view -c ${prefix}.line_probe.raw.bam)
    mapped_records=\$(samtools view -c -F 4 ${prefix}.line_probe.raw.bam)
    mapped_primary_records=\$(samtools view -c -F 2308 ${prefix}.line_probe.raw.bam)
    proper_pairs=\$(samtools view -c -f 2 ${prefix}.line_probe.raw.bam)
    mean_depth=\$(awk '{sum+=\$3; n++} END{if (n > 0) print sum/n; else print 0}' ${prefix}.line_probe.depth.raw.tsv)
    bases_covered_1x=\$(awk '\$3 >= 1 {n++} END{print n+0}' ${prefix}.line_probe.depth.raw.tsv)
    bases_covered_10x=\$(awk '\$3 >= 10 {n++} END{print n+0}' ${prefix}.line_probe.depth.raw.tsv)
    probe_length=\$(awk 'BEGIN{n=0} !/^>/ {n+=length(\$0)} END{print n}' ${probe_fasta})

    samtools view -F 2308 ${prefix}.line_probe.raw.bam \\
        | awk '
            BEGIN {
                OFS="\\t"
            }
            {
                cigar=\$6
                pos=\$4
                aligned=0
                soft=0
                ref_span=0
                c=cigar
                while (match(c, /^[0-9]+[MIDNSHP=X]/)) {
                    len=substr(c, RSTART, RLENGTH - 1) + 0
                    op=substr(c, RSTART + RLENGTH - 1, 1)
                    if (op == "M" || op == "=" || op == "X") aligned += len
                    if (op == "S") soft += len
                    if (op == "M" || op == "D" || op == "N" || op == "=" || op == "X") ref_span += len
                    c=substr(c, RLENGTH + 1)
                }

                end=pos + ref_span - 1
                records++
                total_aligned += aligned
                total_soft += soft
                if (aligned >= 10) aligned_ge_10++
                if (aligned >= 20) aligned_ge_20++
                if (aligned >= 30) aligned_ge_30++
                if (aligned >= 50) aligned_ge_50++
                if (aligned >= 80) aligned_ge_80++
                if (aligned >= 100) aligned_ge_100++
                if (pos <= 1 && end >= 9) cover_1_9++
                if (pos <= 1 && end >= 25) cover_1_25++
                if (pos <= 1 && end >= 31) cover_1_31++
                if (pos <= 1 && end >= 40) cover_1_40++
                if (pos >= 10) start_ge_10++
            }
            END {
                mean_aligned = records > 0 ? total_aligned / records : 0
                mean_soft = records > 0 ? total_soft / records : 0
                print total_aligned+0, total_soft+0, mean_aligned, mean_soft, aligned_ge_10+0, aligned_ge_20+0, aligned_ge_30+0, aligned_ge_50+0, aligned_ge_80+0, aligned_ge_100+0, cover_1_9+0, cover_1_25+0, cover_1_31+0, cover_1_40+0, start_ge_10+0
            }
        ' > probe_alignment_metrics.tsv

    read total_aligned_bases total_softclipped_bases mean_aligned_bases_per_record mean_softclipped_bases_per_record records_aligned_ge_10bp records_aligned_ge_20bp records_aligned_ge_30bp records_aligned_ge_50bp records_aligned_ge_80bp records_aligned_ge_100bp records_covering_probe_1_9 records_covering_probe_1_25 records_covering_probe_1_31 records_covering_probe_1_40 records_starting_at_or_after_10 < probe_alignment_metrics.tsv

    {
        echo -e "sample\\tprobe_fasta\\ttotal_records\\tmapped_records\\tmapped_primary_records\\tproper_pairs\\tmean_depth\\tprobe_length\\tbases_covered_1x\\tbases_covered_10x\\ttotal_aligned_bases\\ttotal_softclipped_bases\\tmean_aligned_bases_per_record\\tmean_softclipped_bases_per_record\\trecords_aligned_ge_10bp\\trecords_aligned_ge_20bp\\trecords_aligned_ge_30bp\\trecords_aligned_ge_50bp\\trecords_aligned_ge_80bp\\trecords_aligned_ge_100bp\\trecords_covering_probe_1_9\\trecords_covering_probe_1_25\\trecords_covering_probe_1_31\\trecords_covering_probe_1_40\\trecords_starting_at_or_after_10"
        echo -e "${prefix}\\t${fasta_name}\\t\${total_records}\\t\${mapped_records}\\t\${mapped_primary_records}\\t\${proper_pairs}\\t\${mean_depth}\\t\${probe_length}\\t\${bases_covered_1x}\\t\${bases_covered_10x}\\t\${total_aligned_bases}\\t\${total_softclipped_bases}\\t\${mean_aligned_bases_per_record}\\t\${mean_softclipped_bases_per_record}\\t\${records_aligned_ge_10bp}\\t\${records_aligned_ge_20bp}\\t\${records_aligned_ge_30bp}\\t\${records_aligned_ge_50bp}\\t\${records_aligned_ge_80bp}\\t\${records_aligned_ge_100bp}\\t\${records_covering_probe_1_9}\\t\${records_covering_probe_1_25}\\t\${records_covering_probe_1_31}\\t\${records_covering_probe_1_40}\\t\${records_starting_at_or_after_10}"
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
    touch ${prefix}.line_probe.raw.bam
    touch ${prefix}.line_probe.raw.bam.bai
    touch ${prefix}.line_probe.depth.raw.tsv
    touch ${prefix}.line_probe.summary.raw.tsv

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
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-B 3 -K 100000000 -Y -M'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fasta_name = probe_fasta.getName()

    """
    BWA_INDEX_PREFIX=\$(find -L ./ -name "*.amb" | sed 's/.amb//')

    samtools collate -u -O ${bam} ${prefix}.collate \\
        | samtools fastq \\
            -n \\
            -F 0x900 \\
            -1 /dev/stdout \\
            -2 /dev/stdout \\
            -0 /dev/null \\
            -s /dev/null \\
            - \\
        | bwa mem -p -t ${task.cpus} ${args} -R ${meta.read_group} \$BWA_INDEX_PREFIX - \\
        | samtools sort -@ ${task.cpus} -o ${prefix}.line_probe.con.bam

    samtools index ${prefix}.line_probe.con.bam

    samtools depth -aa ${prefix}.line_probe.con.bam > ${prefix}.line_probe.depth.con.tsv

    total_records=\$(samtools view -c ${prefix}.line_probe.con.bam)
    mapped_records=\$(samtools view -c -F 4 ${prefix}.line_probe.con.bam)
    mapped_primary_records=\$(samtools view -c -F 2308 ${prefix}.line_probe.con.bam)
    proper_pairs=\$(samtools view -c -f 2 ${prefix}.line_probe.con.bam)
    mean_depth=\$(awk '{sum+=\$3; n++} END{if (n > 0) print sum/n; else print 0}' ${prefix}.line_probe.depth.con.tsv)
    bases_covered_1x=\$(awk '\$3 >= 1 {n++} END{print n+0}' ${prefix}.line_probe.depth.con.tsv)
    bases_covered_10x=\$(awk '\$3 >= 10 {n++} END{print n+0}' ${prefix}.line_probe.depth.con.tsv)
    probe_length=\$(awk 'BEGIN{n=0} !/^>/ {n+=length(\$0)} END{print n}' ${probe_fasta})

    samtools view -F 2308 ${prefix}.line_probe.con.bam \\
        | awk '
            BEGIN {
                OFS="\\t"
            }
            {
                cigar=\$6
                pos=\$4
                aligned=0
                soft=0
                ref_span=0
                c=cigar
                while (match(c, /^[0-9]+[MIDNSHP=X]/)) {
                    len=substr(c, RSTART, RLENGTH - 1) + 0
                    op=substr(c, RSTART + RLENGTH - 1, 1)
                    if (op == "M" || op == "=" || op == "X") aligned += len
                    if (op == "S") soft += len
                    if (op == "M" || op == "D" || op == "N" || op == "=" || op == "X") ref_span += len
                    c=substr(c, RLENGTH + 1)
                }

                end=pos + ref_span - 1
                records++
                total_aligned += aligned
                total_soft += soft
                if (aligned >= 10) aligned_ge_10++
                if (aligned >= 20) aligned_ge_20++
                if (aligned >= 30) aligned_ge_30++
                if (aligned >= 50) aligned_ge_50++
                if (aligned >= 80) aligned_ge_80++
                if (aligned >= 100) aligned_ge_100++
                if (pos <= 1 && end >= 9) cover_1_9++
                if (pos <= 1 && end >= 25) cover_1_25++
                if (pos <= 1 && end >= 31) cover_1_31++
                if (pos <= 1 && end >= 40) cover_1_40++
                if (pos >= 10) start_ge_10++
            }
            END {
                mean_aligned = records > 0 ? total_aligned / records : 0
                mean_soft = records > 0 ? total_soft / records : 0
                print total_aligned+0, total_soft+0, mean_aligned, mean_soft, aligned_ge_10+0, aligned_ge_20+0, aligned_ge_30+0, aligned_ge_50+0, aligned_ge_80+0, aligned_ge_100+0, cover_1_9+0, cover_1_25+0, cover_1_31+0, cover_1_40+0, start_ge_10+0
            }
        ' > probe_alignment_metrics.tsv

    read total_aligned_bases total_softclipped_bases mean_aligned_bases_per_record mean_softclipped_bases_per_record records_aligned_ge_10bp records_aligned_ge_20bp records_aligned_ge_30bp records_aligned_ge_50bp records_aligned_ge_80bp records_aligned_ge_100bp records_covering_probe_1_9 records_covering_probe_1_25 records_covering_probe_1_31 records_covering_probe_1_40 records_starting_at_or_after_10 < probe_alignment_metrics.tsv

    {
        echo -e "sample\\tprobe_fasta\\ttotal_records\\tmapped_records\\tmapped_primary_records\\tproper_pairs\\tmean_depth\\tprobe_length\\tbases_covered_1x\\tbases_covered_10x\\ttotal_aligned_bases\\ttotal_softclipped_bases\\tmean_aligned_bases_per_record\\tmean_softclipped_bases_per_record\\trecords_aligned_ge_10bp\\trecords_aligned_ge_20bp\\trecords_aligned_ge_30bp\\trecords_aligned_ge_50bp\\trecords_aligned_ge_80bp\\trecords_aligned_ge_100bp\\trecords_covering_probe_1_9\\trecords_covering_probe_1_25\\trecords_covering_probe_1_31\\trecords_covering_probe_1_40\\trecords_starting_at_or_after_10"
        echo -e "${prefix}\\t${fasta_name}\\t\${total_records}\\t\${mapped_records}\\t\${mapped_primary_records}\\t\${proper_pairs}\\t\${mean_depth}\\t\${probe_length}\\t\${bases_covered_1x}\\t\${bases_covered_10x}\\t\${total_aligned_bases}\\t\${total_softclipped_bases}\\t\${mean_aligned_bases_per_record}\\t\${mean_softclipped_bases_per_record}\\t\${records_aligned_ge_10bp}\\t\${records_aligned_ge_20bp}\\t\${records_aligned_ge_30bp}\\t\${records_aligned_ge_50bp}\\t\${records_aligned_ge_80bp}\\t\${records_aligned_ge_100bp}\\t\${records_covering_probe_1_9}\\t\${records_covering_probe_1_25}\\t\${records_covering_probe_1_31}\\t\${records_covering_probe_1_40}\\t\${records_starting_at_or_after_10}"
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
    touch ${prefix}.line_probe.con.bam
    touch ${prefix}.line_probe.con.bam.bai
    touch ${prefix}.line_probe.depth.con.tsv
    touch ${prefix}.line_probe.summary.con.tsv

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
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-B 3 -K 100000000 -Y -M'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fasta_name = probe_fasta.getName()

    """
    BWA_INDEX_PREFIX=\$(find -L ./ -name "*.amb" | sed 's/.amb//')

    samtools collate -u -O ${bam} ${prefix}.collate \\
        | samtools fastq \\
            -n \\
            -F 0x900 \\
            -1 /dev/stdout \\
            -2 /dev/stdout \\
            -0 /dev/null \\
            -s /dev/null \\
            - \\
        | bwa mem -p -t ${task.cpus} ${args} -R ${meta.read_group} \$BWA_INDEX_PREFIX - \\
        | samtools sort -@ ${task.cpus} -o ${prefix}.line_probe.dup.bam

    samtools index ${prefix}.line_probe.dup.bam

    samtools depth -aa ${prefix}.line_probe.dup.bam > ${prefix}.line_probe.depth.dup.tsv

    total_records=\$(samtools view -c ${prefix}.line_probe.dup.bam)
    mapped_records=\$(samtools view -c -F 4 ${prefix}.line_probe.dup.bam)
    mapped_primary_records=\$(samtools view -c -F 2308 ${prefix}.line_probe.dup.bam)
    proper_pairs=\$(samtools view -c -f 2 ${prefix}.line_probe.dup.bam)
    mean_depth=\$(awk '{sum+=\$3; n++} END{if (n > 0) print sum/n; else print 0}' ${prefix}.line_probe.depth.dup.tsv)
    bases_covered_1x=\$(awk '\$3 >= 1 {n++} END{print n+0}' ${prefix}.line_probe.depth.dup.tsv)
    bases_covered_10x=\$(awk '\$3 >= 10 {n++} END{print n+0}' ${prefix}.line_probe.depth.dup.tsv)
    probe_length=\$(awk 'BEGIN{n=0} !/^>/ {n+=length(\$0)} END{print n}' ${probe_fasta})

    samtools view -F 2308 ${prefix}.line_probe.dup.bam \\
        | awk '
            BEGIN {
                OFS="\\t"
            }
            {
                cigar=\$6
                pos=\$4
                aligned=0
                soft=0
                ref_span=0
                c=cigar
                while (match(c, /^[0-9]+[MIDNSHP=X]/)) {
                    len=substr(c, RSTART, RLENGTH - 1) + 0
                    op=substr(c, RSTART + RLENGTH - 1, 1)
                    if (op == "M" || op == "=" || op == "X") aligned += len
                    if (op == "S") soft += len
                    if (op == "M" || op == "D" || op == "N" || op == "=" || op == "X") ref_span += len
                    c=substr(c, RLENGTH + 1)
                }

                end=pos + ref_span - 1
                records++
                total_aligned += aligned
                total_soft += soft
                if (aligned >= 10) aligned_ge_10++
                if (aligned >= 20) aligned_ge_20++
                if (aligned >= 30) aligned_ge_30++
                if (aligned >= 50) aligned_ge_50++
                if (aligned >= 80) aligned_ge_80++
                if (aligned >= 100) aligned_ge_100++
                if (pos <= 1 && end >= 9) cover_1_9++
                if (pos <= 1 && end >= 25) cover_1_25++
                if (pos <= 1 && end >= 31) cover_1_31++
                if (pos <= 1 && end >= 40) cover_1_40++
                if (pos >= 10) start_ge_10++
            }
            END {
                mean_aligned = records > 0 ? total_aligned / records : 0
                mean_soft = records > 0 ? total_soft / records : 0
                print total_aligned+0, total_soft+0, mean_aligned, mean_soft, aligned_ge_10+0, aligned_ge_20+0, aligned_ge_30+0, aligned_ge_50+0, aligned_ge_80+0, aligned_ge_100+0, cover_1_9+0, cover_1_25+0, cover_1_31+0, cover_1_40+0, start_ge_10+0
            }
        ' > probe_alignment_metrics.tsv

    read total_aligned_bases total_softclipped_bases mean_aligned_bases_per_record mean_softclipped_bases_per_record records_aligned_ge_10bp records_aligned_ge_20bp records_aligned_ge_30bp records_aligned_ge_50bp records_aligned_ge_80bp records_aligned_ge_100bp records_covering_probe_1_9 records_covering_probe_1_25 records_covering_probe_1_31 records_covering_probe_1_40 records_starting_at_or_after_10 < probe_alignment_metrics.tsv

    {
        echo -e "sample\\tprobe_fasta\\ttotal_records\\tmapped_records\\tmapped_primary_records\\tproper_pairs\\tmean_depth\\tprobe_length\\tbases_covered_1x\\tbases_covered_10x\\ttotal_aligned_bases\\ttotal_softclipped_bases\\tmean_aligned_bases_per_record\\tmean_softclipped_bases_per_record\\trecords_aligned_ge_10bp\\trecords_aligned_ge_20bp\\trecords_aligned_ge_30bp\\trecords_aligned_ge_50bp\\trecords_aligned_ge_80bp\\trecords_aligned_ge_100bp\\trecords_covering_probe_1_9\\trecords_covering_probe_1_25\\trecords_covering_probe_1_31\\trecords_covering_probe_1_40\\trecords_starting_at_or_after_10"
        echo -e "${prefix}\\t${fasta_name}\\t\${total_records}\\t\${mapped_records}\\t\${mapped_primary_records}\\t\${proper_pairs}\\t\${mean_depth}\\t\${probe_length}\\t\${bases_covered_1x}\\t\${bases_covered_10x}\\t\${total_aligned_bases}\\t\${total_softclipped_bases}\\t\${mean_aligned_bases_per_record}\\t\${mean_softclipped_bases_per_record}\\t\${records_aligned_ge_10bp}\\t\${records_aligned_ge_20bp}\\t\${records_aligned_ge_30bp}\\t\${records_aligned_ge_50bp}\\t\${records_aligned_ge_80bp}\\t\${records_aligned_ge_100bp}\\t\${records_covering_probe_1_9}\\t\${records_covering_probe_1_25}\\t\${records_covering_probe_1_31}\\t\${records_covering_probe_1_40}\\t\${records_starting_at_or_after_10}"
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
    touch ${prefix}.line_probe.dup.bam
    touch ${prefix}.line_probe.dup.bam.bai
    touch ${prefix}.line_probe.depth.dup.tsv
    touch ${prefix}.line_probe.summary.dup.tsv

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
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-B 3 -K 100000000 -Y -M'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fasta_name = probe_fasta.getName()

    """
    BWA_INDEX_PREFIX=\$(find -L ./ -name "*.amb" | sed 's/.amb//')

    samtools collate -u -O ${bam} ${prefix}.collate \\
        | samtools fastq \\
            -n \\
            -F 0x900 \\
            -1 /dev/stdout \\
            -2 /dev/stdout \\
            -0 /dev/null \\
            -s /dev/null \\
            - \\
        | bwa mem -p -t ${task.cpus} ${args} -R ${meta.read_group} \$BWA_INDEX_PREFIX - \\
        | samtools sort -@ ${task.cpus} -o ${prefix}.line_probe.sim.bam

    samtools index ${prefix}.line_probe.sim.bam

    samtools depth -aa ${prefix}.line_probe.sim.bam > ${prefix}.line_probe.depth.sim.tsv

    total_records=\$(samtools view -c ${prefix}.line_probe.sim.bam)
    mapped_records=\$(samtools view -c -F 4 ${prefix}.line_probe.sim.bam)
    mapped_primary_records=\$(samtools view -c -F 2308 ${prefix}.line_probe.sim.bam)
    proper_pairs=\$(samtools view -c -f 2 ${prefix}.line_probe.sim.bam)
    mean_depth=\$(awk '{sum+=\$3; n++} END{if (n > 0) print sum/n; else print 0}' ${prefix}.line_probe.depth.sim.tsv)
    bases_covered_1x=\$(awk '\$3 >= 1 {n++} END{print n+0}' ${prefix}.line_probe.depth.sim.tsv)
    bases_covered_10x=\$(awk '\$3 >= 10 {n++} END{print n+0}' ${prefix}.line_probe.depth.sim.tsv)
    probe_length=\$(awk 'BEGIN{n=0} !/^>/ {n+=length(\$0)} END{print n}' ${probe_fasta})

    samtools view -F 2308 ${prefix}.line_probe.sim.bam \\
        | awk '
            BEGIN {
                OFS="\\t"
            }
            {
                cigar=\$6
                pos=\$4
                aligned=0
                soft=0
                ref_span=0
                c=cigar
                while (match(c, /^[0-9]+[MIDNSHP=X]/)) {
                    len=substr(c, RSTART, RLENGTH - 1) + 0
                    op=substr(c, RSTART + RLENGTH - 1, 1)
                    if (op == "M" || op == "=" || op == "X") aligned += len
                    if (op == "S") soft += len
                    if (op == "M" || op == "D" || op == "N" || op == "=" || op == "X") ref_span += len
                    c=substr(c, RLENGTH + 1)
                }

                end=pos + ref_span - 1
                records++
                total_aligned += aligned
                total_soft += soft
                if (aligned >= 10) aligned_ge_10++
                if (aligned >= 20) aligned_ge_20++
                if (aligned >= 30) aligned_ge_30++
                if (aligned >= 50) aligned_ge_50++
                if (aligned >= 80) aligned_ge_80++
                if (aligned >= 100) aligned_ge_100++
                if (pos <= 1 && end >= 9) cover_1_9++
                if (pos <= 1 && end >= 25) cover_1_25++
                if (pos <= 1 && end >= 31) cover_1_31++
                if (pos <= 1 && end >= 40) cover_1_40++
                if (pos >= 10) start_ge_10++
            }
            END {
                mean_aligned = records > 0 ? total_aligned / records : 0
                mean_soft = records > 0 ? total_soft / records : 0
                print total_aligned+0, total_soft+0, mean_aligned, mean_soft, aligned_ge_10+0, aligned_ge_20+0, aligned_ge_30+0, aligned_ge_50+0, aligned_ge_80+0, aligned_ge_100+0, cover_1_9+0, cover_1_25+0, cover_1_31+0, cover_1_40+0, start_ge_10+0
            }
        ' > probe_alignment_metrics.tsv

    read total_aligned_bases total_softclipped_bases mean_aligned_bases_per_record mean_softclipped_bases_per_record records_aligned_ge_10bp records_aligned_ge_20bp records_aligned_ge_30bp records_aligned_ge_50bp records_aligned_ge_80bp records_aligned_ge_100bp records_covering_probe_1_9 records_covering_probe_1_25 records_covering_probe_1_31 records_covering_probe_1_40 records_starting_at_or_after_10 < probe_alignment_metrics.tsv

    {
        echo -e "sample\\tprobe_fasta\\ttotal_records\\tmapped_records\\tmapped_primary_records\\tproper_pairs\\tmean_depth\\tprobe_length\\tbases_covered_1x\\tbases_covered_10x\\ttotal_aligned_bases\\ttotal_softclipped_bases\\tmean_aligned_bases_per_record\\tmean_softclipped_bases_per_record\\trecords_aligned_ge_10bp\\trecords_aligned_ge_20bp\\trecords_aligned_ge_30bp\\trecords_aligned_ge_50bp\\trecords_aligned_ge_80bp\\trecords_aligned_ge_100bp\\trecords_covering_probe_1_9\\trecords_covering_probe_1_25\\trecords_covering_probe_1_31\\trecords_covering_probe_1_40\\trecords_starting_at_or_after_10"
        echo -e "${prefix}\\t${fasta_name}\\t\${total_records}\\t\${mapped_records}\\t\${mapped_primary_records}\\t\${proper_pairs}\\t\${mean_depth}\\t\${probe_length}\\t\${bases_covered_1x}\\t\${bases_covered_10x}\\t\${total_aligned_bases}\\t\${total_softclipped_bases}\\t\${mean_aligned_bases_per_record}\\t\${mean_softclipped_bases_per_record}\\t\${records_aligned_ge_10bp}\\t\${records_aligned_ge_20bp}\\t\${records_aligned_ge_30bp}\\t\${records_aligned_ge_50bp}\\t\${records_aligned_ge_80bp}\\t\${records_aligned_ge_100bp}\\t\${records_covering_probe_1_9}\\t\${records_covering_probe_1_25}\\t\${records_covering_probe_1_31}\\t\${records_covering_probe_1_40}\\t\${records_starting_at_or_after_10}"
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
    touch ${prefix}.line_probe.sim.bam
    touch ${prefix}.line_probe.sim.bam.bai
    touch ${prefix}.line_probe.depth.sim.tsv
    touch ${prefix}.line_probe.summary.sim.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | tail -n 1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
