process COLLECTHSMETRICS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.3.0--hdfd78af_0' :
        'quay.io/biocontainers/picard:3.3.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(dict)
    path(bait_intervals, stageAs: "baits/*")
    path(target_intervals, stageAs: 'targets/*')
    val(library)

    output:
    tuple val(meta), path("*.HSmetricss_MQ00.txt")         , emit: hsmetrics00
    tuple val(meta), path("*.HSmetricss_MQ20.txt")         , emit: hsmetrics20
    tuple val(meta), path("*.target_${library}s_MQ00.covg"), emit: coverage00
    tuple val(meta), path("*.target_${library}s_MQ20.covg"), emit: coverage20
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--REFERENCE_SEQUENCE ${fasta}" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard CollectHsMetrics] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    def bait_interval_list = bait_intervals
    def bait_intervallist_cmd = ""
    if (bait_intervals =~ /.(bed|bed.gz)$/){
        bait_interval_list = bait_intervals.toString().replaceAll(/.(bed|bed.gz)$/, ".interval_list")
        bait_intervallist_cmd = "picard -Xmx${avail_mem}M  BedToIntervalList --INPUT ${bait_intervals} --OUTPUT ${bait_interval_list} --SEQUENCE_DICTIONARY ${dict} --TMP_DIR ."
    }

    def target_interval_list = target_intervals
    def target_intervallist_cmd = ""
    if (target_intervals =~ /.(bed|bed.gz)$/){
        target_interval_list = target_intervals.toString().replaceAll(/.(bed|bed.gz)$/, ".interval_list")
        target_intervallist_cmd = "picard -Xmx${avail_mem}M  BedToIntervalList --INPUT ${target_intervals} --OUTPUT ${target_interval_list} --SEQUENCE_DICTIONARY ${dict} --TMP_DIR ."
    }


    """

    $bait_intervallist_cmd
    $target_intervallist_cmd

    picard \\
        -Xmx${avail_mem}M \\
        CollectHsMetrics \\
        $args \\
        $reference \\
        --BAIT_INTERVALS $bait_interval_list \\
        --TARGET_INTERVALS $target_interval_list \\
        --INPUT $bam \\
        --OUTPUT ${prefix}.HSmetrics_MQ20.txt \\
        --PER_TARGET_COVERAGE ${prefix}.target_${library}s_MQ20.covg

    picard \\
        -Xmx${avail_mem}M \\
        CollectHsMetrics \\
        $args \\
        $reference \\
        --BAIT_INTERVALS $bait_interval_list \\
        --TARGET_INTERVALS $target_interval_list \\
        --INPUT $bam \\
        --MINIMUM_MAPPING_QUALITY 0 \\
        --OUTPUT ${prefix}.HSmetrics_MQ00.txt \\
        --PER_TARGET_COVERAGE ${prefix}.target_${library}s_MQ00.covg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard CollectHsMetrics --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.HSmetrics_MQ00.txt
    touch ${prefix}.HSmetrics_MQ20.txt
    touch ${prefix}.target_${library}_MQ00.covg
    touch ${prefix}.target_${library}_MQ20.covg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard CollectHsMetrics --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}

process COLLECTHSMETRICS_CON {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.3.0--hdfd78af_0' :
        'quay.io/biocontainers/picard:3.3.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(dict)
    path(bait_intervals, stageAs: "baits/*")
    path(target_intervals, stageAs: 'targets/*')
    val(library)

    output:
    tuple val(meta), path("*.con.HSmetricss_MQ00.txt")         , emit: hsmetrics00
    tuple val(meta), path("*.con.HSmetricss_MQ20.txt")         , emit: hsmetrics20
    tuple val(meta), path("*.con.target_${library}s_MQ00.covg"), emit: coverage00
    tuple val(meta), path("*.con.target_${library}s_MQ20.covg"), emit: coverage20
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--REFERENCE_SEQUENCE ${fasta}" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard CollectHsMetrics] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    def bait_interval_list = bait_intervals
    def bait_intervallist_cmd = ""
    if (bait_intervals =~ /.(bed|bed.gz)$/){
        bait_interval_list = bait_intervals.toString().replaceAll(/.(bed|bed.gz)$/, ".interval_list")
        bait_intervallist_cmd = "picard -Xmx${avail_mem}M  BedToIntervalList --INPUT ${bait_intervals} --OUTPUT ${bait_interval_list} --SEQUENCE_DICTIONARY ${dict} --TMP_DIR ."
    }

    def target_interval_list = target_intervals
    def target_intervallist_cmd = ""
    if (target_intervals =~ /.(bed|bed.gz)$/){
        target_interval_list = target_intervals.toString().replaceAll(/.(bed|bed.gz)$/, ".interval_list")
        target_intervallist_cmd = "picard -Xmx${avail_mem}M  BedToIntervalList --INPUT ${target_intervals} --OUTPUT ${target_interval_list} --SEQUENCE_DICTIONARY ${dict} --TMP_DIR ."
    }


    """

    $bait_intervallist_cmd
    $target_intervallist_cmd

    picard \\
        -Xmx${avail_mem}M \\
        CollectHsMetrics \\
        $args \\
        $reference \\
        --BAIT_INTERVALS $bait_interval_list \\
        --TARGET_INTERVALS $target_interval_list \\
        --INPUT $bam \\
        --OUTPUT ${prefix}.con.HSmetrics_MQ20.txt \\
        --PER_TARGET_COVERAGE ${prefix}.con.target_${library}s_MQ20.covg

    picard \\
        -Xmx${avail_mem}M \\
        CollectHsMetrics \\
        $args \\
        $reference \\
        --BAIT_INTERVALS $bait_interval_list \\
        --TARGET_INTERVALS $target_interval_list \\
        --INPUT $bam \\
        --MINIMUM_MAPPING_QUALITY 0 \\
        --OUTPUT ${prefix}.con.HSmetrics_MQ00.txt \\
        --PER_TARGET_COVERAGE ${prefix}.con.target_${library}s_MQ00.covg


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard CollectHsMetrics --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.con.HSmetrics_MQ00.txt
    touch ${prefix}.con.HSmetrics_MQ20.txt
    touch ${prefix}.con.target_${library}_MQ00.covg
    touch ${prefix}.con.target_${library}_MQ20.covg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard CollectHsMetrics --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}

process COLLECTHSMETRICS_RAW {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.3.0--hdfd78af_0' :
        'quay.io/biocontainers/picard:3.3.0--hdfd78af_0' }"

    input:
    tuple val(meta),  path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(dict)
    path(bait_intervals, stageAs: "baits/*")
    path(target_intervals, stageAs: 'targets/*')
    val(library)

    output:
    tuple val(meta), path("*.raw.HSmetricss_MQ00.txt")         , emit: hsmetrics00
    tuple val(meta), path("*.raw.HSmetricss_MQ20.txt")         , emit: hsmetrics20
    tuple val(meta), path("*.raw.target_${library}s_MQ00.covg"), emit: coverage00
    tuple val(meta), path("*.raw.target_${library}s_MQ20.covg"), emit: coverage20
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--REFERENCE_SEQUENCE ${fasta}" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard CollectHsMetrics] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    def bait_interval_list = bait_intervals
    def bait_intervallist_cmd = ""
    if (bait_intervals =~ /.(bed|bed.gz)$/){
        bait_interval_list = bait_intervals.toString().replaceAll(/.(bed|bed.gz)$/, ".interval_list")
        bait_intervallist_cmd = "picard -Xmx${avail_mem}M  BedToIntervalList --INPUT ${bait_intervals} --OUTPUT ${bait_interval_list} --SEQUENCE_DICTIONARY ${dict} --TMP_DIR ."
    }

    def target_interval_list = target_intervals
    def target_intervallist_cmd = ""
    if (target_intervals =~ /.(bed|bed.gz)$/){
        target_interval_list = target_intervals.toString().replaceAll(/.(bed|bed.gz)$/, ".interval_list")
        target_intervallist_cmd = "picard -Xmx${avail_mem}M  BedToIntervalList --INPUT ${target_intervals} --OUTPUT ${target_interval_list} --SEQUENCE_DICTIONARY ${dict} --TMP_DIR ."
    }


    """

    $bait_intervallist_cmd
    $target_intervallist_cmd

    picard \\
        -Xmx${avail_mem}M \\
        CollectHsMetrics \\
        $args \\
        $reference \\
        --BAIT_INTERVALS $bait_interval_list \\
        --TARGET_INTERVALS $target_interval_list \\
        --INPUT $bam \\
        --OUTPUT ${prefix}.raw.HSmetrics_MQ20.txt \\
        --PER_TARGET_COVERAGE ${prefix}.raw.target_${library}s_MQ20.covg

    picard \\
        -Xmx${avail_mem}M \\
        CollectHsMetrics \\
        $args \\
        $reference \\
        --BAIT_INTERVALS $bait_interval_list \\
        --TARGET_INTERVALS $target_interval_list \\
        --INPUT $bam \\
        --MINIMUM_MAPPING_QUALITY 0 \\
        --OUTPUT ${prefix}.raw.HSmetrics_MQ00.txt \\
        --PER_TARGET_COVERAGE ${prefix}.raw.target_${library}s_MQ00.covg


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard CollectHsMetrics --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.raw.HSmetrics_MQ00.txt
    touch ${prefix}.raw.HSmetrics_MQ20.txt
    touch ${prefix}.raw.target_${library}_MQ00.covg
    touch ${prefix}.raw.target_${library}_MQ20.covg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard CollectHsMetrics --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}

process COLLECTHSMETRICS_DUP {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.3.0--hdfd78af_0' :
        'quay.io/biocontainers/picard:3.3.0--hdfd78af_0' }"

    input:
    tuple val(meta),  path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(dict)
    path(bait_intervals, stageAs: "baits/*")
    path(target_intervals, stageAs: 'targets/*')
    val(library)

    output:
    tuple val(meta), path("*.dup.HSmetricss_MQ00.txt")         , emit: hsmetrics00
    tuple val(meta), path("*.dup.HSmetricss_MQ20.txt")         , emit: hsmetrics20
    tuple val(meta), path("*.dup.target_${library}s_MQ00.covg"), emit: coverage00
    tuple val(meta), path("*.dup.target_${library}s_MQ20.covg"), emit: coverage20
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--REFERENCE_SEQUENCE ${fasta}" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard CollectHsMetrics] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    def bait_interval_list = bait_intervals
    def bait_intervallist_cmd = ""
    if (bait_intervals =~ /.(bed|bed.gz)$/){
        bait_interval_list = bait_intervals.toString().replaceAll(/.(bed|bed.gz)$/, ".interval_list")
        bait_intervallist_cmd = "picard -Xmx${avail_mem}M  BedToIntervalList --INPUT ${bait_intervals} --OUTPUT ${bait_interval_list} --SEQUENCE_DICTIONARY ${dict} --TMP_DIR ."
    }

    def target_interval_list = target_intervals
    def target_intervallist_cmd = ""
    if (target_intervals =~ /.(bed|bed.gz)$/){
        target_interval_list = target_intervals.toString().replaceAll(/.(bed|bed.gz)$/, ".interval_list")
        target_intervallist_cmd = "picard -Xmx${avail_mem}M  BedToIntervalList --INPUT ${target_intervals} --OUTPUT ${target_interval_list} --SEQUENCE_DICTIONARY ${dict} --TMP_DIR ."
    }


    """

    $bait_intervallist_cmd
    $target_intervallist_cmd

    picard \\
        -Xmx${avail_mem}M \\
        CollectHsMetrics \\
        $args \\
        $reference \\
        --BAIT_INTERVALS $bait_interval_list \\
        --TARGET_INTERVALS $target_interval_list \\
        --INPUT $bam \\
        --OUTPUT ${prefix}.dup.HSmetrics_MQ20.txt \\
        --PER_TARGET_COVERAGE ${prefix}.dup.target_${library}s_MQ20.covg

    picard \\
        -Xmx${avail_mem}M \\
        CollectHsMetrics \\
        $args \\
        $reference \\
        --BAIT_INTERVALS $bait_interval_list \\
        --TARGET_INTERVALS $target_interval_list \\
        --INPUT $bam \\
        --MINIMUM_MAPPING_QUALITY 0 \\
        --OUTPUT ${prefix}.dup.HSmetrics_MQ00.txt \\
        --PER_TARGET_COVERAGE ${prefix}.dup.target_${library}s_MQ00.covg


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard CollectHsMetrics --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.dup.HSmetrics_MQ00.txt
    touch ${prefix}.dup.HSmetrics_MQ20.txt
    touch ${prefix}.dup.target_${library}_MQ00.covg
    touch ${prefix}.dup.target_${library}_MQ20.covg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard CollectHsMetrics --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}

process COLLECTHSMETRICS_SIM {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.3.0--hdfd78af_0' :
        'quay.io/biocontainers/picard:3.3.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(dict)
    path(bait_intervals, stageAs: "baits/*")
    path(target_intervals, stageAs: 'targets/*')
    val(library)

    output:
    tuple val(meta), path("*.sim.HSmetricss_MQ00.txt")         , emit: hsmetrics00
    tuple val(meta), path("*.sim.HSmetricss_MQ20.txt")         , emit: hsmetrics20
    tuple val(meta), path("*.sim.target_${library}s_MQ00.covg"), emit: coverage00
    tuple val(meta), path("*.sim.target_${library}s_MQ20.covg"), emit: coverage20
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--REFERENCE_SEQUENCE ${fasta}" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard CollectHsMetrics] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    def bait_interval_list = bait_intervals
    def bait_intervallist_cmd = ""
    if (bait_intervals =~ /.(bed|bed.gz)$/){
        bait_interval_list = bait_intervals.toString().replaceAll(/.(bed|bed.gz)$/, ".interval_list")
        bait_intervallist_cmd = "picard -Xmx${avail_mem}M  BedToIntervalList --INPUT ${bait_intervals} --OUTPUT ${bait_interval_list} --SEQUENCE_DICTIONARY ${dict} --TMP_DIR ."
    }

    def target_interval_list = target_intervals
    def target_intervallist_cmd = ""
    if (target_intervals =~ /.(bed|bed.gz)$/){
        target_interval_list = target_intervals.toString().replaceAll(/.(bed|bed.gz)$/, ".interval_list")
        target_intervallist_cmd = "picard -Xmx${avail_mem}M  BedToIntervalList --INPUT ${target_intervals} --OUTPUT ${target_interval_list} --SEQUENCE_DICTIONARY ${dict} --TMP_DIR ."
    }


    """
    $bait_intervallist_cmd
    $target_intervallist_cmd

    picard \\
        -Xmx${avail_mem}M \\
        CollectHsMetrics \\
        $args \\
        $reference \\
        --BAIT_INTERVALS $bait_interval_list \\
        --TARGET_INTERVALS $target_interval_list \\
        --INPUT $bam \\
        --OUTPUT ${prefix}.sim.HSmetrics_MQ20.txt \\
        --PER_TARGET_COVERAGE ${prefix}.sim.target_${library}s_MQ20.covg

    picard \\
        -Xmx${avail_mem}M \\
        CollectHsMetrics \\
        $args \\
        $reference \\
        --BAIT_INTERVALS $bait_interval_list \\
        --TARGET_INTERVALS $target_interval_list \\
        --INPUT $bam \\
        --MINIMUM_MAPPING_QUALITY 0 \\
        --OUTPUT ${prefix}.sim.HSmetrics_MQ00.txt \\
        --PER_TARGET_COVERAGE ${prefix}.sim.target_${library}s_MQ00.covg


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard CollectHsMetrics --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sim.HSmetrics_MQ00.txt
    touch ${prefix}.sim.HSmetrics_MQ20.txt
    touch ${prefix}.sim.target_${library}_MQ00.covg
    touch ${prefix}.sim.target_${library}_MQ20.covg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard CollectHsMetrics --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}

process COLLECTHSMETRICS_DR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.3.0--hdfd78af_0' :
        'quay.io/biocontainers/picard:3.3.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(dict)
    path(bait_intervals, stageAs: "baits/*")
    path(target_intervals, stageAs: 'targets/*')
    val(library)

    output:
    tuple val(meta), path("*.dr.HSmetricss_MQ00.txt")         , emit: hsmetrics00
    tuple val(meta), path("*.dr.HSmetricss_MQ20.txt")         , emit: hsmetrics20
    tuple val(meta), path("*.dr.target_${library}s_MQ00.covg"), emit: coverage00
    tuple val(meta), path("*.dr.target_${library}s_MQ20.covg"), emit: coverage20
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--REFERENCE_SEQUENCE ${fasta}" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard CollectHsMetrics] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    def bait_interval_list = bait_intervals
    def bait_intervallist_cmd = ""
    if (bait_intervals =~ /.(bed|bed.gz)$/){
        bait_interval_list = bait_intervals.toString().replaceAll(/.(bed|bed.gz)$/, ".interval_list")
        bait_intervallist_cmd = "picard -Xmx${avail_mem}M  BedToIntervalList --INPUT ${bait_intervals} --OUTPUT ${bait_interval_list} --SEQUENCE_DICTIONARY ${dict} --TMP_DIR ."
    }

    def target_interval_list = target_intervals
    def target_intervallist_cmd = ""
    if (target_intervals =~ /.(bed|bed.gz)$/){
        target_interval_list = target_intervals.toString().replaceAll(/.(bed|bed.gz)$/, ".interval_list")
        target_intervallist_cmd = "picard -Xmx${avail_mem}M  BedToIntervalList --INPUT ${target_intervals} --OUTPUT ${target_interval_list} --SEQUENCE_DICTIONARY ${dict} --TMP_DIR ."
    }


    """
    $bait_intervallist_cmd
    $target_intervallist_cmd

    picard \\
        -Xmx${avail_mem}M \\
        CollectHsMetrics \\
        $args \\
        $reference \\
        --BAIT_INTERVALS $bait_interval_list \\
        --TARGET_INTERVALS $target_interval_list \\
        --INPUT $bam \\
        --OUTPUT ${prefix}.dr.HSmetrics_MQ20.txt \\
        --PER_TARGET_COVERAGE ${prefix}.dr.target_${library}s_MQ20.covg

    picard \\
        -Xmx${avail_mem}M \\
        CollectHsMetrics \\
        $args \\
        $reference \\
        --BAIT_INTERVALS $bait_interval_list \\
        --TARGET_INTERVALS $target_interval_list \\
        --INPUT $bam \\
        --MINIMUM_MAPPING_QUALITY 0 \\
        --OUTPUT ${prefix}.dr.HSmetrics_MQ00.txt \\
        --PER_TARGET_COVERAGE ${prefix}.dr.target_${library}s_MQ00.covg


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard CollectHsMetrics --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.dr.HSmetrics_MQ00.txt
    touch ${prefix}.dr.HSmetrics_MQ20.txt
    touch ${prefix}.dr.target_${library}_MQ00.covg
    touch ${prefix}.dr.target_${library}_MQ20.covg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard CollectHsMetrics --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
