process ONCOCNV {
    tag "$meta.patient"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1f/1f04dcb708961e66aade70ebd0f9faa1a2980e2eff66bdbd8c947575c8ba4f31/data':
        'community.wave.seqera.io/library/oncocnv_r-pscbs:ffe800b4a6974d56' }"

    input:
    tuple val(meta) , path(tumor), path(tumor_index), path(normal), path(normal_index)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(bed)
    tuple val(meta4), path(fai)

    output:
    tuple val(meta), path("*.profile.png"), emit: png
    tuple val(meta), path("*.profile.txt"), emit: profile
    tuple val(meta), path("*.summary.txt"), emit: summary
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def cghseg = task.ext.args2 ?: 'cghseg'
    def mode = task.ext.args ?: '-m Ampli'
    def normal_id = normal.join(',')
    def tumor_id = tumor.join(',')
    def VERSION = '7.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    export R_USER_CACHE_DIR=\$(pwd)/R_cache

    awk 'BEGIN {OFS="\\t"} {printf "%s\\t%s\\t%s\\tID%04d\\t0\\t%s\\n", \$1, \$2, \$3, ++i, \$5}' ${bed} > targets.bed

    perl \$(which ONCOCNV_getCounts.pl) \\
        getControlStats \\
        $mode \\
        -b targets.bed \\
        -c $normal_id \\
        -o ControlStats.txt

    perl \$(which ONCOCNV_getCounts.pl) \\
        getSampleStats \\
        $mode \\
        -c ControlStats.txt \\
        -s $tumor_id \\
        -o SampleStats.txt

    cat ControlStats.txt \\
        | grep -v start \\
        | awk '{print \$1,\$2,\$3}' \\
        | sed "s/ /\t/g" > target.bed

    perl \$(which createTargetGC.pl) \\
        -bed target.bed \\
        -fi ${fasta} \\
        -od . \\
        -of TargetGC.txt

    cat \$(which processControl.R) \\
        | R \\
        --slave \\
        --args ControlStats.txt ControlStatsProcessed.txt TargetGC.txt

    cat \$(which processSamples.R) \\
        | R \\
        --slave \\
        --args SampleStats.txt ControlStatsProcessed.txt Output.log $cghseg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oncocnv: $VERSION
        perl: \$(perl --version | grep 'This is perl' | sed -E 's/.*v([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/')
        r: \$(R --version | grep "R version" | sed -E 's/R version ([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/')
    END_VERSIONS
    """

    stub:
    def VERSION = '7.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.profile.png
    touch ${prefix}.profile.txt
    touch ${prefix}.summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oncocnv: $VERSION
        perl: \$(perl --version | grep 'This is perl' | sed -E 's/.*v([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/')
        r: \$(R --version | grep "R version" | sed -E 's/R version ([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/')
    END_VERSIONS
    """
}
