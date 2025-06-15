process ONCOCNV {
    tag "$meta.patient"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/oncocnv:7.0.0' :
        'blancojmskcc/oncocnv:7.0.0' }"

    input:
    tuple val(meta) , path(tumor), path(tumor_index), path(normal), path(normal_index)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(bed)
    tuple val(meta4), path(fai)

    output:
    tuple val(meta), path("*.profile.png")  ,emit: png
    tuple val(meta), path("*.profile.txt")  ,emit: profile
    tuple val(meta), path("*.summary.txt")  ,emit: summary
    path "versions.yml"                     ,emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def cghseg = task.ext.args2 ?: 'cghseg'
    def mode = task.ext.args ?: '-m Ampli'
    def normal = normal.join(',')
    def tumor = tumor.join(',')
    def TOOLDIR = '/usr/local/bin/'
    def VERSION = '7.0' 
    """
    export R_USER_CACHE_DIR=\$(pwd)

    awk 'BEGIN {OFS="\\t"} {printf "%s\\t%s\\t%s\\tID%04d\\t0\\t%s\\n", \$1, \$2, \$3, ++i, \$5}' ${bed} > targets.bed

    perl ${TOOLDIR}ONCOCNV_getCounts.pl \\
        getControlStats \\
        ${mode} \\
        -b targets.bed \\
        -c ${normal} \\
        -o ControlStats.txt

    perl ${TOOLDIR}ONCOCNV_getCounts.pl \\
        getSampleStats \\
        ${mode} \\
        -c ControlStats.txt \\
        -s ${tumor} \\
        -o SampleStats.txt

    cat ControlStats.txt \\
        | grep -v start \\
        | awk '{print \$1,\$2,\$3}' \\
        | sed "s/ /\t/g" > target.bed

    perl ${TOOLDIR}createTargetGC.pl \\
        -bed target.bed \\
        -fi ${fasta} \\
        -od . \\
        -of TargetGC.txt

    cat ${TOOLDIR}processControl.R \\
        | R \\
        --slave \\
        --args ControlStats.txt ControlStatsProcessed.txt TargetGC.txt

    cat ${TOOLDIR}processSamples.R \\
        | R \\
        --slave \\
        --args SampleStats.txt ControlStatsProcessed.txt Output.log ${cghseg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oncocnv: ${VERSION}
        perl: \$(perl --version | grep 'This is perl' | sed 's/.*(v//g' | sed 's/)//g')
        r: \$(R --version | grep "R version" | sed 's/R version //g')
    END_VERSIONS
    """
}
