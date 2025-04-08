process CONTROLFREEC_OT_FREEC {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/control-freec:11.6b--hdbdd923_0' :
        'quay.io/biocontainers/control-freec:11.6b--hdbdd923_0' }"

    input:
    tuple val(meta), path(mpileup_tumor)
    tuple val(meta1), path(bam), path(bai)
    path fasta
    path fai
    path known_snps
    path known_snps_tbi
    path chr_directory
    path mappability
    path target_bed
    path mateFile
    val cf_coeff
    val cf_contamination
    val cf_contamination_adjustment
    val cf_ploidy

    output:
    tuple val(meta), path("*_ratio.BedGraph")   , emit: bedgraph, optional: true
    tuple val(meta), path("*_control.cpn")      , emit: control_cpn, optional: true
    tuple val(meta), path("*_sample.cpn")       , emit: sample_cpn
    tuple val(meta), path("GC_profile.*.cpn")   , emit: gcprofile_cpn, optional:true
    tuple val(meta), path("*_BAF.txt")          , emit: BAF
    tuple val(meta), path("*_CNVs")             , emit: CNV
    tuple val(meta), path("*_info.txt")         , emit: info
    tuple val(meta), path("*_ratio.txt")        , emit: ratio
    tuple val(meta), path("config.txt")         , emit: config
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    //"General" configurations
    def degree                      = task.ext.args?["general"]?["degree"]                      ? "degree = ${task.ext.args["general"]["degree"]}"                                              : ""
    def intercept                   = task.ext.args?["general"]?["intercept"]                   ? "intercept = ${task.ext.args["general"]["intercept"]}"                                        : ""
    def mincnalength                = task.ext.args?["general"]?["mincnalength"]                ? "minCNAlength = ${task.ext.args["general"]["mincnalength"]}"                                  : ""
    def minmappabilityperwindow     = task.ext.args?["general"]?["minmappabilityperwindow"]     ? "minMappabilityPerWindow = ${task.ext.args["general"]["minmappabilityperwindow"]}"            : ""
    def minexpectedgc               = task.ext.args?["general"]?["minexpectedgc"]               ? "minExpectedGC = ${task.ext.args["general"]["minexpectedgc"]}"                                : ""
    def maxexpectedgc               = task.ext.args?["general"]?["maxexpectedgc"]               ? "maxExpectedGC = ${task.ext.args["general"]["maxexpectedgc"]}"                                : ""
    def ploidy                      = task.ext.args?["general"]?["ploidy"]                      ? "ploidy = ${task.ext.args["general"]["ploidy"]}"                                              : ""
    def sex                         = task.ext.args?["general"]?["sex"]                         ? "sex = ${task.ext.args["general"]["sex"]}"                                                    : ""
    def step                        = task.ext.args?["general"]?["step"]                        ? "step = ${task.ext.args["general"]["step"]}"                                                  : ""
    def telocentromeric             = task.ext.args?["general"]?["telocentromeric"]             ? "telocentromeric = ${task.ext.args["general"]["telocentromeric"]} "                           : ""
    def uniquematch                 = task.ext.args?["general"]?["uniquematch"]                 ? "uniqueMatch = ${task.ext.args["general"]["uniquematch"]}"                                    : ""
    def window                      = task.ext.args?["general"]?["window"]                      ? "window = ${task.ext.args["general"]["window"]}"                                              : ""

    //"Sample" configuration
    def matefile_tumor             = mpileup_tumor                                              ? "mateFile = \${PWD}/${mpileup_tumor}"                                                         : ""
    def inputformat_tumor          = task.ext.args?["sample"]?["inputformat"]                   ? "inputFormat = ${task.ext.args["sample"]["inputformat"]}"                                     : ""
    def mateorientation_tumor      = task.ext.args?["sample"]?["mateorientation"]               ? "mateOrientation = ${task.ext.args["sample"]["mateorientation"]}"                             : ""

    //"BAF" configuration
    def fastafile                  = fasta                                                      ? "fastaFile = \${PWD}/${fasta}"                                                                : ""
    def minimalcoverageperposition = task.ext.args?["BAF"]?["minimalcoverageperposition"]       ? "minimalCoveragePerPosition = ${task.ext.args["BAF"]["minimalcoverageperposition"]}"          : ""
    def minimalqualityperposition  = task.ext.args?["BAF"]?["minimalqualityperposition"]        ? "minimalQualityPerPosition = ${task.ext.args["BAF"]["minimalqualityperposition"]}"            : ""
    def shiftinquality             = task.ext.args?["BAF"]?["shiftinquality"]                   ? "shiftInQuality = ${task.ext.args["BAF"]["shiftinquality"]}"                                  : ""
    def snpfile                    = known_snps                                                 ? "SNPfile = \$PWD/${known_snps}"                                                               : ""

    //"Target" configuration
    def target_bed                 = target_bed                                                 ? "captureRegions = ${target_bed}"                                                              : ""

    def VERSION = '11.6b' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch config.txt

    echo "[general]" >> config.txt
    echo "noisyData = TRUE" >> config.txt
    echo "BedGraphOutput = TRUE" >> config.txt
    echo "breakPointThreshold = 0.8"
    echo "breakPointType = 2"
    echo "chrFiles =\${PWD}/chrs/" >> config.txt
    echo "chrLenFile = \${PWD}/${fai}"  >> config.txt
    echo "coefficientOfVariation = ${cf_coeff}" >> config.txt
    echo "contamination = ${cf_contamination} >> config.txt
    echo "contaminationAdjustment = ${cf_contamination_adjustment} >> config.txt
    echo ${degree} >> config.txt
    echo "forceGCcontentNormalization = 1" >> config.txt
    echo "gemMappabilityFile = \${PWD}/${mappability}" >> config.txt
    echo ${intercept} >> config.txt
    echo ${mincnalength} >> config.txt
    echo ${minmappabilityperwindow} >> config.txt
    echo ${minexpectedgc} >> config.txt
    echo ${maxexpectedgc} >> config.txt
    echo "minimalSubclonePresence = 20" >> config.txt
    echo "maxThreads = ${task.cpus}" >> config.txt
    echo "outputDir = \${PWD}/"  >> config.txt
    echo "ploidy = ${cf_ploidy}" >> config.txt
    echo "printNA = TRUE" >> config.txt
    echo "readCountThreshold = 10" >> config.txt

    echo ${sex} >> config.txt
    echo ${step} >> config.txt
    echo ${telocentromeric} >> config.txt
    echo ${uniquematch} >> config.txt
    echo ${window} >> config.txt
    echo "[control]" >> config.txt
    echo "mateFile = ${mateFile}" >> config.txt
    echo "inputFormat = pileup" >> config.txt
    echo "mateOrientation = FR" >> config.txt
    echo "[sample]" >> config.txt
    echo ${matefile_tumor} >> config.txt
    echo "inputFormat = pileup" >> config.txt
    echo "mateOrientation = FR" >> config.txt
    echo "[BAF]" >> config.txt
    echo ${fastafile} >> config.txt
    echo ${minimalcoverageperposition} >> config.txt
    echo ${minimalqualityperposition} >> config.txt
    echo ${shiftinquality} >> config.txt
    echo ${snpfile} >> config.txt
    echo "[target]" >> config.txt
    echo ${target_bed} >> config.txt

    mkdir -p chrs/

    for file in $chr_directory/Homo_sapiens.GRCh38.dna.chromosome.*.fa; do
        base_name=\$(basename "\$file")
        new_name="\${base_name#Homo_sapiens.GRCh38.dna.chromosome.}"
        cp "\$file" "chrs/\$new_name"
    done

    freec -conf config.txt -sample $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        controlfreec: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '11.6b' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}_ratio.BedGraph
    touch ${prefix}_sample.cpn
    touch GC_profile.${prefix}.cpn
    touch ${prefix}_BAF.txt
    touch ${prefix}_CNVs
    touch ${prefix}_info.txt
    touch ${prefix}_ratio.txt
    touch config.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        controlfreec: $VERSION
    END_VERSIONS
    """
}
