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
    def degree                     = task.ext.args?["general"]?["degree"]                       ? "degree = ${task.ext.args["general"]["degree"]}"                                              : ""
    def intercept                  = task.ext.args?["general"]?["intercept"]                    ? "intercept = ${task.ext.args["general"]["intercept"]}"                                        : ""
    def mincnalength               = task.ext.args?["general"]?["mincnalength"]                 ? "minCNAlength = ${task.ext.args["general"]["mincnalength"]}"                                  : ""
    def minmappabilityperwindow    = task.ext.args?["general"]?["minmappabilityperwindow"]      ? "minMappabilityPerWindow = ${task.ext.args["general"]["minmappabilityperwindow"]}"            : ""
    def minexpectedgc              = task.ext.args?["general"]?["minexpectedgc"]                ? "minExpectedGC = ${task.ext.args["general"]["minexpectedgc"]}"                                : ""
    def maxexpectedgc              = task.ext.args?["general"]?["maxexpectedgc"]                ? "maxExpectedGC = ${task.ext.args["general"]["maxexpectedgc"]}"                                : ""
    def ploidy                     = task.ext.args?["general"]?["ploidy"]                       ? "ploidy = ${task.ext.args["general"]["ploidy"]}"                                              : ""
    def sex                        = task.ext.args?["general"]?["sex"]                          ? "sex = ${task.ext.args["general"]["sex"]}"                                                    : ""
    def step                       = task.ext.args?["general"]?["step"]                         ? "step = ${task.ext.args["general"]["step"]}"                                                  : ""
    def telocentromeric            = task.ext.args?["general"]?["telocentromeric"]              ? "telocentromeric = ${task.ext.args["general"]["telocentromeric"]} "                           : ""
    def uniquematch                = task.ext.args?["general"]?["uniquematch"]                  ? "uniqueMatch = ${task.ext.args["general"]["uniquematch"]}"                                    : ""
    def window                     = task.ext.args?["general"]?["window"]                       ? "window = ${task.ext.args["general"]["window"]}"                                              : ""
    def matefile_tumor             = mpileup_tumor                                              ? "mateFile = \${PWD}/${mpileup_tumor}"                                                         : ""
    def inputformat_tumor          = task.ext.args?["sample"]?["inputformat"]                   ? "inputFormat = ${task.ext.args["sample"]["inputformat"]}"                                     : ""
    def mateorientation_tumor      = task.ext.args?["sample"]?["mateorientation"]               ? "mateOrientation = ${task.ext.args["sample"]["mateorientation"]}"                             : ""
    def fastafile                  = fasta                                                      ? "fastaFile = \${PWD}/${fasta}"                                                                : ""
    def minimalcoverageperposition = task.ext.args?["BAF"]?["minimalcoverageperposition"]       ? "minimalCoveragePerPosition = ${task.ext.args["BAF"]["minimalcoverageperposition"]}"          : ""
    def minimalqualityperposition  = task.ext.args?["BAF"]?["minimalqualityperposition"]        ? "minimalQualityPerPosition = ${task.ext.args["BAF"]["minimalqualityperposition"]}"            : ""
    def shiftinquality             = task.ext.args?["BAF"]?["shiftinquality"]                   ? "shiftInQuality = ${task.ext.args["BAF"]["shiftinquality"]}"                                  : ""
    def snpfile                    = known_snps                                                 ? "SNPfile = \$PWD/${known_snps}"                                                               : ""
    def target_bed                 = target_bed                                                 ? "captureRegions = ${target_bed}"                                                              : ""
    def VERSION                    = "11.6b"
    """
    cat <<-EOF > config.txt
[general]
noisyData = TRUE
BedGraphOutput = TRUE
breakPointThreshold = 0.8
breakPointType = 2
chrFiles = \${PWD}/chrs/
chrLenFile = \${PWD}/${fai}
coefficientOfVariation = ${cf_coeff}
contamination = ${cf_contamination}
contaminationAdjustment = ${cf_contamination_adjustment}
${degree}
forceGCcontentNormalization = 1
gemMappabilityFile = \${PWD}/${mappability}
${intercept}
${mincnalength}
${minmappabilityperwindow}
${minexpectedgc}
${maxexpectedgc}
minimalSubclonePresence = 20
maxThreads = ${task.cpus}
outputDir = \${PWD}/
ploidy = ${cf_ploidy}
printNA = TRUE
readCountThreshold = 10
${sex}
${step}
${telocentromeric}
${uniquematch}
${window}

[control]
mateFile = ${mateFile}
inputFormat = pileup
mateOrientation = FR

[sample]
${matefile_tumor}
${inputformat_tumor}
${mateorientation_tumor}

[BAF]
${fastafile}
${minimalcoverageperposition}
${minimalqualityperposition}
${shiftinquality}
${snpfile}

[target]
${target_bed}
EOF

    mkdir -p chrs/

    for file in ${chr_directory}/Homo_sapiens.GRCh38.dna.chromosome.*.fa; do
        base_name=\$(basename "\$file")
        new_name="\${base_name#Homo_sapiens.GRCh38.dna.chromosome.}"
        cp "\$file" "chrs/\$new_name"
    done

    freec -conf config.txt -sample ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        controlfreec: "${VERSION}"
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '11.6b'
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
        controlfreec: "${VERSION}"
    END_VERSIONS
    """
}
