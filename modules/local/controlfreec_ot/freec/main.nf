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
    def args_cfg                   = task.ext.args ?: [:]
    def general_cfg                = args_cfg.general ?: [:]
    def sample_cfg                 = args_cfg.sample ?: [:]
    def baf_cfg                    = args_cfg.BAF ?: [:]
    def degree                     = general_cfg.degree                    ? "degree = ${general_cfg.degree}"                                          : ""
    def intercept                  = general_cfg.intercept                 ? "intercept = ${general_cfg.intercept}"                                    : ""
    def mincnalength               = general_cfg.mincnalength              ? "minCNAlength = ${general_cfg.mincnalength}"                              : ""
    def minmappabilityperwindow    = general_cfg.minmappabilityperwindow   ? "minMappabilityPerWindow = ${general_cfg.minmappabilityperwindow}"        : ""
    def minexpectedgc              = general_cfg.minexpectedgc             ? "minExpectedGC = ${general_cfg.minexpectedgc}"                            : ""
    def maxexpectedgc              = general_cfg.maxexpectedgc             ? "maxExpectedGC = ${general_cfg.maxexpectedgc}"                            : ""
    def ploidy                     = general_cfg.ploidy                    ? "ploidy = ${general_cfg.ploidy}"                                          : ""
    def sex                        = general_cfg.sex                       ? "sex = ${general_cfg.sex}"                                                : ""
    def step                       = general_cfg.step                      ? "step = ${general_cfg.step}"                                              : ""
    def telocentromeric            = general_cfg.telocentromeric           ? "telocentromeric = ${general_cfg.telocentromeric} "                       : ""
    def uniquematch                = general_cfg.uniquematch               ? "uniqueMatch = ${general_cfg.uniquematch}"                                : ""
    def window                     = general_cfg.window                    ? "window = ${general_cfg.window}"                                          : ""
    def matefile_tumor             = mpileup_tumor                                              ? "mateFile = \${PWD}/${mpileup_tumor}"                                                         : ""
    def inputformat_tumor          = sample_cfg.inputformat                   ? "inputFormat = ${sample_cfg.inputformat}"                                 : ""
    def mateorientation_tumor      = sample_cfg.mateorientation               ? "mateOrientation = ${sample_cfg.mateorientation}"                         : ""
    def fastafile                  = fasta                                                      ? "fastaFile = \${PWD}/${fasta}"                                                                : ""
    def minimalcoverageperposition = baf_cfg.minimalcoverageperposition       ? "minimalCoveragePerPosition = ${baf_cfg.minimalcoverageperposition}"      : ""
    def minimalqualityperposition  = baf_cfg.minimalqualityperposition        ? "minimalQualityPerPosition = ${baf_cfg.minimalqualityperposition}"        : ""
    def shiftinquality             = baf_cfg.shiftinquality                   ? "shiftInQuality = ${baf_cfg.shiftinquality}"                              : ""
    def snpfile                    = known_snps                                                 ? "SNPfile = \${PWD}/${known_snps}"                                                             : ""
    def capture_regions            = target_bed                                                 ? "captureRegions = ${target_bed}"                                                              : ""
    def VERSION                    = '11.6'
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
readCountThreshold = 50
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
${capture_regions}
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
        controlfreec: ${VERSION}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '11.6'
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
        controlfreec: ${VERSION}
    END_VERSIONS
    """
}
