/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                       DEDUPANDRECAL SUBWORKFLOW                                                    
*******************************************************************************************************************************************************************************************************

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                             IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTP                                                                                                                     } from '../../modules/nf-core/fastp/main'   
include { MOSDEPTH                                                                                                                  } from '../../modules/local/mosdepth/main'
include { FASTQC_DR                                                                                                                 } from '../../modules/local/fastqc/main'
include { BWAMEM2_MEM                                                                                                               } from '../../modules/nf-core/bwamem2/mem/main' 
include { MOSDEPTH_DR                                                                                                               } from '../../modules/local/mosdepth/main'
include { MSISENSORPRO                                                                                                              } from '../../modules/local/msisensorpro/pro/main'   
include { SAMTOOLS_SORT                                                                                                             } from '../../modules/nf-core/samtools/sort/main' 
include { SAMTOOLS_INDEX                                                                                                            } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_STATS                                                                                                            } from '../../modules/local/samtools/stats/main'
include { GATK4_APPLYBQSR                                                                                                           } from '../../modules/local/gatk4/applybqsr/main'
include { MSISENSORPRO_DR                                                                                                           } from '../../modules/local/msisensorpro/pro/main'   
include { COLLECTHSMETRICS                                                                                                          } from '../../modules/local/picard/collecthsmetrics/main'
include { SAMTOOLS_STATS_DR                                                                                                         } from '../../modules/local/samtools/stats/main'
include { COLLECTHSMETRICS_DR                                                                                                       } from '../../modules/local/picard/collecthsmetrics/main'
include { SAMTOOLS_SORT_INDEX                                                                                                       } from '../../modules/local/samtools/sort_index/main'
include { SURVIVOR_SCAN_READS                                                                                                       } from '../../modules/local/survivor/scanreads/main'
include { GATK4_MARKDUPLICATES          	                                                                                        } from '../../modules/local/gatk4/markduplicates/main'
include { SAMTOOLS_COLLATEFASTQ                                                                                                     } from '../../modules/nf-core/samtools/collatefastq/main'   
include { GATK4_BASERECALIBRATOR                                                                                                    } from '../../modules/local/gatk4/baserecalibrator/main'
include { SURVIVOR_SCAN_READS_DR                                                                                                    } from '../../modules/local/survivor/scanreads/main'
include { FGBIO_ERRORRATEBYREADPOSITION                                                                                             } from '../../modules/local/fgbio/errorratebyreadposition/main'
include { PICARD_COLLECTMULTIPLEMETRICS                                                                                             } from '../../modules/local/picard/collectmultiplemetrics/main'
include { FGBIO_ERRORRATEBYREADPOSITION_DR                                                                                          } from '../../modules/local/fgbio/errorratebyreadposition/main'
include { PICARD_COLLECTMULTIPLEMETRICS_DR                                                                                          } from '../../modules/local/picard/collectmultiplemetrics/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                          IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { softwareVersionsToYAML                                                                                                    } from '../nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                            RUN SUBWORKFLOW 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DEDUPANDRECAL {

    take:
    ch_fai
    ch_bwa2
    ch_dict
    ch_fasta
    ch_msi_f
    ch_fastqs

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: BWA-MEM2 mapping
    //
    save_merged = 'false'
    save_trimmed_fail = 'false'
    discard_trimmed_pass = 'true'
    FASTP(ch_fastqs, discard_trimmed_pass, save_trimmed_fail, save_merged)
    ch_fastqs_fastp = FASTP.out.reads
    ch_versions = ch_versions.mix(FASTP.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json)

    //
    // MODULE: BWA-MEM2 mapping
    //
    sort_bam = 'sort'
    BWAMEM2_MEM(ch_fastqs_fastp, ch_bwa2, ch_fasta, ch_fai, sort_bam)
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions)
    ch_bam = BWAMEM2_MEM.out.bam

    //
    // MODULE: Run SamToools Sort & Index
    //
    SAMTOOLS_SORT_INDEX(ch_bam, ch_fasta, params.fai)
    ch_versions = ch_versions.mix(SAMTOOLS_SORT_INDEX.out.versions.first())
    ch_bam_raw = SAMTOOLS_SORT_INDEX.out.bam
    ch_bai_raw = SAMTOOLS_SORT_INDEX.out.bai
    ch_bam_bai_raw = SAMTOOLS_SORT_INDEX.out.bam_bai 

    //
    // MODULE: Run SAMtools Stats
    //
    SAMTOOLS_STATS(ch_bam_raw, ch_fasta)
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.stats)
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())

    //
    // MODULE: Run MosDepth
    //
    MOSDEPTH(ch_bam_raw, ch_bai_raw, ch_fasta, params.fai, params.intervals_bed_gunzip, params.intervals_bed_gunzip_index)
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.summary_txt)

    //
    // MODULE: Run ErrorRateByReadPosition 
    //
    FGBIO_ERRORRATEBYREADPOSITION(ch_bam_bai_raw, ch_fasta, ch_fai, ch_dict, params.known_sites, params.known_sites_tbi, params.interval_list)
    ch_versions = ch_versions.mix(FGBIO_ERRORRATEBYREADPOSITION.out.versions.first())

    //
    // MODULE: Run Survivor ScanReads to get Error Profiles
    //
    SURVIVOR_SCAN_READS(ch_bam_bai_raw, params.read_length)
    ch_versions = ch_versions.mix(SURVIVOR_SCAN_READS.out.versions.first())

    //
    // MODULE: Run MSI Sensor PRO
    ///
    MSISENSORPRO(ch_bam_bai_raw, ch_msi_f)
    ch_versions = ch_versions.mix(MSISENSORPRO.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(MSISENSORPRO.out.summary.map{it[1]}.collect())
    ch_multiqc_files = ch_multiqc_files.mix(MSISENSORPRO.out.msi_uns.map{it[1]}.collect())
    ch_multiqc_files = ch_multiqc_files.mix(MSISENSORPRO.out.msi_all.map{it[1]}.collect())

    //
    // MODULE: Run Picard's Collect HS Metrics for raw BAM files
    //
    COLLECTHSMETRICS(ch_bam_bai_raw, ch_fasta, ch_fai, ch_dict, params.baits, params.targets, params.seq_library)
    ch_versions = ch_versions.mix(COLLECTHSMETRICS.out.versions.first())
    ch_coverage_raw  = COLLECTHSMETRICS.out.coverage
    ch_hsmetrics_raw = COLLECTHSMETRICS.out.hsmetrics

    //
    // MODULE: Run Picard Tool CollectMultipleMetrics
    //
    PICARD_COLLECTMULTIPLEMETRICS(ch_bam_raw, ch_fasta, ch_fai)
    ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions.first())

    //
    // MODULE: Run GATK4 MarkDuplicates
    //
    GATK4_MARKDUPLICATES(ch_bam_raw, params.fasta, params.fai)
    ch_multiqc_files = ch_multiqc_files.mix(GATK4_MARKDUPLICATES.out.metrics.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(GATK4_MARKDUPLICATES.out.complex_metrics.collect{it[1]})
    ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES.out.versions.first())
    ch_bam_dedup =  GATK4_MARKDUPLICATES.out.bam

    //
    // MODULE: Run BaseRecalibrator
    //
    GATK4_BASERECALIBRATOR(ch_bam_dedup, params.fasta, params.fai, params.dict, params.known_sites, params.known_sites_tbi)
    ch_versions = ch_versions.mix(GATK4_BASERECALIBRATOR.out.versions.first())
    ch_bqsr_table = GATK4_BASERECALIBRATOR.out.table

    //
    // MODULE: Run ApplyBaseRecalibrator
    //
    GATK4_APPLYBQSR(ch_bam_dedup, ch_bqsr_table, params.fasta, params.fai, params.dict)
    ch_versions = ch_versions.mix(GATK4_APPLYBQSR.out.versions.first())
    ch_bam_recal = GATK4_APPLYBQSR.out.bam

    //
    // MODULE: Run SAMtools Index
    //
    SAMTOOLS_INDEX(ch_bam_recal)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())
    ch_bam_bai_dr = ch_bam_recal.join(SAMTOOLS_INDEX.out.bai, by: 0)

    //
    // MODULE: Run SAMtools Stats
    //
    SAMTOOLS_STATS_DR(ch_bam_bai_dr, ch_fasta)
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS_DR.out.stats)
    ch_versions = ch_versions.mix(SAMTOOLS_STATS_DR.out.versions.first())

    //
    // MODULE: Run MosDepth after MarkDup & ReCal
    //
    MOSDEPTH_DR(ch_bam_bai_dr, ch_fasta, params.fai, params.intervals_bed_gunzip, params.intervals_bed_gunzip_index)
    ch_versions = ch_versions.mix(MOSDEPTH_DR.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH_DR.out.summary_txt)

    //
    // MODULE: Run ErrorRateByReadPosition 
    //
    FGBIO_ERRORRATEBYREADPOSITION_DR(ch_bam_bai_dr, ch_fasta, ch_fai, ch_dict, params.known_sites, params.known_sites_tbi, params.interval_list)
    ch_versions = ch_versions.mix(FGBIO_ERRORRATEBYREADPOSITION_DR.out.versions.first())

    //
    // MODULE: Run Survivor ScanReads to get Error Profiles
    //
    SURVIVOR_SCAN_READS_DR(ch_bam_bai_dr, params.read_length)
    ch_versions = ch_versions.mix(SURVIVOR_SCAN_READS_DR.out.versions.first())

    //
    // MODULE: Run MSI Sensor PRO
    ///
    MSISENSORPRO_DR(ch_bam_bai_dr, ch_msi_f)
    ch_versions = ch_versions.mix(MSISENSORPRO_DR.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(MSISENSORPRO_DR.out.summary.map{it[1]}.collect())
    ch_multiqc_files = ch_multiqc_files.mix(MSISENSORPRO_DR.out.msi_uns.map{it[1]}.collect())
    ch_multiqc_files = ch_multiqc_files.mix(MSISENSORPRO_DR.out.msi_all.map{it[1]}.collect())

    //
    // MODULE: Run Picard's Collect HS Metrics for consensus BAM files
    //
    COLLECTHSMETRICS_DR(ch_bam_bai_dr, ch_fasta, ch_fai, ch_dict, params.baits, params.targets, params.seq_library)
    ch_versions = ch_versions.mix(COLLECTHSMETRICS_DR.out.versions.first())
    ch_coverage_con  = COLLECTHSMETRICS_DR.out.coverage
    ch_hsmetrics_con = COLLECTHSMETRICS_DR.out.hsmetrics

    //
    // MODULE: Run Picard Tool CollectMultipleMetrics
    //
    PICARD_COLLECTMULTIPLEMETRICS_DR(ch_bam_bai_dr, ch_fasta, ch_fai)
    ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS_DR.out.versions.first())

    //
    // MODULE: Extract FastQ reads from BAM
    //
    SAMTOOLS_COLLATEFASTQ(ch_bam_bai_dr, ch_fasta, [])
    ch_versions = ch_versions.mix(SAMTOOLS_COLLATEFASTQ.out.versions)
    ch_reads_dandr = SAMTOOLS_COLLATEFASTQ.out.fastq

    //
    // MODULE: Run FastQC
    //
    FASTQC_DR(ch_reads_dandr)
    ch_versions = ch_versions.mix(FASTQC_DR.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_DR.out.zip.collect{it[1]})

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'software_versions.yml', sort: true, newLine: true)
        .set { ch_vers_coll }

    emit:

    raw_bam         = ch_bam_raw
    raw_bai         = ch_bai_raw
    versions        = ch_vers_coll
    bam_final       = ch_bam_bai_dr
    reads_final     = ch_reads_dandr
    multiqc_files   = ch_multiqc_files 

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                                THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
