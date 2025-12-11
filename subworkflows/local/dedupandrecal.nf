/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                       DEDUPANDRECAL SUBWORKFLOW                                                    
*******************************************************************************************************************************************************************************************************

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                             IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MOSDEPTH                                                                                                                  } from '../../modules/local/mosdepth/main'
include { BWAMEM2_MEM                                                                                                               } from '../../modules/nf-core/bwamem2/mem/main' 
include { MOSDEPTH_DR                                                                                                               } from '../../modules/local/mosdepth/main'
include { SAMTOOLS_SORT                                                                                                             } from '../../modules/nf-core/samtools/sort/main' 
include { SAMTOOLS_INDEX                                                                                                            } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_STATS                                                                                                            } from '../../modules/nf-core/samtools/stats/main'
include { FASTQ_CONSENSUS                                                                                                           } from '../../modules/local/fastqc_consensus/main'
include { GATK4_APPLYBQSR                                                                                                           } from '../../modules/local/gatk4/applybqsr/main'
include { MSISENSORPRO_CON                                                                                                          } from '../../modules/local/msisensorpro/pro/main'   
include { MSISENSORPRO_RAW                                                                                                          } from '../../modules/local/msisensorpro/pro/main'   
include { SAMTOOLS_STATS_DR                                                                                                         } from '../../modules/local/samtools/stats_dr/main'
include { SAMTOOLS_SORT_INDEX                                                                                                       } from '../../modules/local/samtools/sort_index/main'
include { SURVIVOR_SCAN_READS                                                                                                       } from '../../modules/local/survivor/scanreads/main'
include { COLLECTHSMETRICS_CON                                                                                                      } from '../../modules/local/picard/collecthsmetrics/main'
include { COLLECTHSMETRICS_RAW                                                                                                      } from '../../modules/local/picard/collecthsmetrics/main'
include { GATK4_MARKDUPLICATES          	                                                                                        } from '../../modules/local/gatk4/markduplicates/main'
include { SAMTOOLS_COLLATEFASTQ                                                                                                     } from '../../modules/nf-core/samtools/collatefastq/main'   
include { GATK4_BASERECALIBRATOR                                                                                                    } from '../../modules/local/gatk4/baserecalibrator/main'
include { FGBIO_ERRORRATEBYREADPOSITION                                                                                             } from '../../modules/local/fgbio/errorratebyreadposition/main'

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
    ch_fastp_fastq

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: BWA-MEM2 mapping
    //
    sort_bam = 'sort'
    BWAMEM2_MEM(ch_fastp_fastq, ch_bwa2, ch_fasta, ch_fai, sort_bam)
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions)
    ch_bam = BWAMEM2_MEM.out.bam

    //
    // MODULE: Run SamToools Sort & Index
    //
    SAMTOOLS_SORT_INDEX(ch_bam, ch_fasta, params.fai)
    ch_versions = ch_versions.mix(SAMTOOLS_SORT_INDEX.out.versions.first())
    ch_bam_raw = SAMTOOLS_SORT_INDEX.out.bam
    ch_bam_bai = SAMTOOLS_SORT_INDEX.out.bam_bai 
    ch_bam_raw_index = SAMTOOLS_SORT_INDEX.out.bai

    //
    // MODULE: Run Survivor ScanReads to get Error Profiles
    //
    SURVIVOR_SCAN_READS(ch_bam_bai, params.read_length)
    ch_versions = ch_versions.mix(SURVIVOR_SCAN_READS.out.versions.first())

    //
    // MODULE: Run Picard's Collect HS Metrics for raw BAM files
    //
    COLLECTHSMETRICS_RAW(ch_bam_raw, ch_bam_raw_index, ch_fasta, ch_fai, ch_dict, params.blocklist_bed, params.targets)
    ch_versions = ch_versions.mix(COLLECTHSMETRICS_RAW.out.versions.first())
    ch_coverage_raw  = COLLECTHSMETRICS_RAW.out.coverage
    ch_hsmetrics_raw = COLLECTHSMETRICS_RAW.out.hsmetrics

    //
    // MODULE: Run MSI Sensor PRO
    ///
    MSISENSORPRO_RAW(ch_bam_bai, ch_msi_f)
    ch_versions = ch_versions.mix(MSISENSORPRO_RAW.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(MSISENSORPRO_RAW.out.summary.map{it[1]}.collect())
    ch_multiqc_files = ch_multiqc_files.mix(MSISENSORPRO_RAW.out.msi_uns.map{it[1]}.collect())
    ch_multiqc_files = ch_multiqc_files.mix(MSISENSORPRO_RAW.out.msi_all.map{it[1]}.collect())

    //
    // MODULE: Run MosDepth
    //
    MOSDEPTH(ch_bam_raw, ch_bam_raw_index, ch_fasta, params.fai, params.intervals_bed_gunzip, params.intervals_bed_gunzip_index)
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.summary_txt)

    //
    // MODULE: Run GATK4 MarkDuplicates
    //
    GATK4_MARKDUPLICATES(ch_bam_raw, params.fasta, params.fai)
    ch_multiqc_files = ch_multiqc_files.mix(GATK4_MARKDUPLICATES.out.metrics.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(GATK4_MARKDUPLICATES.out.complex_metrics.collect{it[1]})
    ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES.out.versions.first())
    ch_bam_dedup =  GATK4_MARKDUPLICATES.out.bam

    //
    // MODULE: Run SAMtools Stats
    //
    SAMTOOLS_STATS(ch_bam_dedup, ch_fasta)
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.stats)
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())

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
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.stats)
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())

    //
    // MODULE: Run MosDepth after MarkDup & ReCal
    //
    MOSDEPTH_DR(ch_bam_bai_dr, ch_fasta, params.fai, params.intervals_bed_gunzip, params.intervals_bed_gunzip_index)
    ch_versions = ch_versions.mix(MOSDEPTH_DR.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH_DR.out.summary_txt)

    //
    // MODULE: Run Picard's Collect HS Metrics for consensus BAM files
    //
    COLLECTHSMETRICS_CON(ch_bam_bai_dr, ch_fasta, ch_fai, ch_dict, params.blocklist_bed, params.baits, params.targets)
    ch_versions = ch_versions.mix(COLLECTHSMETRICS_CON.out.versions.first())
    ch_coverage_con  = COLLECTHSMETRICS_CON.out.coverage
    ch_hsmetrics_con = COLLECTHSMETRICS_CON.out.hsmetrics

    //
    // MODULE: Run MSI Sensor PRO
    ///
    MSISENSORPRO_CON(ch_bam_bai_dr, ch_msi_f)
    ch_versions = ch_versions.mix(MSISENSORPRO_CON.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(MSISENSORPRO_CON.out.summary.map{it[1]}.collect())
    ch_multiqc_files = ch_multiqc_files.mix(MSISENSORPRO_CON.out.msi_uns.map{it[1]}.collect())
    ch_multiqc_files = ch_multiqc_files.mix(MSISENSORPRO_CON.out.msi_all.map{it[1]}.collect())

    //
    // MODULE: Extract FastQ reads from BAM
    //
    SAMTOOLS_COLLATEFASTQ(ch_bam_bai_dr, ch_fasta, [])
    ch_versions = ch_versions.mix(SAMTOOLS_COLLATEFASTQ.out.versions)
    ch_reads_dr = SAMTOOLS_COLLATEFASTQ.out.fastq

    //
    // MODULE: Run FastQC
    //
    FASTQ_CONSENSUS(ch_reads_dr)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_CONSENSUS.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQ_CONSENSUS.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'software_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    emit:

    raw_bam         = ch_bam_raw
    raw_bai         = ch_bai_raw
    versions        = ch_collated_versions
    bam_final       = ch_bam_bai_dr
    reads_final     = ch_reads_dr
    multiqc_files   = ch_multiqc_files 

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                                THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
