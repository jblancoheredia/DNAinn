/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                       PREPROCESSING SUBWORKFLOW                                                    
*******************************************************************************************************************************************************************************************************

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                             IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTP                                                                                                                     } from '../../modules/local/fastp/main'
include { FASTQC                                                                                                                    } from '../../modules/nf-core/fastqc/main'
include { SPADES                                                                                                                    } from '../../modules/local/spades/main'
include { CAT_FASTQ                                                                                                                 } from '../../modules/local/cat/fastq/main'
include { FGBIO_SORT                                                                                                                } from '../../modules/local/fgbio/sort/main'
include { SAMBLASTER                                                                                                                } from '../../modules/local/samblaster/main'
include { ALIGN_CONBAM                                                                                                              } from '../../modules/local/align_con_bam/main'
include { ALIGN_RAWBAM                                                                                                              } from '../../modules/local/align_raw_bam/main'
include { SAMTOOLS_SORT                                                                                                             } from '../../modules/nf-core/samtools/sort/main' 
include { FGBIO_SORTUBAM                                                                                                            } from '../../modules/local/fgbio/sortubam/main'
include { FILTER_CONTIGS                                                                                                            } from '../../modules/local/filter_contigs/main'
include { FGBIO_ZIPFILBAM                                                                                                           } from '../../modules/local/fgbio/zipfilbam/main'
include { FASTQ_CONSENSUS                                                                                                           } from '../../modules/local/fastqc_consensus/main'
include { FGBIO_ZIPRAWBAM                                                                                                           } from '../../modules/local/fgbio/ziprawbam/main'
include { FGBIO_FASTQTOBAM                                                                                                          } from '../../modules/nf-core/fgbio/fastqtobam/main'
include { FGBIO_SORTCONBAM                                                                                                          } from '../../modules/local/fgbio/sortconbam/main.nf'
include { MSISENSORPRO_FIN                                                                                                          } from '../../modules/local/msisensorpro/pro/main'   
include { MSISENSORPRO_RAW                                                                                                          } from '../../modules/local/msisensorpro/pro/main'   
include { SAMTOOLS_SORT_INDEX                                                                                                       } from '../../modules/local/samtools/sort_index/main'   
include { COLLECTHSMETRICS_FIN                                                                                                      } from '../../modules/local/picard/collecthsmetrics/main'
include { COLLECTHSMETRICS_RAW                                                                                                      } from '../../modules/local/picard/collecthsmetrics/main'
include { GATK4_MARKDUPLICATES          	                                                                                        } from '../../modules/local/gatk4/markduplicates/main'
include { FGBIO_GROUPREADSBYUMI                                                                                                     } from '../../modules/local/fgbio/groupreadsbyumi/main'
include { SAMTOOLS_COLLATEFASTQ                                                                                                     } from '../../modules/nf-core/samtools/collatefastq/main'   
include { FGBIO_COLLECTDUPLEXSEQMETRICS                                                                                             } from '../../modules/local/fgbio/collectduplexseqmetrics/main'
include { FGBIO_CALLDUPLEXCONSENSUSREADS                                                                                            } from '../../modules/nf-core/fgbio/callduplexconsensusreads/main'

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

workflow PREPROCESSING {

    take:
    ch_bwa
    ch_fai
    ch_dict
    ch_fasta
    ch_msi_f
    ch_samplesheet
    ch_fastq_single
    ch_fastq_multiple

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (ch_fastq_multiple)
    .reads
    .mix(ch_fastq_single)
    .set { ch_cat_fastq }
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)

    //
    // MODULE: Run FastQC
    //
    FASTQC(ch_cat_fastq)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    //
    // MODULE: Run FastP 
    //
    FASTP(ch_cat_fastq)
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.html.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{it[1]})
    ch_versions = ch_versions.mix(FASTP.out.versions)
    ch_fastp_fastq = FASTP.out.reads

    //
    // MODULE: Run fgbio FastqToBam
    //
    FGBIO_FASTQTOBAM(ch_fastp_fastq)
    ch_versions = ch_versions.mix(FGBIO_FASTQTOBAM.out.versions.first())
    ch_ubam = FGBIO_FASTQTOBAM.out.bam

    //
    // MODULE: Run fgbio SortBam
    //
    FGBIO_SORTUBAM(ch_ubam)
    ch_versions = ch_versions.mix(FGBIO_SORTUBAM.out.versions.first())
    ch_ubam_sorted = FGBIO_SORTUBAM.out.ubam

    //
    // MODULE: Align with bwa mem but avoid sort
    //
    sort = false
    ALIGN_RAWBAM(ch_ubam_sorted, ch_fasta, ch_fai, ch_dict, ch_bwa, sort) //"template-coordinate")
    ch_versions = ch_versions.mix(ALIGN_RAWBAM.out.versions.first())
    ch_raw_bam = ALIGN_RAWBAM.out.bam

    //
    // MODULE: Run SamToools Sort & Index
    //
    SAMTOOLS_SORT_INDEX(ch_raw_bam, ch_fasta, params.fai)
    ch_versions = ch_versions.mix(SAMTOOLS_SORT_INDEX.out.versions.first())
    ch_bam_raw_sorted = SAMTOOLS_SORT_INDEX.out.bam
    ch_bam_raw_indxed = SAMTOOLS_SORT_INDEX.out.bai

    //
    // MODULE: Run Picard's Collect HS Metrics for raw BAM files
    //
    COLLECTHSMETRICS_ORI(ch_bam_raw_sorted, ch_bam_raw_indxed, ch_fasta, ch_fai, ch_dict, params.blocklist_bed, params.targets_Av1)
    ch_versions = ch_versions.mix(COLLECTHSMETRICS_ORI.out.versions.first())
    ch_coverage_raw  = COLLECTHSMETRICS_ORI.out.coverage
    ch_hsmetrics_raw = COLLECTHSMETRICS_ORI.out.hsmetrics

    //
    // MODULE: Run MSI Sensor PRO
    ///
    MSISENSORPRO_PRO_ORI(ch_bam_raw_sorted, ch_bam_raw_indxed, ch_msi_f)
    ch_versions = ch_versions.mix(MSISENSORPRO_PRO_ORI.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(MSISENSORPRO_PRO_ORI.out.summary.map{it[1]}.collect())
    ch_multiqc_files = ch_multiqc_files.mix(MSISENSORPRO_PRO_ORI.out.msi_uns.map{it[1]}.collect())
    ch_multiqc_files = ch_multiqc_files.mix(MSISENSORPRO_PRO_ORI.out.msi_all.map{it[1]}.collect())

    //
    // MODULE: Run SamBlaster
    //
    SAMBLASTER(ch_raw_bam)
    ch_versions = ch_versions.mix(SAMBLASTER.out.versions.first())

    //
    // MODULE: Run fgbio GroupReadsByUmi
    //
    FGBIO_GROUPREADSBYUMI(SAMBLASTER.out.bam, params.group_strategy, params.group_edits, params.group_include_secondary, params.group_include_supplementary, params.group_min_map_q)
    ch_multiqc_files = ch_multiqc_files.mix(FGBIO_GROUPREADSBYUMI.out.histogram.map{it[1]}.collect())
    ch_versions = ch_versions.mix(FGBIO_GROUPREADSBYUMI.out.versions.first())
    ch_bam_grouped = FGBIO_GROUPREADSBYUMI.out.bam

    //
    // MODULE: Run SPAdes
    //
    SPADES(SAMBLASTER.out.split_reads)
    ch_versions = ch_versions.mix(SPADES.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(SPADES.out.log.map{it[1]}.collect())

    //
    // MODULE: Run FilterContigs custom script
    //
    FILTER_CONTIGS(SPADES.out.contigs, '100', '1')
    ch_versions = ch_versions.mix(FILTER_CONTIGS.out.versions.first())
    ch_split_contigs = FILTER_CONTIGS.out.fasta

    //
    // MODULE: Run fgbio CollectDuplexSeqMetrics
    //
    FGBIO_COLLECTDUPLEXSEQMETRICS(ch_bam_grouped, params.intervals_Av1)
    ch_multiqc_files = ch_multiqc_files.mix(FGBIO_COLLECTDUPLEXSEQMETRICS.out.metrics.map{it[1]}.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FGBIO_COLLECTDUPLEXSEQMETRICS.out.pdf.map{it[1]}.collect())   
    ch_versions = ch_versions.mix(FGBIO_COLLECTDUPLEXSEQMETRICS.out.versions.first())

    //
    // MODULE: Run fgbio CallDuplexConsensusReads
    //
    FGBIO_CALLDUPLEXCONSENSUSREADS(ch_bam_grouped, params.call_min_reads, params.call_min_baseq)
    ch_versions = ch_versions.mix(FGBIO_CALLDUPLEXCONSENSUSREADS.out.versions.first())
    FGBIO_CALLDUPLEXCONSENSUSREADS.out.bam.set { ch_consensus_bam }

    //
    // MODULE: Run fgbio SortBam
    //
    FGBIO_SORTCONBAM(ch_consensus_bam)
    ch_versions = ch_versions.mix(FGBIO_SORTCONBAM.out.versions.first())
    ch_consensus_bam_sorted = FGBIO_SORTCONBAM.out.bam

    //
    // MODULE: Align with BWA mem
    //
    ALIGN_CONBAM(ch_consensus_bam_sorted, ch_fasta, ch_fai, ch_dict, ch_bwa, "none")
    ch_versions = ch_versions.mix(ALIGN_CONBAM.out.versions.first())

    //
    // MODULE: Run SamTools Sort
    //
    SAMTOOLS_SORT(ALIGN_CONBAM.out.bam, ch_fasta, params.fai)
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    //
    // MODULE: Run GATK4 MarkDuplicates
    //
    GATK4_MARKDUPLICATES(SAMTOOLS_SORT.out.bam, params.fasta, params.fai)
    ch_multiqc_files = ch_multiqc_files.mix(GATK4_MARKDUPLICATES.out.metrics.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(GATK4_MARKDUPLICATES.out.complex_metrics.collect{it[1]})
    ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES.out.versions.first())
    ch_bam_consensus_dedup =  GATK4_MARKDUPLICATES.out.bam

    //
    // MODULE: Run Picard's Collect HS Metrics for consensus BAM files
    //
    COLLECTHSMETRICS_FIN(ch_bam_consensus_dedup, ch_fasta, ch_fai, ch_dict, params.blocklist_bed, params.targets_Av1)
    ch_versions = ch_versions.mix(COLLECTHSMETRICS_FIN.out.versions.first())
    ch_coverage_con  = COLLECTHSMETRICS_FIN.out.coverage
    ch_hsmetrics_con = COLLECTHSMETRICS_FIN.out.hsmetrics

    //
    // MODULE: Run MSI Sensor PRO
    ///
    MSISENSORPRO_PRO_FIN(ch_bam_consensus_dedup, ch_msi_f)
    ch_versions = ch_versions.mix(MSISENSORPRO_PRO_FIN.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(MSISENSORPRO_PRO_FIN.out.summary.map{it[1]}.collect())
    ch_multiqc_files = ch_multiqc_files.mix(MSISENSORPRO_PRO_FIN.out.msi_uns.map{it[1]}.collect())
    ch_multiqc_files = ch_multiqc_files.mix(MSISENSORPRO_PRO_FIN.out.msi_all.map{it[1]}.collect())

    //
    // MODULE: Extract FastQ reads from BAM
    //
    SAMTOOLS_COLLATEFASTQ(ch_bam_consensus_dedup, ch_fasta, [])
    ch_versions = ch_versions.mix(SAMTOOLS_COLLATEFASTQ.out.versions)
    ch_consensus_reads = SAMTOOLS_COLLATEFASTQ.out.fastq

    //
    // MODULE: Run FastQC
    //
    FASTQ_CONSENSUS(ch_consensus_reads)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_CONSENSUS.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQ_CONSENSUS.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'software_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    emit:

    ubam            = ch_ubam
    raw_bam         = ch_raw_bam
    versions        = ch_collated_versions
    raw_reads       = ch_cat_fastq
    group_bam       = ch_bam_grouped
    multiqc_files   = ch_multiqc_files 
    consensus_bam   = ch_consensus_bam
    reads_consens   = ch_consensus_reads
    split_contigs   = ch_split_contigs
    cons_bam_dedup  = ch_bam_consensus_dedup
    split_reads     = SAMBLASTER.out.split_reads
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                                THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
