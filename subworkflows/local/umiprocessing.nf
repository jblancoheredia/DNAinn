/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                       UMIPROCESSING SUBWORKFLOW                                                    
*******************************************************************************************************************************************************************************************************

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                             IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BAMCUT                                                                                                                    } from '../../modules/local/bamcut/main' // <- In Beta
include { SPADES                                                                                                                    } from '../../modules/local/spades/main' // <- In use
include { REPEATSEQ                                                                                                                 } from '../../modules/local/repeatseq/main' // <- In Beta
include { MERGE_REPS                                                                                                                } from '../../modules/local/merge_reps/main' // <- In Beta
include { SAMBLASTER                                                                                                                } from '../../modules/local/samblaster/main' // <- In use
include { MOSDEPTH_DUP                                                                                                              } from '../../modules/local/mosdepth/main' // <- In use
include { MOSDEPTH_CON                                                                                                              } from '../../modules/local/mosdepth/main' // <- In use
include { MOSDEPTH_RAW                                                                                                              } from '../../modules/local/mosdepth/main' // <- In use
include { MOSDEPTH_SIM                                                                                                              } from '../../modules/local/mosdepth/main' // <- In use
include { ALIGN_BAM_CON                                                                                                             } from '../../modules/local/umi_align_bam/main' // <- In use
include { ALIGN_BAM_RAW                                                                                                             } from '../../modules/local/umi_align_bam/main' // <- In use
include { PRESEQ_CCURVE                                                                                                             } from '../../modules/local/preseq/ccurve/main' // <- In use
include { FILTER_CONTIGS                                                                                                            } from '../../modules/local/filter_contigs/main' // <- In use
include { PRESEQ_LCEXTRAP                                                                                                           } from '../../modules/local/preseq/lcextrap/main' // <- In use
include { UMI_READ_COUNTS                                                                                                           } from '../../modules/local/umi_read_counts/main'
include { FASTQC_CONSENSUS                                                                                                          } from '../../modules/local/fastqc_consensus/main' // <- In use
include { FGBIO_FASTQTOBAM                                                                                                          } from '../../modules/nf-core/fgbio/fastqtobam/main' // <- In use
include { FGBIO_SORTCONBAM                                                                                                          } from '../../modules/local/fgbio/sortconbam/main.nf' // <- In use
include { MSISENSORPRO_DUP                                                                                                          } from '../../modules/local/msisensorpro/pro/main' // <- In use
include { MSISENSORPRO_CON                                                                                                          } from '../../modules/local/msisensorpro/pro/main' // <- In use
include { MSISENSORPRO_RAW                                                                                                          } from '../../modules/local/msisensorpro/pro/main' // <- In use
include { MSISENSORPRO_SIM                                                                                                          } from '../../modules/local/msisensorpro/pro/main' // <- In use
include { FGBIO_CORRECTUMIS                                                                                                         } from '../../modules/local/fgbio/correctumis/main' // <- New in use
include { SAMTOOLS_STATS_CON                                                                                                        } from '../../modules/local/samtools/stats/main'
include { SAMTOOLS_STATS_DUP                                                                                                        } from '../../modules/local/samtools/stats/main'
include { SAMTOOLS_STATS_RAW                                                                                                        } from '../../modules/local/samtools/stats/main'
include { SAMTOOLS_STATS_SIM                                                                                                        } from '../../modules/local/samtools/stats/main'
include { COLLECT_UMI_METRICS                                                                                                       } from '../../modules/local/collect_umi_metrics/main'
include { COLLECTHSMETRICS_CON                                                                                                      } from '../../modules/local/picard/collecthsmetrics/main' // <- In use
include { COLLECTHSMETRICS_DUP                                                                                                      } from '../../modules/local/picard/collecthsmetrics/main' // <- In use
include { COLLECTHSMETRICS_RAW                                                                                                      } from '../../modules/local/picard/collecthsmetrics/main' // <- In use
include { COLLECTHSMETRICS_SIM                                                                                                      } from '../../modules/local/picard/collecthsmetrics/main' // <- In use
include { FGBIO_GROUPREADSBYUMI                                                                                                     } from '../../modules/local/fgbio/groupreadsbyumi/main' // <- In use
include { SAMTOOLS_COLLATEFASTQ                                                                                                     } from '../../modules/nf-core/samtools/collatefastq/main' // <- In use
include { SAMTOOLS_SORT_INDEX_CON                                                                                                   } from '../../modules/local/samtools/sort_index/main' // <- In use
include { SAMTOOLS_SORT_INDEX_RAW                                                                                                   } from '../../modules/local/samtools/sort_index/main' // <- In use
include { SURVIVOR_SCAN_READS_CON                                                                                                   } from '../../modules/local/survivor/scanreads/main' // <- In use
include { SURVIVOR_SCAN_READS_DUP                                                                                                   } from '../../modules/local/survivor/scanreads/main' // <- In use
include { SURVIVOR_SCAN_READS_RAW                                                                                                   } from '../../modules/local/survivor/scanreads/main' // <- In use
include { SURVIVOR_SCAN_READS_SIM                                                                                                   } from '../../modules/local/survivor/scanreads/main' // <- In use
include { FGBIO_FILTERCONSENSUSREADS                                                                                                } from '../../modules/local/fgbio/filterconsensusreads/main' // <- New in use
include { FGBIO_COLLECTDUPLEXSEQMETRICS                                                                                             } from '../../modules/local/fgbio/collectduplexseqmetrics/main' // <- In use
include { PICARD_COLLECTMULTIPLEMETRICS                                                                                             } from '../../modules/local/picard/collectmultiplemetrics/main' // <- In use
include { FGBIO_CALLDUPLEXCONSENSUSREADS                                                                                            } from '../../modules/nf-core/fgbio/callduplexconsensusreads/main' // <- In use
include { FGBIO_ERRORRATEBYREADPOSITION_CON                                                                                         } from '../../modules/local/fgbio/errorratebyreadposition/main' // <- In use
include { FGBIO_ERRORRATEBYREADPOSITION_DUP                                                                                         } from '../../modules/local/fgbio/errorratebyreadposition/main' // <- In use
include { FGBIO_ERRORRATEBYREADPOSITION_RAW                                                                                         } from '../../modules/local/fgbio/errorratebyreadposition/main' // <- In use
include { FGBIO_ERRORRATEBYREADPOSITION_SIM                                                                                         } from '../../modules/local/fgbio/errorratebyreadposition/main' // <- In use

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

workflow UMIPROCESSING {

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
    // MODULE: Run fgbio FastqToBam
    //
    FGBIO_FASTQTOBAM(ch_fastqs)
    ch_versions = ch_versions.mix(FGBIO_FASTQTOBAM.out.versions.first())
    ch_ubam = FGBIO_FASTQTOBAM.out.bam

    //
    // MODULE: Align with bwa mem but avoid sort
    //
    sort = false
    ALIGN_BAM_RAW(ch_ubam, ch_fasta, ch_fai, ch_dict, ch_bwa2, sort)
    ch_versions = ch_versions.mix(ALIGN_BAM_RAW.out.versions.first())
    ch_raw_bam = ALIGN_BAM_RAW.out.bam
    ch_raw_sort_bam = ALIGN_BAM_RAW.out.sort_bam
    ch_raw_sort_bai = ALIGN_BAM_RAW.out.sort_bai

    //
    // MODULE: Run fgbio correctumis
    //
    FGBIO_CORRECTUMIS(ch_raw_bam, params.correct_max_mismatch, params.correct_min_distance, params.correct_min_corrected)
    ch_multiqc_files = ch_multiqc_files.mix(FGBIO_CORRECTUMIS.out.metrics)
    ch_versions = ch_versions.mix(FGBIO_CORRECTUMIS.out.versions.first())
    ch_bam_fcu = FGBIO_CORRECTUMIS.out.bam

    //
    // MODULE: Run SamToools Sort & Index
    //
    SAMTOOLS_SORT_INDEX_RAW(ch_bam_fcu, ch_fasta, params.fai)
    ch_versions = ch_versions.mix(SAMTOOLS_SORT_INDEX_RAW.out.versions.first())
    ch_bam_fcu_sort = SAMTOOLS_SORT_INDEX_RAW.out.bam
    ch_bam_fcu_indx = SAMTOOLS_SORT_INDEX_RAW.out.bai
    ch_bam_fcu_stix = SAMTOOLS_SORT_INDEX_RAW.out.bam_bai

    //
    // MODULE: Run Picard Tool CollectMultipleMetrics
    //
    PICARD_COLLECTMULTIPLEMETRICS(ch_bam_fcu_stix, ch_fasta, ch_fai)
    ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions.first())

    //
    // MODULE: Run ErrorRateByReadPosition 
    //
    FGBIO_ERRORRATEBYREADPOSITION_RAW(ch_bam_fcu_sort, ch_fasta, ch_fai, ch_dict, params.known_sites, params.known_sites_tbi, params.interval_list)
    ch_versions = ch_versions.mix(FGBIO_ERRORRATEBYREADPOSITION_RAW.out.versions.first())

    //
    // MODULE: Run MosDepth
    //
    MOSDEPTH_RAW(ch_bam_fcu_stix, ch_fasta, params.fai, params.intervals_bed_gunzip, params.intervals_bed_gunzip_index)
    ch_versions = ch_versions.mix(MOSDEPTH_RAW.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH_RAW.out.summary_txt)

    //
    // MODULE: Run Survivor ScanReads to get Error Profiles
    //
    SURVIVOR_SCAN_READS_RAW(ch_bam_fcu_stix, params.read_length)
    ch_versions = ch_versions.mix(SURVIVOR_SCAN_READS_RAW.out.versions.first())

    //
    // MODULE: Run Picard's Collect HS Metrics for raw BAM files
    //
    COLLECTHSMETRICS_RAW(ch_bam_fcu_stix, ch_fasta, ch_fai, ch_dict, params.hsmetrics_baits, params.hsmetrics_trgts, params.seq_library)
    ch_versions = ch_versions.mix(COLLECTHSMETRICS_RAW.out.versions.first())
    ch_coverage_raw  = COLLECTHSMETRICS_RAW.out.coverage
    ch_hsmetrics_raw = COLLECTHSMETRICS_RAW.out.hsmetrics

    //
    // MODULE: Run MSI Sensor PRO
    ///
    MSISENSORPRO_RAW(ch_bam_fcu_stix, ch_msi_f)
    ch_versions = ch_versions.mix(MSISENSORPRO_RAW.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(MSISENSORPRO_RAW.out.summary.map{it[1]}.collect())
    ch_multiqc_files = ch_multiqc_files.mix(MSISENSORPRO_RAW.out.msi_uns.map{it[1]}.collect())
    ch_multiqc_files = ch_multiqc_files.mix(MSISENSORPRO_RAW.out.msi_all.map{it[1]}.collect())

//    //
//    // MODULE: Run BAMCUT to split by Chromosome
//    //
//    BAMCUT(ch_bam_fcu_stix)
//    ch_versions = ch_versions.mix(BAMCUT.out.versions)
//    ch_bams = BAMCUT.out.bams
//    ch_bais = BAMCUT.out.bais
//
//    ch_bam_by_chrom = ch_bams.flatMap { meta, bamList ->
//        bamList.collect { bam ->
//            def chrom = bam.name.replace(meta.id + "_", "").replace(".bam", "")
//            tuple(meta, chrom, bam)
//        }
//    }
//    
//    ch_bai_by_chrom = ch_bais.flatMap { meta, baiList ->
//        baiList.collect { bai ->
//            def chrom = bai.name.replace(meta.id + "_", "").replace(".bai", "")
//            tuple(meta, chrom, bai)
//        }
//    }
//
//    //
//    // MODULE: Run RepeatSeq
//    //
//    REPEATSEQ(ch_bam_by_chrom, ch_bai_by_chrom, ch_fasta, ch_fai, params.rep_regions)
//    ch_repeats_vcfs = REPEATSEQ.out.vcf.collect()
//    ch_repeats_calls = REPEATSEQ.out.calls.collect()
//    ch_repeats_repeatseqs = REPEATSEQ.out.repeatseq.collect()
//    ch_versions = ch_versions.mix(REPEATSEQ.out.versions.first())
//    REPEATSEQ.out.vcf
//        .map { meta, file -> [meta.id, meta, file] }  
//        .groupTuple()
//        .map { id, metas, files -> [metas[0], files] }
//        .branch {
//            meta, files -> 
//                complete: files.size() > 0
//                    return [meta, files]
//                failed: true
//                    return meta
//        }
//        .set { grouped_vcfs }
//    REPEATSEQ.out.calls
//        .map { meta, file -> [meta.id, meta, file] }  
//        .groupTuple()
//        .map { id, metas, files -> [metas[0], files] }
//        .branch {
//            meta, files -> 
//                complete: files.size() > 0
//                    return [meta, files]
//                failed: true
//                    return meta
//        }
//        .set { grouped_calls }
//    REPEATSEQ.out.repeatseq
//        .map { meta, file -> [meta.id, meta, file] }  
//        .groupTuple()
//        .map { id, metas, files -> [metas[0], files] }
//        .branch {
//            meta, files -> 
//                complete: files.size() > 0
//                    return [meta, files]
//                failed: true
//                    return meta
//        }
//        .set { grouped_repeatseq }
//
//    //
//    // MODULE: Run MergeBams to collect the outputs from RepeatSeq
//    //
//    grouped_vcfs.complete
//        .join(grouped_calls.complete)
//        .join(grouped_repeatseq.complete)
//        .set { complete_sets }
//    MERGE_REPS(
//        complete_sets.map { meta, vcfs, calls, repeatseqs -> [meta, vcfs] },
//        complete_sets.map { meta, vcfs, calls, repeatseqs -> [meta, calls] },
//        complete_sets.map { meta, vcfs, calls, repeatseqs -> [meta, repeatseqs] }
//    )

    //
    // MODULE: Run SamBlaster
    //
    SAMBLASTER(ch_bam_fcu)
    ch_versions = ch_versions.mix(SAMBLASTER.out.versions.first())
    ch_split_reads = SAMBLASTER.out.split_reads

    //
    // MODULE: Run SPAdes
    //
    SPADES(ch_split_reads)
    ch_versions = ch_versions.mix(SPADES.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(SPADES.out.log.map{it[1]}.collect())

    //
    // MODULE: Run FilterContigs custom script
    //
    FILTER_CONTIGS(SPADES.out.contigs, '100', '1')
    ch_versions = ch_versions.mix(FILTER_CONTIGS.out.versions.first())
    ch_split_contigs = FILTER_CONTIGS.out.fasta
    
    //
    // MODULE: Run fgbio GroupReadsByUmi
    //
    FGBIO_GROUPREADSBYUMI(SAMBLASTER.out.bam, params.group_strategy, params.group_edits, params.group_include_secondary, params.group_allow_inter_contig, params.group_include_supplementary, params.group_min_map_q, params.group_include_non_pf_reads, params.group_mark_duplicates)
    ch_multiqc_files = ch_multiqc_files.mix(FGBIO_GROUPREADSBYUMI.out.histogram.map{it[1]}.collect())
    ch_versions = ch_versions.mix(FGBIO_GROUPREADSBYUMI.out.versions.first())
    ch_grouped_family_sizes = FGBIO_GROUPREADSBYUMI.out.histogram
    ch_bam_grouped = FGBIO_GROUPREADSBYUMI.out.bam

    //
    // MODULE: Run fgbio CollectDuplexSeqMetrics
    //
    FGBIO_COLLECTDUPLEXSEQMETRICS(ch_bam_grouped, params.interval_list)
    ch_multiqc_files = ch_multiqc_files.mix(FGBIO_COLLECTDUPLEXSEQMETRICS.out.metrics.map{it[1]}.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FGBIO_COLLECTDUPLEXSEQMETRICS.out.pdf.map{it[1]}.collect())   
    ch_versions = ch_versions.mix(FGBIO_COLLECTDUPLEXSEQMETRICS.out.versions.first())

    //
    // MODULE: Run fgbio CallDuplexConsensusReads
    //
    FGBIO_CALLDUPLEXCONSENSUSREADS(ch_bam_grouped, params.call_min_reads, params.call_min_baseq)
    ch_versions = ch_versions.mix(FGBIO_CALLDUPLEXCONSENSUSREADS.out.versions.first())
    FGBIO_CALLDUPLEXCONSENSUSREADS.out.bam.set { ch_bam_consensus }

    //
    // MODULE: Run fgbio SortBam
    //
    FGBIO_SORTCONBAM(ch_bam_consensus)
    ch_versions = ch_versions.mix(FGBIO_SORTCONBAM.out.versions.first())
    ch_bam_consensus_sorted = FGBIO_SORTCONBAM.out.bam

    //
    // MODULE: Run FgBIO FilterConsensusReads to produce the "Consensus", "Duplex" & "Simplex" BAM files
    //
    FGBIO_FILTERCONSENSUSREADS(ch_bam_consensus_sorted, params.fasta, params.fai, params.filter_min_reads, params.filter_min_base_quality, params.filter_max_base_error_rate, params.filter_max_read_error_rate, params.filter_max_no_call_fraction)
    ch_versions = ch_versions.mix(FGBIO_FILTERCONSENSUSREADS.out.versions.first())
    ch_bam_bai_con_fil = FGBIO_FILTERCONSENSUSREADS.out.suplex_bam_bai
    ch_bam_bai_duplex_fil = FGBIO_FILTERCONSENSUSREADS.out.duplex_bam_bai
    ch_bam_bai_simplex_fil = FGBIO_FILTERCONSENSUSREADS.out.simplex_bam_bai

    // Combine BAM fils by meta data
	ch_align_bam_con_in = ch_bam_bai_con_fil
	    .join(ch_bam_bai_duplex_fil)
	    .join(ch_bam_bai_simplex_fil)

    //
    // MODULE: Align with BWA mem
    //
    ALIGN_BAM_CON(ch_align_bam_con_in, ch_fasta, ch_fai, ch_dict, ch_bwa2)
    ch_versions = ch_versions.mix(ALIGN_BAM_CON.out.versions.first())
    ch_bam_con = ALIGN_BAM_CON.out.bam
    ch_bam_duplex = ALIGN_BAM_CON.out.duplex_bam
    ch_bam_simplex = ALIGN_BAM_CON.out.simplex_bam

    // Combine BAM fils by meta data
	ch_sort_index_in = ch_bam_con
	    .join(ch_bam_duplex)
	    .join(ch_bam_simplex)

    //
    // MODULE: Run SamToools Sort & Index
    //
    SAMTOOLS_SORT_INDEX_CON(ch_sort_index_in, ch_fasta, ch_fai)
    ch_versions = ch_versions.mix(SAMTOOLS_SORT_INDEX_CON.out.versions)
    ch_bam_dup_stix = SAMTOOLS_SORT_INDEX_CON.out.bam_duplex
    ch_bam_sim_stix = SAMTOOLS_SORT_INDEX_CON.out.bam_simplex
    ch_bam_con_stix = SAMTOOLS_SORT_INDEX_CON.out.bam_consensus

    //
    // MODULE: Run Survivor ScanReads to get Consensus Error Profiles
    //
    SURVIVOR_SCAN_READS_CON(ch_bam_con_stix, params.read_length)
    ch_versions = ch_versions.mix(SURVIVOR_SCAN_READS_CON.out.versions.first())

    //
    // MODULE: Run Survivor ScanReads to get Duplex Error Profiles
    //
    SURVIVOR_SCAN_READS_DUP(ch_bam_dup_stix, params.read_length)
    ch_versions = ch_versions.mix(SURVIVOR_SCAN_READS_DUP.out.versions.first())

    //
    // MODULE: Run Survivor ScanReads to get Simplex Error Profiles
    //
    SURVIVOR_SCAN_READS_SIM(ch_bam_sim_stix, params.read_length)
    ch_versions = ch_versions.mix(SURVIVOR_SCAN_READS_SIM.out.versions.first())

    // Combine BAM fils by meta data
	ch_umi_metrics_in = ch_bam_con_stix
	    .join(ch_bam_dup_stix)
	    .join(ch_bam_sim_stix)

    //
    // MODULE: Run SamTools View to count reads accross the BAM files
    //
    COLLECT_UMI_METRICS(ch_umi_metrics_in)
    ch_versions = ch_versions.mix(COLLECT_UMI_METRICS.out.versions.first())
    ch_con_family_sizes = COLLECT_UMI_METRICS.out.con_family_sizes

    // Combine BAM fils by meta data
	ch_umi_read_counts_in = ch_ubam
	    .join(ch_bam_fcu)
	    .join(ch_bam_grouped)
	    .join(ch_bam_consensus)
	    .join(ch_bam_bai_con_fil)
	    .join(ch_bam_con_stix)
	    .join(ch_bam_dup_stix)
	    .join(ch_bam_sim_stix)

    //
    // MODULE: Run SamTools View to count reads accross the BAM files
    //
    UMI_READ_COUNTS(ch_umi_read_counts_in)
    ch_versions = ch_versions.mix(UMI_READ_COUNTS.out.versions.first())

    //
    // MODULE: Run Preseq CCurve
    //
    PRESEQ_CCURVE(ch_con_family_sizes)
    ch_versions = ch_versions.mix(PRESEQ_CCURVE.out.versions.first())

    //
    // MODULE: Run Preseq LCExtrap
    //
    PRESEQ_LCEXTRAP(ch_grouped_family_sizes)
    ch_versions = ch_versions.mix(PRESEQ_LCEXTRAP.out.versions.first())

    //
    // MODULE: Run MosDepth
    //
    MOSDEPTH_CON(ch_bam_con_stix, ch_fasta, params.fai, params.intervals_bed_gunzip, params.intervals_bed_gunzip_index)
    ch_versions = ch_versions.mix(MOSDEPTH_CON.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH_CON.out.summary_txt)

    //
    // MODULE: Run MosDepth
    //
    MOSDEPTH_DUP(ch_bam_dup_stix, ch_fasta, params.fai, params.intervals_bed_gunzip, params.intervals_bed_gunzip_index)
    ch_versions = ch_versions.mix(MOSDEPTH_DUP.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH_DUP.out.summary_txt)

    //
    // MODULE: Run MosDepth
    //
    MOSDEPTH_SIM(ch_bam_sim_stix, ch_fasta, params.fai, params.intervals_bed_gunzip, params.intervals_bed_gunzip_index)
    ch_versions = ch_versions.mix(MOSDEPTH_SIM.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH_SIM.out.summary_txt)

    //
    // MODULE: Run ErrorRateByReadPosition in Final BAM
    //
    FGBIO_ERRORRATEBYREADPOSITION_CON(ch_bam_con_stix, ch_fasta, ch_fai, ch_dict, params.known_sites, params.known_sites_tbi, params.interval_list)
    ch_versions = ch_versions.mix(FGBIO_ERRORRATEBYREADPOSITION_CON.out.versions)

    //
    // MODULE: Run ErrorRateByReadPosition in Final BAM
    //
    FGBIO_ERRORRATEBYREADPOSITION_DUP(ch_bam_dup_stix, ch_fasta, ch_fai, ch_dict, params.known_sites, params.known_sites_tbi, params.interval_list)
    ch_versions = ch_versions.mix(FGBIO_ERRORRATEBYREADPOSITION_DUP.out.versions)

    //
    // MODULE: Run ErrorRateByReadPosition in Final BAM
    //
    FGBIO_ERRORRATEBYREADPOSITION_SIM(ch_bam_sim_stix, ch_fasta, ch_fai, ch_dict, params.known_sites, params.known_sites_tbi, params.interval_list)
    ch_versions = ch_versions.mix(FGBIO_ERRORRATEBYREADPOSITION_SIM.out.versions)

    //
    // MODULE: Run Picard's Collect HS Metrics for consensus BAM files
    //
    COLLECTHSMETRICS_CON(ch_bam_con_stix, ch_fasta, ch_fai, ch_dict, params.hsmetrics_baits, params.hsmetrics_trgts, params.seq_library)
    ch_versions = ch_versions.mix(COLLECTHSMETRICS_CON.out.versions.first())
    ch_coverage_con  = COLLECTHSMETRICS_CON.out.coverage
    ch_hsmetrics_con = COLLECTHSMETRICS_CON.out.hsmetrics

    //
    // MODULE: Run Picard's Collect HS Metrics for consensus BAM files
    //
    COLLECTHSMETRICS_DUP(ch_bam_dup_stix, ch_fasta, ch_fai, ch_dict, params.hsmetrics_baits, params.hsmetrics_trgts, params.seq_library)
    ch_versions = ch_versions.mix(COLLECTHSMETRICS_DUP.out.versions.first())
    ch_coverage_con  = COLLECTHSMETRICS_DUP.out.coverage
    ch_hsmetrics_con = COLLECTHSMETRICS_DUP.out.hsmetrics

    //
    // MODULE: Run Picard's Collect HS Metrics for consensus BAM files
    //
    COLLECTHSMETRICS_SIM(ch_bam_sim_stix, ch_fasta, ch_fai, ch_dict, params.hsmetrics_baits, params.hsmetrics_trgts, params.seq_library)
    ch_versions = ch_versions.mix(COLLECTHSMETRICS_SIM.out.versions.first())
    ch_coverage_con  = COLLECTHSMETRICS_SIM.out.coverage
    ch_hsmetrics_con = COLLECTHSMETRICS_SIM.out.hsmetrics

    //
    // MODULE: Run MSI Sensor PRO
    ///
    MSISENSORPRO_CON(ch_bam_con_stix, ch_msi_f)
    ch_versions = ch_versions.mix(MSISENSORPRO_CON.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(MSISENSORPRO_CON.out.summary.map{it[1]}.collect())
    ch_multiqc_files = ch_multiqc_files.mix(MSISENSORPRO_CON.out.msi_uns.map{it[1]}.collect())
    ch_multiqc_files = ch_multiqc_files.mix(MSISENSORPRO_CON.out.msi_all.map{it[1]}.collect())

    //
    // MODULE: Run MSI Sensor PRO
    ///
    MSISENSORPRO_DUP(ch_bam_dup_stix, ch_msi_f)
    ch_versions = ch_versions.mix(MSISENSORPRO_DUP.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(MSISENSORPRO_DUP.out.summary.map{it[1]}.collect())
    ch_multiqc_files = ch_multiqc_files.mix(MSISENSORPRO_DUP.out.msi_uns.map{it[1]}.collect())
    ch_multiqc_files = ch_multiqc_files.mix(MSISENSORPRO_DUP.out.msi_all.map{it[1]}.collect())

    //
    // MODULE: Run MSI Sensor PRO
    ///
    MSISENSORPRO_SIM(ch_bam_sim_stix, ch_msi_f)
    ch_versions = ch_versions.mix(MSISENSORPRO_SIM.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(MSISENSORPRO_SIM.out.summary.map{it[1]}.collect())
    ch_multiqc_files = ch_multiqc_files.mix(MSISENSORPRO_SIM.out.msi_uns.map{it[1]}.collect())
    ch_multiqc_files = ch_multiqc_files.mix(MSISENSORPRO_SIM.out.msi_all.map{it[1]}.collect())

    //
    // MODULE: Extract FastQ reads from BAM
    //
    SAMTOOLS_COLLATEFASTQ(ch_bam_dup_stix, ch_fasta, [])
    ch_versions = ch_versions.mix(SAMTOOLS_COLLATEFASTQ.out.versions)
    ch_consensus_reads = SAMTOOLS_COLLATEFASTQ.out.fastq

    //
    // MODULE: Run FastQC
    //
    FASTQC_CONSENSUS(ch_consensus_reads)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_CONSENSUS.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC_CONSENSUS.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'software_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    emit:

    ubam            = ch_ubam
    raw_bam         = ch_raw_sort_bam
    raw_bai         = ch_raw_sort_bai
    raw_baix        = ch_bam_fcu_stix
    versions        = ch_collated_versions
    group_bam       = ch_bam_grouped
    duplex_bam      = ch_bam_dup_stix
    split_reads     = ch_split_reads
    multiqc_files   = ch_multiqc_files
    finalized_bam   = ch_bam_con_stix
    split_contigs   = ch_split_contigs
    reads_finalized = ch_consensus_reads

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                                THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
