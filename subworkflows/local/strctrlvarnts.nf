/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                       STRCTRLVARNTS SUBWORKFLOW                                                    
*******************************************************************************************************************************************************************************************************

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                             IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DELLY                                                                                                                     } from '../../modules/local/delly/main'
include { SVABA                                                                                                                     } from '../../modules/local/svaba/main'
include { DRAWSV                                                                                                                    } from '../../modules/local/drawsv/main'
include { GRIDSS                                                                                                                    } from '../../modules/local/gridss/main'
include { SVSOMF                                                                                                                    } from '../../modules/local/svsomf/main'
include { BWAMEM2                                                                                                                   } from '../../modules/local/bwamem2/main' 
include { RECALL_SV                                                                                                                 } from '../../modules/local/recallsv/main'
include { TIDDIT_SV                                                                                                                 } from '../../modules/local/tiddit/sv/main'
include { IANNOTATESV                                                                                                               } from '../../modules/local/iannotatesv/main'
include { MANTA_SOMATIC                                                                                                             } from '../../modules/local/manta/somatic/main'
include { SURVIVOR_MERGE                                                                                                            } from '../../modules/local/survivor/merge/main'
include { SURVIVOR_STATS                                                                                                            } from '../../modules/local/survivor/stats/main'
include { ANNOTSV_ANNOTSV                                                                                                           } from '../../modules/local/annotsv/annotsv/main' 
include { MANTA_TUMORONLY                                                                                                           } from '../../modules/local/manta/tumoronly/main'
include { SURVIVOR_FILTER                                                                                                           } from '../../modules/local/survivor/filter/main'
include { SERACARE_CHECKUP                                                                                                          } from '../../modules/local/seracare/checkup/main'
include { TABIX_BGZIPTABIX                                                                                                          } from '../../modules/nf-core/tabix/bgziptabix/main'
include { GATK4_BEDTOINTERVALLIST                                                                                                   } from '../../modules/nf-core/gatk4/bedtointervallist/main'
include { ANNOTSV_INSTALLANNOTATIONS                                                                                                } from '../../modules/nf-core/annotsv/installannotations/main'

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

workflow STRCTRLVARNTS {

    take:
    ch_bwa
    ch_fai
    ch_bwa2
    ch_dict
    ch_fasta
    ch_bam_pairs
    ch_known_sites
    ch_bcf_mpileup
    ch_split_reads
    ch_reads_finalized
    ch_intervals_gunzip
    ch_intervals_gunzip_index

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: Run Manta in Only Tumour mode
    //
    MANTA_SOMATIC(ch_bam_pairs, ch_intervals_gunzip_index, ch_intervals_gunzip, ch_fasta, ch_fai, [])
    ch_versions = ch_versions.mix(MANTA_SOMATIC.out.versions)
    ch_manta_vcf = MANTA_SOMATIC.out.vcf

    //
    // MODULE: Run Manta in Only Tumour mode
    //
    MANTA_TUMORONLY(ch_bam_pairs, ch_intervals_gunzip, ch_intervals_gunzip_index, ch_fasta, ch_fai, [])
    ch_versions = ch_versions.mix(MANTA_TUMORONLY.out.versions)
    ch_manta_candidate_small_indels_vcf = MANTA_TUMORONLY.out.candidate_small_indels_vcf
    ch_manta_candidate_small_indels_vcf_tbi = MANTA_TUMORONLY.out.candidate_small_indels_vcf_tbi

    //
    // MODULE: Run TIDDIT in SV mode
    //
    TIDDIT_SV(ch_bam_pairs, ch_fasta, ch_bwa)
    ch_versions = ch_versions.mix(TIDDIT_SV.out.versions)
    ch_tiddit_vcf = TIDDIT_SV.out.vcf
    ch_tiddit_ploidy = TIDDIT_SV.out.ploidy

    //
    // MODULE: Run BGZIP & Tabix
    //
    TABIX_BGZIPTABIX(ch_tiddit_vcf)
    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions)
    ch_tiddit_vcf_zip = TABIX_BGZIPTABIX.out.gz_tbi

    //
    // MODULE: Run Delly Call
    //
    DELLY(ch_bam_pairs, ch_fasta, ch_fai, params.exclude_bed)
    ch_versions = ch_versions.mix(DELLY.out.versions)
    ch_delly_vcf = DELLY.out.vcf

    //
    // MODULE: Run SvABA Note: version 1.2.0
    //
    SVABA(ch_bam_pairs, params.bwa, params.known_sites, params.known_sites_tbi, params.intervals)
    ch_versions = ch_versions.mix(SVABA.out.versions)
    ch_svaba_vcf = SVABA.out.vcf

    //
    // MODULE: Run Gridds (Extract overlapping fragments & calling)
    //
    GRIDSS(ch_bam_pairs, ch_fasta, ch_fai, params.intervals, params.blocklist_bed, params.bwa, params.kraken2db)
    ch_versions = ch_versions.mix(GRIDSS.out.versions)
    ch_gridss_vcf = GRIDSS.out.vcf

    //
    // Combine the vcf by meta key patient
    //
    ch_survivor_merge_input = ch_delly_vcf
        .map { meta, vcf -> [meta.patient, meta, vcf] }
        .join(ch_gridss_vcf.map { meta, vcf -> [meta.patient, meta, vcf] })
        .join(ch_manta_vcf.map  { meta, vcf -> [meta.patient, meta, vcf] })
        .join(ch_svaba_vcf.map  { meta, vcf -> [meta.patient, meta, vcf] })
        .join(ch_tiddit_vcf.map { meta, vcf -> [meta.patient, meta, vcf] })
        .map { patient, meta_delly, delly_vcf, meta_gridss, gridss_vcf, meta_manta, manta_vcf, meta_svaba, svaba_vcf, meta_tiddit, tiddit_vcf ->
            tuple(
                meta_delly,                 //
                meta_delly,     delly_vcf,  //
                meta_gridss,    gridss_vcf, //
                meta_manta,     manta_vcf,  //  
                meta_svaba,     svaba_vcf,  //
                meta_tiddit,    tiddit_vcf  
            )
        }

    //
    // MODULE: Run Survivor to merge Unfiltered VCFs
    //
    SURVIVOR_MERGE(ch_survivor_merge_input, params.chromosomes, 1000, 2, 0, 0, 0, 30)
    ch_versions = ch_versions.mix(SURVIVOR_MERGE.out.versions)
    ch_merged_bed = SURVIVOR_MERGE.out.bed
    ch_merged_vcf = SURVIVOR_MERGE.out.vcf

    //
    // MODULE: Run GATK4 to Convert BED into Interval List
    //
    GATK4_BEDTOINTERVALLIST(ch_merged_bed, ch_dict)
    ch_versions = ch_versions.mix(GATK4_BEDTOINTERVALLIST.out.versions)
    ch_merged_int_list = GATK4_BEDTOINTERVALLIST.out.interval_list

    //
    // Join interval lists with BAM pairs based on patient
    //
    ch_recall_input = ch_bam_pairs
        .map {meta, tbam, tbai, nbam, nbai -> [meta.patient, meta, tbam, tbai, nbam, nbai] }
        .join(ch_merged_int_list.map { meta, interval_list -> [meta.patient, meta, interval_list] })
        .map { patient, meta_bam_pairs, tbam, tbai, nbam, nbai, meta_interval_list, interval_list ->
            tuple(
                meta_bam_pairs, tbam, tbai, nbam, nbai, interval_list
            )
        }

    //
    // MODULE: Run Gridds in ReCall mode
    //
    RECALL_SV(ch_recall_input, ch_fasta, ch_fai, ch_known_sites, params.refflat, params.intervals, params.blocklist_bed, params.bwa, params.kraken2db, params.pon_directory)
    ch_versions = ch_versions.mix(RECALL_SV.out.versions)
    ch_prefilter_vcf = RECALL_SV.out.vcf

    //
    // MODULE: Run SV Somatic Filter
    //
    SVSOMF(ch_prefilter_vcf, ch_fasta, ch_fai, params.pon_directory)
    ch_versions = ch_versions.mix(SVSOMF.out.versions)
    ch_recall_vcf = SVSOMF.out.vcf

    //
    // Combine the vcf by meta key patient
    //
    ch_survivor_filter_input = ch_delly_vcf
        .map { meta, vcf -> [meta.patient, meta, vcf] }
        .join(ch_gridss_vcf.map { meta, vcf -> [meta.patient, meta, vcf] })
        .join(ch_manta_vcf.map  { meta, vcf -> [meta.patient, meta, vcf] })
        .join(ch_recall_vcf.map { meta, vcf -> [meta.patient, meta, vcf] })
        .join(ch_svaba_vcf.map  { meta, vcf -> [meta.patient, meta, vcf] })
        .join(ch_tiddit_vcf.map { meta, vcf -> [meta.patient, meta, vcf] })
        .map { patient, meta_delly, delly_vcf, meta_gridss, gridss_vcf, meta_manta, manta_vcf, meta_recall, recall_vcf, meta_svaba, svaba_vcf, meta_tiddit, tiddit_vcf ->
            tuple(
                meta_delly,                 //
                meta_delly,     delly_vcf,  //
                meta_gridss,    gridss_vcf, //
                meta_manta,     manta_vcf,  //
                meta_recall,    recall_vcf, //
                meta_svaba,     svaba_vcf,  //
                meta_tiddit,    tiddit_vcf  
            )
        }

    //
    // MODULE: Run Survivor to filter Unfiltered VCFs
    //
    SURVIVOR_FILTER(ch_survivor_filter_input, 10000, 3, 0, 0, 0, 1000)
    ch_versions = ch_versions.mix(SURVIVOR_FILTER.out.versions)
    ch_filtered_all = SURVIVOR_FILTER.out.filtered_all
    ch_filtered_vcf = SURVIVOR_FILTER.out.filtered_vcf

    //
    // MODULE: Run Survivor Stats
    //
    SURVIVOR_STATS(ch_filtered_vcf, -1, -1, -1)
    ch_versions = ch_versions.mix(SURVIVOR_STATS.out.versions)
    ch_multiqc_files  = ch_multiqc_files.mix(SURVIVOR_STATS.out.stats.collect{it[1]}.ifEmpty([]))

    //
    // MODULE: Run AnnotSV
    //
    ANNOTSV_ANNOTSV(ch_filtered_vcf, ch_bcf_mpileup, params.chrgtf, params.allow_list_genes, params.genome, params.annotsv_dir)
    ch_versions = ch_versions.mix(ANNOTSV_ANNOTSV.out.versions)
    ch_sv_annotated = ANNOTSV_ANNOTSV.out.tsv

    //
    // MODULE: Run iAnnotateSV 
    //
    IANNOTATESV(ch_filtered_all)
    ch_versions = ch_versions.mix(IANNOTATESV.out.versions)
    ch_annotated_tsv = IANNOTATESV.out.tsv
    ch_annotated_ann = IANNOTATESV.out.ann

    //
    // Join annotated SVs with BAM pairs based on patient
    //
    ch_bam_pairs_by_patient = ch_bam_pairs.map {meta, bam_t, bai_t, bam_n, bai_n -> tuple(meta.patient, meta, bam_t, bai_t, bam_n, bai_n) }
    ch_drawsv_input = ch_bam_pairs_by_patient
        .join(ch_annotated_tsv
        .map {meta, tsv -> tuple(meta.patient, meta, tsv)}
        )
        .map { patient, meta_b, bam_t, bai_t, bam_n, bai_n, meta_t, tsv ->
            tuple(
                meta_b, 
                meta_b, bam_t, bai_t, 
                meta_b, bam_n, bai_n, 
                meta_t, tsv
            )
        }

    //
    // MODULE: Run DrawSV
    //
    DRAWSV(ch_drawsv_input, params.annotations, params.cytobands, params.drawsv_chr, params.protein_domains)
    ch_versions = ch_versions.mix(DRAWSV.out.versions)
    ch_drawsv_pdf = DRAWSV.out.pdf

    //
    // Check-Up for SeraCare samples only
    //
    ch_seracare_sample = ch_annotated_tsv
        .filter { meta, file -> 
            meta.id.contains("SeraCare") 
        }

    //
    // MODULE: Run SeraCare Check-Up
    //
    SERACARE_CHECKUP(ch_seracare_sample)
    ch_versions = ch_versions.mix(SERACARE_CHECKUP.out.versions)
  
    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info", 
            name: 'software_versions.yml', 
            sort: true, 
            newLine: true)
        .set { ch_collated_versions }

    emit:

    versions        = ch_collated_versions
    multiqc_files   = ch_multiqc_files 

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                                THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
