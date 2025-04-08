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
include { GRIDSS                                                                                                                    } from '../../modules/local/gridss/main'
include { BWAMEM2                                                                                                                   } from '../../modules/local/bwamem2/main' 
include { RECALL_SV                                                                                                                 } from '../../modules/local/recallsv/main'
include { TIDDIT_SV                                                                                                                 } from '../../modules/local/tiddit/sv/main'
include { MANTA_SOMATIC                                                                                                             } from '../../modules/local/manta/somatic/main'
include { SURVIVOR_MERGE                                                                                                            } from '../../modules/local/survivor/merge/main'
include { ANNOTSV_ANNOTSV                                                                                                           } from '../../modules/local/annotsv/annotsv/main' 
include { MANTA_TUMORONLY                                                                                                           } from '../../modules/local/manta/tumoronly/main'
include { SURVIVOR_FILTER                                                                                                           } from '../../modules/local/survivor/filter/main'
include { TABIX_BGZIPTABIX                                                                                                          } from '../../modules/nf-core/tabix/bgziptabix/main'
include { STANDARDIZE_SVS_MERGE                                                                                                     } from '../../modules/local/standarize_svs/merge/main'
include { STANDARDIZE_SVS_FILTER                                                                                                    } from '../../modules/local/standarize_svs/filter/main'
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
    ch_normal_bam
    ch_normal_bai
    ch_known_sites
    ch_bcf_mpileup
    ch_split_reads
    ch_reads_finalized
    ch_intervals_gunzip
    ch_intervals_gunzip_index     

    main:
    ch_versions = Channel.empty()
//    ch_multiqc_files = Channel.empty()

    //
    // MODULE: BWA-MEM2 mapping
    //
    sort_bam = 'sort'
    BWAMEM2(ch_reads_finalized, ch_split_reads, ch_bwa2, ch_fasta, ch_fai, sort_bam)
    ch_versions = ch_versions.mix(BWAMEM2.out.versions)
    ch_bam = BWAMEM2.out.bam

    //
    // MODULE: Run Manta in Only Tumour mode
    //
    MANTA_SOMATIC(ch_bam, ch_intervals_gunzip_index, ch_intervals_gunzip, ch_fasta, ch_normal_bam, ch_normal_bai, ch_fai, [])
    ch_versions = ch_versions.mix(MANTA_SOMATIC.out.versions)
    ch_manta_vcf = MANTA_SOMATIC.out.vcf

    //
    // MODULE: Run Manta in Only Tumour mode
    //
    MANTA_TUMORONLY(ch_bam, ch_intervals_gunzip, ch_intervals_gunzip_index, ch_fasta, ch_fai, [])
    ch_versions = ch_versions.mix(MANTA_TUMORONLY.out.versions)
    ch_manta_candidate_small_indels_vcf = MANTA_TUMORONLY.out.candidate_small_indels_vcf
    ch_manta_candidate_small_indels_vcf_tbi = MANTA_TUMORONLY.out.candidate_small_indels_vcf_tbi

    //
    // MODULE: Run TIDDIT in SV mode
    //
    TIDDIT_SV(ch_bam, ch_fasta, ch_bwa)
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
    DELLY(ch_bam, ch_fasta, ch_fai, params.exclude_bed, params.normal_bam, params.normal_bai)
    ch_versions = ch_versions.mix(DELLY.out.versions)
    ch_delly_vcf = DELLY.out.vcf

    //
    // MODULE: Run SvABA Note: version 1.2.0
    //
    SVABA(ch_bam, params.bwa, params.normal_bam, params.normal_bai, params.known_sites, params.known_sites_tbi, params.intervals)
    ch_versions = ch_versions.mix(SVABA.out.versions)
    ch_svaba_vcf = SVABA.out.vcf

    //
    // MODULE: Run Gridds (Extract overlapping fragments & calling)
    //
    GRIDSS(ch_bam, ch_fasta, ch_fai, params.normal_bam, params.normal_bai, params.intervals, params.blocklist_bed, params.bwa, params.kraken2db)
    ch_versions = ch_versions.mix(GRIDSS.out.versions)
    ch_gridss_vcf = GRIDSS.out.vcf

    //
    // MODULE: Run in-house module to standarize VCFs before merging step
    //
    STANDARDIZE_SVS_MERGE(ch_delly_vcf, ch_svaba_vcf, ch_manta_vcf, ch_tiddit_vcf, ch_gridss_vcf)
    ch_versions = ch_versions.mix(STANDARDIZE_SVS_MERGE.out.versions)
	ch_delly_std_vcf  = STANDARDIZE_SVS_MERGE.out.delly_vcf
	ch_svaba_std_vcf  = STANDARDIZE_SVS_MERGE.out.svaba_vcf
	ch_manta_std_vcf  = STANDARDIZE_SVS_MERGE.out.manta_vcf
	ch_tiddit_std_vcf = STANDARDIZE_SVS_MERGE.out.tiddit_vcf
	ch_gridss_std_vcf = STANDARDIZE_SVS_MERGE.out.gridss_vcf

    //
    // MODULE: Run Survivor to merge Unfiltered VCFs
    //
    SURVIVOR_MERGE(ch_delly_std_vcf, ch_svaba_std_vcf, ch_manta_std_vcf, ch_tiddit_std_vcf, ch_gridss_std_vcf, params.dict, 250, 0, 0, 1, 1000, 1)
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
    // MODULE: Run Gridds in ReCall mode
    //
    RECALL_SV(ch_bam, ch_fasta, ch_fai, ch_merged_int_list, ch_known_sites, params.refflat, params.normal_bam, params.normal_bai, params.intervals, params.blocklist_bed, params.bwa, params.kraken2db, params.pon_directory)
    ch_versions = ch_versions.mix(RECALL_SV.out.versions)
    ch_recall_vcf = RECALL_SV.out.vcf

    //
    // MODULE: Run in-house module to standarize VCFs before filtering step
    //
    STANDARDIZE_SVS_FILTER(ch_recall_vcf)
    ch_versions = ch_versions.mix(STANDARDIZE_SVS_FILTER.out.versions)
	ch_recall_std_vcf  = STANDARDIZE_SVS_FILTER.out.recall_vcf

    //
    // MODULE: Run Survivor to filter Unfiltered VCFs
    //
    SURVIVOR_FILTER(ch_delly_std_vcf, ch_svaba_std_vcf, ch_manta_std_vcf, ch_tiddit_std_vcf, ch_gridss_std_vcf, ch_recall_std_vcf, 1000, 2, 0, 0, 0, 1000)
    ch_versions = ch_versions.mix(SURVIVOR_FILTER.out.versions)
    ch_filtered_vcf = SURVIVOR_FILTER.out.filtered_vcf

//    // 
//    // MODULE: Run AnnotSV InstallAnnotations mode. To be run only once!
//    //
//    ANNOTSV_INSTALLANNOTATIONS()
//    ch_versions = ch_versions.mix(ANNOTSV_INSTALLANNOTATIONS.out.versions)

//    //
//    // MODULE: Run in-house module to standarize VCFs before filtering step
//    //
//    STANDARDIZE_SVS_ANNOTE(ch_filtered_vcf)
//    ch_versions = ch_versions.mix(STANDARDIZE_SVS_ANNOTE.out.versions)
//	ch_survivor_filt_std_vcf  = STANDARDIZE_SVS_ANNOTE.out.vcf

    //
    // MODULE: Run AnnotSV
    //
    ANNOTSV_ANNOTSV(ch_filtered_vcf, ch_bcf_mpileup, params.chrgtf, params.allow_list_genes, params.genome, params.annotsv_dir)
    ch_versions = ch_versions.mix(ANNOTSV_ANNOTSV.out.versions)
    ch_sv_annotated = ANNOTSV_ANNOTSV.out.tsv
    
    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'software_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    emit:

    versions        = ch_collated_versions
    sv_annotated    = ch_sv_annotated
//    draw_sv_pdf     = ch_draw_sv_pdf
//    multiqc_files   = ch_multiqc_files 

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                                THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

