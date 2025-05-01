/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                       VARIANTDSCVRY SUBWORKFLOW                                                    
*******************************************************************************************************************************************************************************************************

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                             IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { LOFREQ                                                                                                                    } from '../../modules/local/lofreq/main'
include { FREEBAYES                                                                                                                 } from '../../modules/local/freebayes/main'
include { VARDICTJAVA                                                                                                               } from '../../modules/local/vardictjava/main'  
include { VCFCALLS2TSV                                                                                                              } from '../../modules/local/vcfcalls2tsv/main'
include { FGBIO_CLIPBAM                                                                                                             } from '../../modules/local/fgbio/clipbam/main'
include { GATK4_MUTECT2                                                                                                             } from '../../modules/local/gatk4/mutect2/main' 
include { BCFTOOLS_MERGE                                                                                                            } from '../../modules/local/bcftools/merge/main' 
include { GATK4_FILTERMUTECTCALLS                                                                                                   } from '../../modules/local/gatk4/filtermutectcalls/main'   
include { GATK4_GETPILEUPSUMMARIES                                                                                                  } from '../../modules/nf-core/gatk4/getpileupsummaries/main'
include { GETBASECOUNTS_MULTISAMPLE                                                                                                 } from '../../modules/local/getbasecounts/multisample/main'
include { GATK4_CALCULATECONTAMINATION                                                                                              } from '../../modules/nf-core/gatk4/calculatecontamination/main' 
include { GATK4_LEARNREADORIENTATIONMODEL                                                                                           } from '../../modules/nf-core/gatk4/learnreadorientationmodel/main'

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

workflow VARIANTDSCVRY {

    take:
    ch_fai
    ch_dict
    ch_fasta
    ch_targets
    ch_intervals
    ch_consensus_bam

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: Run FgBio ClipBAM 
    //
    FGBIO_CLIPBAM(ch_consensus_bam, ch_fasta, ch_fai)
    ch_versions = ch_versions.mix(FGBIO_CLIPBAM.out.versions)
    ch_bam_clipped = FGBIO_CLIPBAM.out.bam
// TO-DO:    multi_qc_files = FGBIO_CLIPBAM.out.metrics

    //
    // MODULE: Run Mutect2
    //
    GATK4_MUTECT2(ch_bam_clipped, ch_intervals, ch_fasta, ch_fai, ch_dict, params.germres, params.germres_tbi, params.pon, params.pon_tbi, params.normal_bam, params.normal_bai)
    ch_versions = ch_versions.mix(GATK4_MUTECT2.out.versions)
    ch_mutect2_vcf = GATK4_MUTECT2.out.vcf
    ch_mutect2_f1r2 = GATK4_MUTECT2.out.f1r2
    ch_mutect2_stats = GATK4_MUTECT2.out.stats

    //
    // MODULE: Run GATK4 Learn Orientation Model
    //
    GATK4_LEARNREADORIENTATIONMODEL(ch_mutect2_f1r2)
    ch_versions = ch_versions.mix(GATK4_LEARNREADORIENTATIONMODEL.out.versions)
    ch_orientation = GATK4_LEARNREADORIENTATIONMODEL.out.artifactprior

    //
    // MODULE: Run GATK4 Get PileUp Summaries
    //
    GATK4_GETPILEUPSUMMARIES(ch_bam_clipped, ch_intervals, ch_fasta, ch_fai, ch_dict, params.germres, params.germres_tbi)
    ch_versions = ch_versions.mix(GATK4_GETPILEUPSUMMARIES.out.versions)
    ch_pileups = GATK4_GETPILEUPSUMMARIES.out.table

    //
    // MODULE: Run GATK4 Calculate Contamination
    //
    GATK4_CALCULATECONTAMINATION(ch_pileups)
    ch_versions = ch_versions.mix(GATK4_CALCULATECONTAMINATION.out.versions)
    ch_segmentation = GATK4_CALCULATECONTAMINATION.out.segmentation
    ch_contamination = GATK4_CALCULATECONTAMINATION.out.contamination

    //
    // MODULE: Run Filter MUTECT2 Calls
    //
    GATK4_FILTERMUTECTCALLS(ch_mutect2_vcf, ch_fasta, ch_fai, ch_dict, ch_mutect2_stats, ch_orientation, ch_segmentation, ch_contamination, [])
    ch_versions = ch_versions.mix(GATK4_FILTERMUTECTCALLS.out.versions)
    ch_mutect2_filtered_vcf = GATK4_FILTERMUTECTCALLS.out.vcf
    ch_mutect2_filtered_stats = GATK4_FILTERMUTECTCALLS.out.stats

    //
    // MODULE: Run VarDictJava
    //
    VARDICTJAVA(ch_bam_clipped, ch_fasta, ch_fai, params.vardict_bed, params.normal_bam, params.normal_bai)
    ch_versions = ch_versions.mix(VARDICTJAVA.out.versions)
    ch_vardict_vcf = VARDICTJAVA.out.vcf

    //
    // MODULE: Run FreeBayes
    //
    FREEBAYES(ch_bam_clipped, ch_fasta, ch_fai, params.targets_bed)
    ch_versions = ch_versions.mix(FREEBAYES.out.versions)
    ch_freebayes_vcf = FREEBAYES.out.vcf

    //
    // MODULE: Run LoFreq to predict variants using multiple processors
    //
    LOFREQ(ch_bam_clipped, ch_intervals, ch_fasta, ch_fai, params.known_sites_tbi, params.known_sites, params.normal_bam, params.normal_bai)
    ch_versions = ch_versions.mix(LOFREQ.out.versions)
    ch_lofreq_snvs_vcf = LOFREQ.out.vcf_snvs
    ch_lofreq_snvs_indels = LOFREQ.out.vcf_indels

    //
    // Pre-Merge: collect all callers and reshape to proper tuple
    //
    ch_all_vcfs = ch_lofreq_snvs_vcf.map { meta, vcf, idx -> [meta.id, [meta: meta, caller: 'lofreq', vcf: vcf, idx: idx]] }
        .mix(ch_lofreq_snvs_indels.map { meta, vcf, idx -> [meta.id, [meta: meta, caller: 'lofreq', vcf: vcf, idx: idx]] })
        .mix(ch_mutect2_filtered_vcf.map { meta, vcf, idx -> [meta.id, [meta: meta, caller: 'mutect2', vcf: vcf, idx: idx]] })
        .mix(ch_freebayes_vcf.map { meta, vcf, idx -> [meta.id, [meta: meta, caller: 'freebayes', vcf: vcf, idx: idx]] })
        .mix(ch_vardict_vcf.map { meta, vcf, idx -> [meta.id, [meta: meta, caller: 'vardict', vcf: vcf, idx: idx]] })
        .groupTuple()
        .map { id, calls ->
            def meta = calls[0].meta
            def vcfs = calls.collect { it.vcf }
            def tbis = calls.collect { it.idx }
            tuple(meta, vcfs, tbis)
        }
        .set { ch_pre_merge }

    //
    // MODULE: From the Callers VCF files produce a single TSV 
    //
    VCFCALLS2TSV(ch_pre_merge)
    ch_versions = ch_versions.mix(VCFCALLS2TSV.out.versions)
    ch_variants_tsv = VCFCALLS2TSV.out.tsv

    //
    // MODULE: Run BCFtools to merge the VCFs
    //
    BCFTOOLS_MERGE(ch_pre_merge, ch_fasta, ch_fai, params.targets_bed)
    ch_versions = ch_versions.mix(BCFTOOLS_MERGE.out.versions)
    ch_merged_vcf = BCFTOOLS_MERGE.out.vcf


    //
    // MODULES: Run GetBaseCountsMultiSample
    //
    GETBASECOUNTS_MULTISAMPLE(ch_bam_clipped, ch_fasta, ch_fai, ch_merged_vcf)
    ch_versions = ch_versions.mix(GETBASECOUNTS_MULTISAMPLE.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'software_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    emit:

    versions        = ch_collated_versions
    variants        = ch_variants_tsv
//    multiqc_files   = ch_multiqc_files 

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                                THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
