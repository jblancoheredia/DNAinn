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
include { SNPEFF_SNPEFF                                                                                                             } from '../../modules/nf-core/snpeff/snpeff/main'
include { BCFTOOLS_MERGE                                                                                                            } from '../../modules/local/bcftools/merge/main' 
include { ENSEMBLVEP_VEP                                                                                                            } from '../../modules/nf-core/ensemblvep/vep/main'
include { GATK4_HAPLOTYPECALLER                                                                                                     } from '../../modules/local/gatk4/haplotypecaller/main'
include { GATK4_FILTERMUTECTCALLS                                                                                                   } from '../../modules/local/gatk4/filtermutectcalls/main'   
include { GATK4_GETPILEUPSUMMARIES                                                                                                  } from '../../modules/local/gatk4/getpileupsummaries/main'
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
    ch_bam_pairs
    ch_known_sites
    ch_consensus_bam
    ch_gatk_interval_list

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: Run Haplotypecaller
    //
    GATK4_HAPLOTYPECALLER(ch_bam_pairs, ch_known_sites, ch_gatk_interval_list, ch_fasta, ch_dict, ch_fai)
    ch_versions  = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions.ifEmpty(null))
    ch_haplotypecaller_raw = GATK4_HAPLOTYPECALLER.out.vcf
    ch_haplotypecaller_tbi = GATK4_HAPLOTYPECALLER.out.tbi
    ch_haplotypecaller_raw_tbi_combined = ch_haplotypecaller_raw.join(ch_haplotypecaller_tbi)
    ch_haplotypecaller_vcf_tbi = ch_haplotypecaller_raw_tbi_combined.map { meta, vcf, tbi ->
        [meta, vcf, tbi]
    }

    //
    // MODULE: Run FgBio ClipBAM 
    //
    FGBIO_CLIPBAM(ch_bam_pairs, ch_fasta, ch_fai)
    ch_versions = ch_versions.mix(FGBIO_CLIPBAM.out.versions)
    ch_bam_clipped = FGBIO_CLIPBAM.out.bam
    ch_txt = FGBIO_CLIPBAM.out.txt
// TO-DO:    multi_qc_files = FGBIO_CLIPBAM.out.metrics

    //
    // MODULE: Run Mutect2
    //
    GATK4_MUTECT2(ch_bam_clipped, ch_intervals, ch_fasta, ch_fai, ch_dict, ch_txt, params.germres, params.germres_tbi, params.pon, params.pon_tbi)
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
    ch_filtermutect_joined = ch_mutect2_vcf
        .join(ch_mutect2_stats, by: 0)
        .join(ch_orientation  , by: 0)
        .join(ch_segmentation , by: 0)
        .join(ch_contamination, by: 0)
    ch_filtermutect_vcf           = ch_filtermutect_joined.map { meta, vcf, tbi, stats, orientation, segmentation, contamination -> tuple(meta, vcf, tbi)      }
    ch_filtermutect_stats         = ch_filtermutect_joined.map { meta, vcf, tbi, stats, orientation, segmentation, contamination -> tuple(meta, stats)         }
    ch_filtermutect_orientation   = ch_filtermutect_joined.map { meta, vcf, tbi, stats, orientation, segmentation, contamination -> tuple(meta, orientation)   }
    ch_filtermutect_segmentation  = ch_filtermutect_joined.map { meta, vcf, tbi, stats, orientation, segmentation, contamination -> tuple(meta, segmentation)  }
    ch_filtermutect_contamination = ch_filtermutect_joined.map { meta, vcf, tbi, stats, orientation, segmentation, contamination -> tuple(meta, contamination) }
    GATK4_FILTERMUTECTCALLS(ch_filtermutect_vcf, ch_fasta, ch_fai, ch_dict, ch_filtermutect_stats, ch_filtermutect_orientation, ch_filtermutect_segmentation, ch_filtermutect_contamination, [])
    ch_versions = ch_versions.mix(GATK4_FILTERMUTECTCALLS.out.versions)
    ch_mutect2_filtered_vcf = GATK4_FILTERMUTECTCALLS.out.vcf
    ch_mutect2_filtered_stats = GATK4_FILTERMUTECTCALLS.out.stats

    //
    // MODULE: Run VarDictJava
    //
    VARDICTJAVA(ch_bam_clipped, ch_fasta, ch_fai, params.vardict_bed)
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
    LOFREQ(ch_bam_clipped, ch_intervals, ch_fasta, ch_fai, params.known_sites_tbi, params.known_sites)
    ch_versions = ch_versions.mix(LOFREQ.out.versions)
    ch_lofreq_snvs_vcf = LOFREQ.out.vcf_snvs
    ch_lofreq_snvs_indels = LOFREQ.out.vcf_indels

    //
    // Build per-caller channel and annotate with SnpEff + VEP
    //
    ch_callers_vcf = ch_lofreq_snvs_vcf.map      { meta, vcf, idx -> tuple(meta + [caller: 'lofreq.somatic_final.snvs'], vcf, idx) }
        .mix(ch_lofreq_snvs_indels.map           { meta, vcf, idx -> tuple(meta + [caller: 'lofreq.somatic_final.indels'], vcf, idx) })
        .mix(ch_mutect2_filtered_vcf.map         { meta, vcf, idx -> tuple(meta + [caller: 'mutect2.filtered'], vcf, idx) })
        .mix(ch_freebayes_vcf.map                { meta, vcf, idx -> tuple(meta + [caller: 'freebayes'], vcf, idx) })
        .mix(ch_vardict_vcf.map                  { meta, vcf, idx -> tuple(meta + [caller: 'vardict'], vcf, idx) })

    ch_snpeff_in = ch_callers_vcf.map { meta, vcf, idx -> tuple(meta, vcf) }
    ch_snpeff_cache = ch_callers_vcf.map { meta, vcf, idx -> tuple(meta, file(params.snpeff_cache)) }
    SNPEFF_SNPEFF(ch_snpeff_in, params.snpeff_db, ch_snpeff_cache)
    ch_multiqc_files = ch_multiqc_files.mix(SNPEFF_SNPEFF.out.report.map{ it[3] }.collect())
    ch_multiqc_files = ch_multiqc_files.mix(SNPEFF_SNPEFF.out.summary_html.map{ it[3] }.collect())
    ch_multiqc_files = ch_multiqc_files.mix(SNPEFF_SNPEFF.out.genes_txt.map{ it[3] }.collect())

    ch_vep_in = SNPEFF_SNPEFF.out.vcf.map { meta, vcf ->
        tuple(meta, vcf, [])
    }
    ch_vep_cache = SNPEFF_SNPEFF.out.vcf.map { meta, vcf -> tuple(meta, file(params.vep_cache)) }
    ch_vep_fasta = SNPEFF_SNPEFF.out.vcf.map { meta, vcf -> tuple(meta, file(params.fasta)) }
    ENSEMBLVEP_VEP(ch_vep_in, params.genome, params.vep_species, params.ensembl_version, ch_vep_cache, ch_vep_fasta, [])
    ch_multiqc_files = ch_multiqc_files.mix(ENSEMBLVEP_VEP.out.report.map{ it[3] }.collect())

    //
    // Pre-Merge: collect all callers and reshape to proper tuple
    //
    def caller_order = ['freebayes', 'lofreq.somatic_final.snvs', 'lofreq.somatic_final.indels', 'mutect2.filtered', 'vardict']
    ch_all_vcfs = ENSEMBLVEP_VEP.out.vcf
        .join(ENSEMBLVEP_VEP.out.tbi, by: 0)
        .map { meta, vcf, tbi -> [meta.id, [meta: [id: meta.id], caller: meta.caller, vcf: vcf, idx: tbi]] }
        .groupTuple()
        .map { id, calls ->
            def sorted_calls = calls.sort { caller_order.indexOf(it.caller) }
            def meta = sorted_calls[0].meta
            def vcfs = sorted_calls.collect { it.vcf }
            def tbis = sorted_calls.collect { it.idx }
            tuple(meta, vcfs, tbis)
        }
        .set { ch_pre_merge }

    //
    // MODULE: From the Callers VCF files produce a single TSV 
    //
    VCFCALLS2TSV(ch_pre_merge)
    ch_versions = ch_versions.mix(VCFCALLS2TSV.out.versions)
    ch_variants_tsv = VCFCALLS2TSV.out.tsv

//    //
//    // MODULES: Run GetBaseCountsMultiSample
//    //
//    GETBASECOUNTS_MULTISAMPLE(ch_pre_merge, ch_bam_clipped, ch_fasta, ch_fai)
//    ch_versions = ch_versions.mix(GETBASECOUNTS_MULTISAMPLE.out.versions)

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
