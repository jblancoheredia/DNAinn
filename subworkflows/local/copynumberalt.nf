/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                       COPYNUMBERALT SUBWORKFLOW                                                    
*******************************************************************************************************************************************************************************************************

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                             IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ONCOCNV                                                                                                                   } from '../../modules/nf-core/oncocnv/main'
include { COPYNCAT                                                                                                                  } from '../../modules/local/copyncat/main'
include { FACETS_CNV                                                                                                                } from '../../modules/local/facets/main'
include { CNVKIT_CALL                                                                                                               } from '../../modules/nf-core/cnvkit/call/main'
include { CNVKIT_BATCH                                                                                                              } from '../../modules/local/cnvkit/batch/main'
include { CNVKIT_EXPORT                                                                                                             } from '../../modules/nf-core/cnvkit/export/main'
include { SEQUENZA_FITS                                                                                                             } from '../../modules/local/sequenza/fits/main'
include { BCFTOOLS_MPILEUP                                                                                                          } from '../../modules/local/bcftools/mpileup/main'
include { SAMTOOLS_MPILEUP                                                                                                          } from '../../modules/local/samtools/mpileup/main'
include { CNVKIT_REFERENCE                                                                                                          } from '../../modules/nf-core/cnvkit/reference/main'
include { CNVKIT_ANTITARGET                                                                                                         } from '../../modules/nf-core/cnvkit/antitarget/main'
include { CNVKIT_GENEMETRICS                                                                                                        } from '../../modules/nf-core/cnvkit/genemetrics/main'
include { CONTROLFREEC_OT_FREEC                                                                                                     } from '../../modules/local/controlfreec_ot/freec/main'
include { SEQUENZAUTILS_BAM2SEQZ                                                                                                    } from '../../modules/local/sequenzautils/bam2seqz/main'
include { CONTROLFREEC_OT_FREEC2BED                                                                                                 } from '../../modules/local/controlfreec_ot/freec2bed/main' 
include { CONTROLFREEC_OT_MAKEGRAPH2                                                                                                } from '../../modules/local/controlfreec_ot/makegraph2/main'
include { CONTROLFREEC_OT_FREEC2CIRCOS                                                                                              } from '../../modules/local/controlfreec_ot/freec2circos/main'
include { CONTROLFREEC_OT_ASSESSSIGNIFICANCE                                                                                        } from '../../modules/local/controlfreec_ot/assesssignificance/main'

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

workflow COPYNUMBERALT {

    take:
    ch_fai
    ch_fasta
    ch_raw_bam
    ch_raw_bai
    ch_bam_pairs
    ch_intervals
    ch_targets_bed
    ch_consensus_bam
    ch_normal_con_bam
    ch_normal_con_bai
    ch_cnvkit_reference
    ch_cnvkit_antitarget

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

//    //
//    // MODULE: Run OncoCNV
//    //
//    ONCOCNV(ch_bam_pairs, ch_fasta, ch_targets_bed, ch_fai)
//    ch_versions = ch_versions.mix(ONCOCNV.out.versions)
//    ch_oncocnv_png      = ONCOCNV.out.png
//    ch_oncocnv_profile  = ONCOCNV.out.profile
//    ch_oncocnv_summary  = ONCOCNV.out.summary

    //
    // MODULE: Run BCFtools Mpileup
    //
    BCFTOOLS_MPILEUP(ch_consensus_bam, ch_fai, ch_fasta, ch_intervals)
    ch_versions = ch_versions.mix(BCFTOOLS_MPILEUP.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(BCFTOOLS_MPILEUP.out.stats.collect{it[1]})
    ch_bcf_mpileup = BCFTOOLS_MPILEUP.out.vcf
    ch_bcf_mpileup_tbi = BCFTOOLS_MPILEUP.out.tbi

    //
    // MODULE: Run SamTools Mpileup
    //
    SAMTOOLS_MPILEUP(ch_consensus_bam, ch_intervals, params.fasta)
    ch_versions = ch_versions.mix(SAMTOOLS_MPILEUP.out.versions)
    ch_sam_mpileup = SAMTOOLS_MPILEUP.out.pup
    ch_sam_mpileup_tbi = SAMTOOLS_MPILEUP.out.tbi
  
//    //
//    // MODULE: Run ControlFreec Freec 
//    //
//    CONTROLFREEC_OT_FREEC(ch_sam_mpileup, ch_consensus_bam, params.fasta, params.chr_fai, params.known_sites, params.known_sites_tbi, params.by_chr_dir, params.cf_mappability, params.intervals, params.cf_control_mpileup, params.cf_coeff, params.cf_contamination, params.cf_contamination_adjustment, params.cf_ploidy)
//    ch_versions = ch_versions.mix(CONTROLFREEC_OT_FREEC.out.versions)
//    ch_cfot_baf   = CONTROLFREEC_OT_FREEC.out.BAF
//    ch_cfot_cnvs  = CONTROLFREEC_OT_FREEC.out.CNV
//    ch_cfot_ratio = CONTROLFREEC_OT_FREEC.out.ratio
//
//    //
//    // MODULE: Run ControlFreec Assess Significance
//    //
//    CONTROLFREEC_OT_ASSESSSIGNIFICANCE(ch_cfot_cnvs, ch_cfot_ratio)
//    ch_versions = ch_versions.mix(CONTROLFREEC_OT_ASSESSSIGNIFICANCE.out.versions)
//
//    //
//    // MODULE: Run ControlFreec Freec2BED
//    //
//    CONTROLFREEC_OT_FREEC2BED(ch_cfot_ratio)
//    ch_versions = ch_versions.mix(CONTROLFREEC_OT_FREEC2BED.out.versions)
//
//    //
//    // MODULE: Run ControlFreec Freec2Circos
//    //
//    CONTROLFREEC_OT_FREEC2CIRCOS(ch_cfot_ratio)
//    ch_versions = ch_versions.mix(CONTROLFREEC_OT_FREEC2CIRCOS.out.versions)
//
//    //
//    // MODULE: Run ControlFreec MakeGraph2
//    //
//    CONTROLFREEC_OT_MAKEGRAPH2(CONTROLFREEC_OT_FREEC.out.ratio.join(CONTROLFREEC_OT_FREEC.out.BAF, failOnDuplicate: true, failOnMismatch: true))
//    ch_versions = ch_versions.mix(CONTROLFREEC_OT_MAKEGRAPH2.out.versions)

//    //
//    // MODULE: Run CNVkit AntiTarget Module just once per panel
//    //
//    CNVKIT_ANTITARGET(ch_targets_bed)
//    ch_versions = ch_versions.mix(CNVKIT_ANTITARGET.out.versions)
//    ch_antitargets = CNVKIT_ANTITARGET.out.bed

//    //
//    // MODULE: Run CNVkit Reference Module just once per panel
//    //
//    CNVKIT_REFERENCE(ch_fasta, ch_targets_bed, ch_cnvkit_antitarget)
//    ch_versions = ch_versions.mix(CNVKIT_REFERENCE.out.versions)
//    ch_cnvkit_reference = CNVKIT_REFERENCE.out.cnn

//    //
//    // MODULE: Run CNVKIT Batch
//    //
//    generate_pon = false
//    CNVKIT_BATCH(ch_bam_pairs, ch_fasta, ch_fai, ch_cnvkit_antitarget, ch_cnvkit_reference, generate_pon)
//    ch_versions = ch_versions.mix(CNVKIT_BATCH.out.versions.first())
//    ch_cnvkit_call_input = CNVKIT_BATCH.out.cns.map{ meta, cns -> [meta, cns[2], []]}
//
//    //
//    // MODULE: Run CNVkit Call
//    //
//    CNVKIT_CALL(ch_cnvkit_call_input)
//    ch_versions = ch_versions.mix(CNVKIT_CALL.out.versions.first())
//    ch_cnvkit_call = CNVKIT_CALL.out.cns
//
//    //
//    // MODULE: Run CNVkit Export
//    //
//    CNVKIT_EXPORT(CNVKIT_CALL.out.cns)
//    ch_versions = ch_versions.mix(CNVKIT_EXPORT.out.versions.first())
//    ch_cnvkit_cns = CNVKIT_EXPORT.out.output
//
//    //
//    // MODULE: Run CNVkit GeneMetrics
//    //
//    ch_genemetrics_input = CNVKIT_BATCH.out.cnr.join(CNVKIT_BATCH.out.cns).map{ meta, cnr, cns -> [meta, cnr, cns[2]]}
//    CNVKIT_GENEMETRICS(ch_genemetrics_input)
//    ch_versions = ch_versions.mix(CNVKIT_GENEMETRICS.out.versions.first())
//    ch_multiqc_files = ch_multiqc_files.mix(CNVKIT_GENEMETRICS.out.tsv)
//
//    //
//    // MODULE: Run FACETS 
//    //
//    FACETS_CNV(ch_bam_pairs, ch_intervals, params.common_vcf, params.common_vcf_tbi)
//    ch_versions = ch_versions.mix(FACETS_CNV.out.versions.first())
//
//    //
//    // MODULE: Run Sequenzautils BAM2seqz
//    //
//    SEQUENZAUTILS_BAM2SEQZ(ch_bam_pairs, ch_fasta, ch_fai, params.wigfile20)
//    ch_seqz = SEQUENZAUTILS_BAM2SEQZ.out.seqz
//    ch_versions = ch_versions.mix(SEQUENZAUTILS_BAM2SEQZ.out.versions)
//
//    //
//    // MODULE: Run Sequenza Fits
//    //
//    SEQUENZA_FITS(ch_seqz)
//    ch_versions = ch_versions.mix(SEQUENZA_FITS.out.versions.first())

    // Build CopyNcat input
//    ch_copyncat_input = ch_cnvkit_call
//        .join(cnvkit_cns)
//        .join(CNVKIT_EXPORT.out.output).map { l, r -> tuple(l[0], l[1], l[2], r[1]) }
//        .join(CONTROLFREEC_OT_FREEC.out.CNV).map { l, r -> tuple(l[0], l[1], l[2], l[3], r[1]) }
//        .join(CONTROLFREEC_OT_FREEC.out.config).map { l, r -> tuple(l[0], l[1], l[2], l[3], l[4], r[1]) }
//        .join(SEQUENZA_FITS.out.segments).map { l, r -> tuple(l[0], l[1], l[2], l[3], l[4], l[5], r[1]) }
//        .join(SEQUENZA_FITS.out.confints_cp).map { l, r -> tuple(l[0], l[1], l[2], l[3], l[4], l[5], l[6], r[1]) }
//        .join(SEQUENZA_FITS.out.alt_solutions).map { l, r -> tuple(l[0], l[1], l[2], l[3], l[4], l[5], l[6], l[7], r[1]) }
//        .join(ONCOCNV.out.profile).map { l, r -> tuple(l[0], l[1], l[2], l[3], l[4], l[5], l[6], l[7], l[8], r[1]) }
//        .join(FACETS_CNV.out.vcf).map { l, r -> tuple(l[0], l[1], l[2], l[3], l[4], l[5], l[6], l[7], l[8], l[9], r[1]) }
//    ch_copyncat_input = CNVKIT_CALL.out.cns
//        .join(CNVKIT_BATCH.out.cns)
//        .map { meta, call_cns, batch_cns ->
//            tuple(meta, call_cns, batch_cns)
//        }
//        .join(CNVKIT_EXPORT.out.output)
//        .map { meta, call_cns, batch_cns, export_out ->
//            tuple(meta, call_cns, batch_cns, export_out)
//        }
//        .join(CONTROLFREEC_OT_FREEC.out.CNV)
//        .map { meta, call_cns, batch_cns, export_out, freec_cnv ->
//            tuple(meta, call_cns, batch_cns, export_out, freec_cnv)
//        }
//        .join(CONTROLFREEC_OT_FREEC.out.config)
//        .map { meta, call_cns, batch_cns, export_out, freec_cnv, freec_cfg ->
//            tuple(meta, call_cns, batch_cns, export_out, freec_cnv, freec_cfg)
//        }
//        .join(SEQUENZA_FITS.out.segments)
//        .map { meta, call_cns, batch_cns, export_out, freec_cnv, freec_cfg, seq_segments ->
//            tuple(meta, call_cns, batch_cns, export_out, freec_cnv, freec_cfg, seq_segments)
//        }
//        .join(SEQUENZA_FITS.out.confints_cp)
//        .map { meta, call_cns, batch_cns, export_out, freec_cnv, freec_cfg, seq_segments, seq_confints ->
//            tuple(meta, call_cns, batch_cns, export_out, freec_cnv, freec_cfg, seq_segments, seq_confints)
//        }
//        .join(SEQUENZA_FITS.out.alt_solutions)
//        .map { meta, call_cns, batch_cns, export_out, freec_cnv, freec_cfg, seq_segments, seq_confints, seq_alt ->
//            tuple(meta, call_cns, batch_cns, export_out, freec_cnv, freec_cfg, seq_segments, seq_confints, seq_alt)
//        }
//        .join(ONCOCNV.out.profile)
//        .map { meta, call_cns, batch_cns, export_out, freec_cnv, freec_cfg, seq_segments, seq_confints, seq_alt, oncocnv_profile ->
//            tuple(meta, call_cns, batch_cns, export_out, freec_cnv, freec_cfg, seq_segments, seq_confints, seq_alt, oncocnv_profile)
//        }
//        .join(FACETS_CNV.out.vcf)
//        .map { meta, call_cns, batch_cns, export_out, freec_cnv, freec_cfg, seq_segments, seq_confints, seq_alt, oncocnv_profile, facets_vcf ->
//            tuple(meta, call_cns, batch_cns, export_out, freec_cnv, freec_cfg, seq_segments, seq_confints, seq_alt, oncocnv_profile, facets_vcf)
//        }
//    //
//    // MODULE: Run CopyNcat (merge CNV calls from CNVkit, ControlFreec, Sequenza, OncoCNV, FACETS)
//    //
//    COPYNCAT(ch_copyncat_input)
//    ch_versions = ch_versions.mix(COPYNCAT.out.versions.first())

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'software_versions.yml', sort: true, newLine: true)
        .set { ch_col_vers }

    emit:

    versions            = ch_col_vers
    sam_mpileup         = ch_sam_mpileup
    bcf_mpileup         = ch_bcf_mpileup
    multiqc_files       = ch_multiqc_files
//    copyncat_tsv        = COPYNCAT.out.tsv
//    copyncat_summary    = COPYNCAT.out.summary

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                                THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
