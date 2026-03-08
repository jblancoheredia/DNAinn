/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                    IMPORT PLUGINS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap                                                              } from 'plugin/nf-validation'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MULTIQC                                                                       } from '../modules/nf-core/multiqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                 IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { COPYNUMBERALT                                                                 } from '../subworkflows/local/copynumberalt'
include { DEDUPANDRECAL                                                                 } from '../subworkflows/local/dedupandrecal'
include { IMMUNONCOLOGY                                                                 } from '../subworkflows/local/immunoncology'
include { PREPROCESSING                                                                 } from '../subworkflows/local/preprocessing'
include { STRCTRLVARNTS                                                                 } from '../subworkflows/local/strctrlvarnts'
include { TELOMEREFEATS                                                                 } from '../subworkflows/local/telomerefeats'
include { UMIPROCESSING                                                                 } from '../subworkflows/local/umiprocessing'
include { VARIANTDSCVRY                                                                 } from '../subworkflows/local/variantdscvry'
include { paramsSummaryMultiqc                                                          } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                                                        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                                                        } from '../subworkflows/local/utils_nfcore_dnainn_pipeline'
include { validateInputSamplesheet                                                      } from '../subworkflows/local/utils_nfcore_dnainn_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                            CREATE CUSTOM CHANNELS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_fai                                      = Channel.fromPath(params.fai).map                              { it -> [[id:it.Name], it] }.collect()
ch_bwa                                      = Channel.fromPath(params.bwa).map                              { it -> [[id:it.Name], it] }.collect()
ch_bwa2                                     = Channel.fromPath(params.bwa2).map                             { it -> [[id:it.Name], it] }.collect()
ch_dict                                     = Channel.fromPath(params.dict).map                             { it -> [[id:it.Name], it] }.collect()
ch_fasta                                    = Channel.fromPath(params.fasta).map                            { it -> [[id:it.Name], it] }.collect()
ch_msi_f                                    = Channel.fromPath(params.msi_list).map                         { it -> [[id:it.Name], it] }.collect()
ch_targets                                  = Channel.fromPath(params.targets).map                          { it -> [[id:it.Name], it] }.collect()
ch_intervals                                = Channel.fromPath(params.intervals).map                        { it -> [[id:it.Name], it] }.collect()
ch_known_sites                              = Channel.fromPath(params.known_sites).map                      { vcf -> def tbi = 
                                                                                                             file("${params.known_sites_tbi}") 
                                                                                                             [[id:vcf.Name], vcf, tbi] }.collect()
ch_targets_bed                              = Channel.fromPath(params.targets_bed).map                      { it -> [[id:it.Name], it] }.collect()
ch_normal_con_bam                           = Channel.fromPath(params.normal_con_bam).map                   { it -> [[id:it.Name], it] }.collect()
ch_normal_con_bai                           = Channel.fromPath(params.normal_con_bai).map                   { it -> [[id:it.Name], it] }.collect()
ch_bwa_dragen_hg38                          = Channel.fromPath(params.bwa_dragen_hg38).map                  { it -> [[id:it.Name], it] }.collect()
ch_cnvkit_reference                         = Channel.fromPath(params.cnvkit_reference).map                 { it -> [[id:it.Name], it] }.collect()
ch_intervals_gunzip                         = Channel.fromPath(params.intervals_bed_gunzip).map             { it -> [[id:it.Name], it] }.collect()
ch_cnvkit_antitarget                        = Channel.fromPath(params.cnvkit_antitarget).map                { it -> [[id:it.Name], it] }.collect()
ch_gatk_interval_list                       = Channel.fromPath(params.gatk_interval_list).map               { it -> [[id:it.Name], it] }.collect()
ch_intervals_gunzip_index                   = Channel.fromPath(params.intervals_bed_gunzip_index).map       { it -> [[id:it.Name], it] }.collect()


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                              RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DNAINN {

    take:
    ch_samplesheet

    main:

    //
    // SUBWORKFLOW: PRE processing
    //
    PREPROCESSING(
        ch_samplesheet
    )
    ch_fastqs                   = PREPROCESSING.out.fastqs
    ch_versions					= PREPROCESSING.out.versions
    ch_multiqc_files			= PREPROCESSING.out.multiqc_files

    if (params.run_umiprocessing) {
        
        //
        // SUBWORKFLOW: UMI processing
        //
        UMIPROCESSING(
            ch_fai,
            ch_bwa2,
            ch_dict,
            ch_fasta,
            ch_msi_f,
            ch_fastqs
        )
        ch_ubam						= UMIPROCESSING.out.ubam
        ch_raw_bam					= UMIPROCESSING.out.raw_bam
        ch_raw_bai					= UMIPROCESSING.out.raw_bai
        ch_raw_baix                 = UMIPROCESSING.out.raw_baix
        ch_versions					= ch_versions.mix(UMIPROCESSING.out.versions)
        ch_bam_duplex               = UMIPROCESSING.out.duplex_bam
        ch_bam_grouped				= UMIPROCESSING.out.group_bam
        ch_split_reads              = UMIPROCESSING.out.split_reads
        ch_bam_bai_dedup            = UMIPROCESSING.out.bam_dedup
        ch_multiqc_files			= ch_multiqc_files.mix(UMIPROCESSING.out.multiqc_files)
        ch_bam_finalized			= UMIPROCESSING.out.finalized_bam
        ch_split_contigs            = UMIPROCESSING.out.split_contigs
        ch_reads_finalized          = UMIPROCESSING.out.reads_finalized


    } else {

        //
        // SUBWORKFLOW: Deduplication & Recalibration
        //
        DEDUPANDRECAL(
            ch_fai,
            ch_bwa2,
            ch_dict,
            ch_fasta,
            ch_msi_f,
            ch_fastqs
        )
        ch_raw_bam					= DEDUPANDRECAL.out.raw_bam
        ch_raw_bai					= DEDUPANDRECAL.out.raw_bai
        ch_versions					= ch_versions.mix(DEDUPANDRECAL.out.versions)
        ch_split_reads              = DEDUPANDRECAL.out.split_reads
        ch_multiqc_files			= ch_multiqc_files.mix(DEDUPANDRECAL.out.multiqc_files)
        ch_split_contigs            = DEDUPANDRECAL.out.split_contigs
        ch_bam_bai_dedup            = DEDUPANDRECAL.out.bam_dedup
        ch_bam_finalized			= DEDUPANDRECAL.out.bam_final
        ch_reads_finalized			= DEDUPANDRECAL.out.reads_final

    }

    //
    // Pair tumour samples with normal samples by patient meta key if only tumour use backup normal
    //
    def control_normal_bam = file(params.normal_con_bam)
    def control_normal_bai = file(params.normal_con_bai)

    ch_bam_finalized
        .map { meta, bam, bai -> [meta.patient, meta, bam, bai] }
        .groupTuple(by: 0)
        .map { patient, meta_list, bam_list, bai_list -> 
            def samples = []
            meta_list.eachWithIndex { meta, i ->
                samples << [meta, bam_list[i], bai_list[i]]
            }

            def tumour = samples.find { it[0].sample_type == 'tumour' }
            def normal = samples.find { it[0].sample_type == 'normal' }

            if (tumour && normal) {
                return [tumour, normal]
            } else if (tumour && tumour[0].matched_normal == 0) {
                return [tumour, null]
            } else {
                return null
            }
        }
        .filter { it != null }
        .map { tumour, normal ->
            def meta_t = tumour[0]
            def bam_t  = tumour[1]
            def bai_t  = tumour[2]

            if (normal) {
                def bam_n = normal[1]
                def bai_n = normal[2]
                return [meta_t, bam_t, bai_t, bam_n, bai_n]
            } else {
                return [meta_t, bam_t, bai_t, control_normal_bam, control_normal_bai]
            }
        }
        .set { ch_bam_pairs }

    ch_bam_bai_dedup
        .map { meta, bam, bai -> [meta.patient, meta, bam, bai] }
        .groupTuple(by: 0)
        .map { patient, meta_list, bam_list, bai_list -> 
            def samples = []
            meta_list.eachWithIndex { meta, i ->
                samples << [meta, bam_list[i], bai_list[i]]
            }

            def tumour = samples.find { it[0].sample_type == 'tumour' }
            def normal = samples.find { it[0].sample_type == 'normal' }

            if (tumour && normal) {
                return [tumour, normal]
            } else if (tumour && tumour[0].matched_normal == 0) {
                return [tumour, null]
            } else {
                return null
            }
        }
        .filter { it != null }
        .map { tumour, normal ->
            def meta_t = tumour[0]
            def bam_t  = tumour[1]
            def bai_t  = tumour[2]

            if (normal) {
                def bam_n = normal[1]
                def bai_n = normal[2]
                return [meta_t, bam_t, bai_t, bam_n, bai_n]
            } else {
                return [meta_t, bam_t, bai_t, control_normal_bam, control_normal_bai]
            }
        }
        .set { ch_bam_dedup_pairs }

    if (params.run_copynumberalt) {

        //
        // SUBWORKFLOW: Run copy number workflow (optional)
        //
        COPYNUMBERALT(
            ch_fai,
            ch_fasta,
            ch_raw_bam,
            ch_raw_bai,
            ch_bam_pairs,
            ch_intervals,
            ch_targets_bed,
            ch_bam_finalized,
            ch_normal_con_bam,
            ch_normal_con_bai,
            ch_cnvkit_reference,
            ch_cnvkit_antitarget
        )
//        ch_versions = ch_versions.mix(COPYNUMBERALT.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(COPYNUMBERALT.out.multiqc_files)
        ch_sam_mpileup = COPYNUMBERALT.out.sam_mpileup
        ch_bcf_mpileup = COPYNUMBERALT.out.bcf_mpileup
    } else {
        ch_sam_mpileup = Channel.empty()
        ch_bcf_mpileup = Channel.empty()
    }

    if (params.run_variantdscvry) {

        //
        // SUBWORKFLOW: Run variant discovery (optional)
        //
        VARIANTDSCVRY(
            ch_fai,
            ch_dict,
            ch_fasta,
            ch_targets,
            ch_intervals,
            ch_bam_pairs,
            ch_known_sites,
            ch_bam_finalized,
            ch_gatk_interval_list
        )
        ch_variants = VARIANTDSCVRY.out.variants
        ch_versions = ch_versions.mix(VARIANTDSCVRY.out.versions)
//    ch_multiqc_files = ch_multiqc_files.mix(VARIANTDSCVRY.out.multiqc_files)
    } else {
        ch_variants = Channel.empty()
    }

    if (params.run_telomerefeats) {

        //
        // SUBWORKFLOW: Run extract telomeric features (optional)
        //
        TELOMEREFEATS(
            ch_fai,
            ch_fasta,
            ch_fastqs,
            ch_bwa_dragen_hg38
        )
        ch_telhun = TELOMEREFEATS.out.telhun
        ch_versions = ch_versions.mix(TELOMEREFEATS.out.versions)
//    ch_multiqc_files = ch_multiqc_files.mix(TELOMEREFEATS.out.multiqc_files)
    } else {
        ch_telhun = Channel.empty()
    }

    if (params.run_immunoncology) {

        //
        // SUBWORKFLOW: Run immunoncology analysis (optional)
        //
        IMMUNONCOLOGY(
            ch_raw_bam,
            ch_raw_bai,
            ch_fastqs,
            ch_reads_finalized
        )
        ch_versions = ch_versions.mix(IMMUNONCOLOGY.out.versions)
//    ch_multiqc_files = ch_multiqc_files.mix(TELOMEREFEATS.out.multiqc_files)
    } else {
        ch_immuonco = Channel.empty()
    }

    if (params.run_strctrlvarnts) {

    //
    // SUBWORKFLOW: Run structural variants discovery
    //
    STRCTRLVARNTS(
        ch_bwa,
        ch_fai,
        ch_bwa2,
        ch_dict,
        ch_fasta,
        ch_bcf_mpileup,
        ch_known_sites,
        ch_split_reads,
        ch_split_contigs,
        ch_bam_dedup_pairs,
        ch_reads_finalized,
        ch_intervals_gunzip,
        ch_intervals_gunzip_index
    )
    ch_versions = ch_versions.mix(STRCTRLVARNTS.out.versions)
    } else {
        ch_sv_annotated = Channel.empty()
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'software_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : 
                                            file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    ch_multiqc_files = ch_multiqc_files
        .flatten() 
        .map { it instanceof Map || it instanceof List ? it[-1] : it } 
        .mix(
            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
            ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false),
            ch_collated_versions
        )
    
    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                       THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
