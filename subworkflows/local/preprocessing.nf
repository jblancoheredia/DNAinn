/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                       PREPROCESSING SUBWORKFLOW                                                    
*******************************************************************************************************************************************************************************************************

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                             IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                                                                                                                    } from '../../modules/nf-core/fastqc/main'
include { CAT_FASTQ                                                                                                                 } from '../../modules/local/cat/fastq/main'
include { DOWNSAMPLINGS_COUNT                                                           } from '../modules/local/downsamplings/count'
include { DOWNSAMPLINGS_SEQTK                                                           } from '../modules/local/downsamplings/seqtk'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                            RUN SUBWORKFLOW 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PREPROCESSING {

    take:
    ch_fastq
    ch_fastq_single
    ch_fastq_multiple

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    //
    // MODULE: Run FastQC
    //
    FASTQC(ch_fastq.mix())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (ch_fastq_multiple)
    .reads
    .mix(ch_fastq_single)
    .set { ch_cat_fastq }
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)

    //
    // MODULE: Run in-house script for counting reads
    //
    DOWNSAMPLINGS_COUNT(ch_cat_fastq)
    ch_versions = ch_versions.mix(DOWNSAMPLINGS_COUNT.out.versions)
    ch_global_min_reads = DOWNSAMPLINGS_COUNT.out.total_reads
        .map { file -> 
            def count = file.text.trim()
            count.toInteger()
        }
        .collect()
        .map { counts -> 
            def min_count = counts.min()
            min_count
        }

    if (params.run_downsamplings) {

        if (params.downsampling_total_reads) {

            //
            // MODULE: Run Downsampling with seqtk
            //
            DOWNSAMPLINGS_SEQTK(ch_cat_fastq, params.downsampling_total_reads)
            ch_versions = ch_versions.mix(DOWNSAMPLINGS_SEQTK.out.versions)
            ch_downsampled_reads = DOWNSAMPLINGS_SEQTK.out.downsampled_reads

            ch_fastqs = ch_downsampled_reads

        } else {

            //
            // MODULE: Run Downsampling with seqtk
            //
            DOWNSAMPLINGS_SEQTK(ch_cat_fastq, ch_global_min_reads)
            ch_versions = ch_versions.mix(DOWNSAMPLINGS_SEQTK.out.versions)
            ch_downsampled_reads = DOWNSAMPLINGS_SEQTK.out.downsampled_reads

            ch_fastqs = ch_downsampled_reads

        }

    } else {

        ch_fastqs = ch_cat_fastq

    }

    emit:

    fastqs          = ch_fastqs
    versions        = ch_versions
    multiqc_files   = ch_multiqc_files 
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                                THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
