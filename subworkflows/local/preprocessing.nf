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
include { DOWNSAMPLINGS_COUNT                                                                                                       } from '../../modules/local/downsamplings/count'
include { DOWNSAMPLINGS_SEQTK                                                                                                       } from '../../modules/local/downsamplings/seqtk'
include { validateInputSamplesheet                                                                                                  } from './utils_nfcore_dnainn_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                            RUN SUBWORKFLOW 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PREPROCESSING {

    take:
    ch_samplesheet
    
    main:
    ch_reports          = Channel.empty()
    ch_versions         = Channel.empty()
    ch_multiqc_files    = Channel.empty()

    //
    // Create channel from input file provided through params.input
    //
    Channel
        .fromSamplesheet("input")
        .map {
            meta, fastq_1, fastq_2 ->
                if (!fastq_2) {
                    return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                } else {
                    return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                }
        }
        .groupTuple()
        .map {
            validateInputSamplesheet(it)
        }
        .branch {
            meta, fastqs ->
                single  : fastqs.size() == 1
                    return [ meta, fastqs.flatten() ]
                multiple: fastqs.size() > 1
                    return [ meta, fastqs.flatten() ]
        }
        .set { ch_fastq }

    ch_fastq_single   = ch_fastq.single.map   {meta, fastqs -> addReadgroupToMeta(meta, fastqs)}
    ch_fastq_multiple = ch_fastq.multiple.map {meta, fastqs -> addReadgroupToMeta(meta, fastqs)}

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                      FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Add readgroup to meta and remove lane
def addReadgroupToMeta(meta, files) {
    def flowcell = flowcellLaneFromFastq(files[0])
    def CN = params.seq_center ? "CN:${params.seq_center}\\t" : ''
    def read_group = "\"@RG\\tID:${meta.id}\\t${CN}PU:${flowcell}\\tSM:${meta.patient}_${meta.id}\\tLB:${meta.id}\\tPL:${params.seq_platform}\""
    meta  = meta + [read_group: read_group.toString()]

    return [ meta, files ]
}

// Parse first line of a FASTQ file, return the flowcell id and lane number.
def flowcellLaneFromFastq(path) {
    def line
    path.withInputStream {
        InputStream gzipStream = new java.util.zip.GZIPInputStream(it)
        Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
        BufferedReader buffered = new BufferedReader(decoder)
        line = buffered.readLine()
    }
    assert line.startsWith('@')
    line = line.substring(1)
    def fields = line.split(':')
    String fcid

    if (fields.size() >= 7) {
        // CASAVA 1.8+ format, from  https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.htm
        // "@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>:<UMI> <read>:<is filtered>:<control number>:<index>"
        fcid = fields[2]
    } else if (fields.size() == 5) {
        fcid = fields[0]
    }
    return fcid
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                                THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
