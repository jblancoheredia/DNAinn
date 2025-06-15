process SEQUENZA_FITS {
    tag "$meta.patient"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/sequenza:3.0.1' :
        'blancojmskcc/sequenza:3.0.1' }"

    input:
    tuple val(meta),  path(seqz)

    output:
	tuple val(meta), path("*_CN_bars.pdf")  			, emit: plot_cn_bars
	tuple val(meta), path("*_gc_plots.pdf")  			, emit: gc_plots
	tuple val(meta), path("*_segments.txt")  			, emit:	segments
	tuple val(meta), path("*_mutations.txt")  			, emit: mutations
	tuple val(meta), path("*_model_fit.pdf")  			, emit: plot_model_fit
	tuple val(meta), path("*_confints_CP.txt")  		, emit: confints_cp
	tuple val(meta), path("*_CP_contours.pdf")  		, emit: cp_contours
	tuple val(meta), path("*_genome_view.pdf")  		, emit: genome_view
	tuple val(meta), path("*_sequenza_log.txt")  		, emit: log
	tuple val(meta), path("*_alternative_fit.pdf")  	, emit:	plot_alt_fit
	tuple val(meta), path("*_chromosome_view.pdf")  	, emit: chr_view
	tuple val(meta), path("*_chromosome_depths.pdf")  	, emit: chr_depth
	tuple val(meta), path("*_sequenza_extract.RData")  	, emit: rdata_extract
	tuple val(meta), path("*_sequenza_cp_table.RData")  , emit: rdata_cp_table
	tuple val(meta), path("*_alternative_solutions.txt"), emit: alt_solutions
	path "versions.yml"             					, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def VERSION = '3.0.0'
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.patient}"
    """ 
    cat << EOF > run_sequenza.${prefix}.R
    Sys.setenv("VROOM_CONNECTION_SIZE" = 1000000000)

    library(sequenza)

    data.file <- "${seqz}"

    output.dir <- "sequenza_results"

    sample.id <- "${prefix}"

    seqz.data <- sequenza.extract(data.file)

    CP <- sequenza.fit(seqz.data)

    sequenza.results(sequenza.extract = seqz.data,
                 cp.table = CP,
                 sample.id = sample.id,
                 out.dir = output.dir)

    EOF

    Rscript run_sequenza.${prefix}.R

    mv sequenza_results/${prefix}* .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenza: ${VERSION}
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.patient}"
    """
    touch ${prefix}.tsv
    touch ${prefix}.pdf
    touch ${prefix}.RData

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenza: ${VERSION}
    END_VERSIONS
    """
}
