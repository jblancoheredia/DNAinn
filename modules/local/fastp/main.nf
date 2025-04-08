process FASTP {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::fastp=0.23.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastp:0.23.4--h5f740d0_0' :
        'quay.io/biocontainers/fastp:0.23.4--h5f740d0_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.fastp.fastq.gz') , emit: reads
    tuple val(meta), path('*.json')           , emit: json
    tuple val(meta), path('*.html')           , emit: html
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    fastp \\
        --in1 ${reads[0]} \\
        --out1 ${prefix}_1.unfiltered.fastq.gz \\
        ${meta.single_end ? "" : "--in2 ${reads[1]}"} \\
        ${meta.single_end ? "" : "--out2 ${prefix}_2.unfiltered.fastq.gz"} \\
        --json ${prefix}.fastp.json \\
        --html ${prefix}.fastp.html \\
        --thread $task.cpus \\
        --umi \\
        --umi_loc per_read \\
        --umi_len 3 \\
        --trim_front1 7 \\
        --trim_front2 7 \\
        $args

    zcat ${prefix}_1.unfiltered.fastq.gz | \\
    awk 'BEGIN{OFS="\\n"} 
    {
        header=\$0; seq=\$2; plus=\$3; qual=\$4;
        if(NR%4==1) {
            if(header !~ /[ACTG]{3}_[ACTG]{3}/) {
                getline; getline; getline;
                next;
            }
        }
        gsub("_", "-", header)
        print header;
        getline; print;
        getline; print;
        getline; print;
    }' | gzip > ${prefix}_1.fastp.fastq.gz

    if [ -f "${prefix}_2.unfiltered.fastq.gz" ]; then
        zcat ${prefix}_2.unfiltered.fastq.gz | \\
        awk 'BEGIN{OFS="\\n"} 
        {
            header=\$0; seq=\$2; plus=\$3; qual=\$4;
            if(NR%4==1) {
                if(header !~ /[ACTG]{3}_[ACTG]{3}/) {
                    getline; getline; getline;
                    next;
                }
            }
            gsub("_", "-", header)
            print header;
            getline; print;
            getline; print;
            getline; print;
        }' | gzip > ${prefix}_2.fastp.fastq.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_1.fastp.fastq.gz
    ${meta.single_end ? "" : "touch ${prefix}_2.fastp.fastq.gz"}
    touch ${prefix}.fastp.json
    touch ${prefix}.fastp.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
    END_VERSIONS
    """
}
