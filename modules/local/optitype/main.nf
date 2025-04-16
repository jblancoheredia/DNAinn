process OPTITYPE {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(reads)
    val(seq_type)
    path(img)

    output:
    tuple val(meta), path("*_coverage_plot.pdf", optional: true), emit: pdf
    tuple val(meta), path("*_result.tsv", optional: true)       , emit: tsv
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args   ?: ''
    prefix    = task.ext.prefix ?: "${meta.id}"
    """

    echo "${task.process}:
             optitype: 1.3" > versions.yml
    echo "OptiType failed for ${prefix}" > ${prefix}_result.tsv
    echo "OptiType failed for ${prefix}" > ${prefix}_coverage_plot.pdf

    mkdir 01_data/
    cp ${reads[0]} 01_data/
    cp ${reads[1]} 01_data/
    mkdir 03_outs/
    mkdir 02_code/
    cp ${img} 02_code/
    cd 02_code/

    (
        singularity \\
            exec \\
            --bind \$(pwd)/../01_data:/data \\
            --bind \$(pwd)/../03_outs:/outs \\
            --writable-tmpfs \\
            ${img} \\
            OptiTypePipeline.py \\
                -i /data/${reads[0]} /data/${reads[1]} \\
                --${seq_type} \\
                --prefix ${prefix} \\
                --outdir /outs/
    ) || echo "OptiType failed for ${prefix}, continuing anyway."

    cd ..

    [[ -f 03_outs/${prefix}_result.tsv ]] && mv 03_outs/${prefix}_result.tsv .
    [[ -f 03_outs/${prefix}_coverage_plot.pdf ]] && mv 03_outs/${prefix}_coverage_plot.pdf .

    rm -rf 01_data/
    rm -rf 02_code/
    rm -rf 03_outs/
    
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_result.tsv
    touch ${prefix}_coverage_plot.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        optitype: \$(cat \$(which OptiTypePipeline.py) | grep -e "Version:" | sed -e "s/Version: //g")
    END_VERSIONS
    """
}
