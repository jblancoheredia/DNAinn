process MERGE_REPS {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta),  path(vcfs)
    tuple val(meta1), path(calls)
    tuple val(meta2), path(repeatseqs)

    output:
    tuple val(meta), path("${meta.id}_repeatseq.vcf")      , emit: vcf
    tuple val(meta), path("${meta.id}_repeatseq.calls")    , emit: calls
    tuple val(meta), path("${meta.id}_repeatseq.repeatseq"), emit: repeatseq
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    head -n 8 ${vcfs[0]} > ${prefix}_repeatseq.vcf
    for vcf in ${vcfs}; do
        grep -v "^#" \$vcf >> ${prefix}_repeatseq.vcf || [ \$? -eq 1 ]
    done

    cat ${calls} > ${prefix}_repeatseq.calls

    cat ${repeatseqs} > ${prefix}_repeatseq.repeatseq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MergeReps: "1.1.0"
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_repeatseq.vcf
    touch ${prefix}_repeatseq.calls
    touch ${prefix}_repeatseq.repeatseq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MergeReps: "1.1.0"
    END_VERSIONS
    """
}

process MERGE_REPS_RAW {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta),  path(vcfs)
    tuple val(meta1), path(calls)
    tuple val(meta2), path(repeatseqs)

    output:
    tuple val(meta), path("${meta.id}_repeatseq.raw.vcf")      , emit: vcf
    tuple val(meta), path("${meta.id}_repeatseq.raw.calls")    , emit: calls
    tuple val(meta), path("${meta.id}_repeatseq.raw.repeatseq"), emit: repeatseq
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    head -n 8 ${vcfs[0]} > ${prefix}_repeatseq.raw.vcf
    for vcf in ${vcfs}; do
        grep -v "^#" \$vcf >> ${prefix}_repeatseq.raw.vcf || [ \$? -eq 1 ]
    done

    cat ${calls} > ${prefix}_repeatseq.raw.calls

    cat ${repeatseqs} > ${prefix}_repeatseq.raw.repeatseq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MergeReps: "1.1.0"
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_repeatseq.raw.vcf
    touch ${prefix}_repeatseq.raw.calls
    touch ${prefix}_repeatseq.raw.repeatseq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MergeReps: "1.1.0"
    END_VERSIONS
    """
}

process MERGE_REPS_CON {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta),  path(vcfs)
    tuple val(meta1), path(calls)
    tuple val(meta2), path(repeatseqs)

    output:
    tuple val(meta), path("${meta.id}_repeatseq.con.vcf")      , emit: vcf
    tuple val(meta), path("${meta.id}_repeatseq.con.calls")    , emit: calls
    tuple val(meta), path("${meta.id}_repeatseq.con.repeatseq"), emit: repeatseq
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    head -n 8 ${vcfs[0]} > ${prefix}_repeatseq.con.vcf
    for vcf in ${vcfs}; do
        grep -v "^#" \$vcf >> ${prefix}_repeatseq.con.vcf || [ \$? -eq 1 ]
    done

    cat ${calls} > ${prefix}_repeatseq.con.calls

    cat ${repeatseqs} > ${prefix}_repeatseq.con.repeatseq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MergeReps: "1.1.0"
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_repeatseq.con.vcf
    touch ${prefix}_repeatseq.con.calls
    touch ${prefix}_repeatseq.con.repeatseq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MergeReps: "1.1.0"
    END_VERSIONS
    """
}

process MERGE_REPS_DUP {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta),  path(vcfs)
    tuple val(meta1), path(calls)
    tuple val(meta2), path(repeatseqs)

    output:
    tuple val(meta), path("${meta.id}_repeatseq.dup.vcf")      , emit: vcf
    tuple val(meta), path("${meta.id}_repeatseq.dup.calls")    , emit: calls
    tuple val(meta), path("${meta.id}_repeatseq.dup.repeatseq"), emit: repeatseq
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    head -n 8 ${vcfs[0]} > ${prefix}_repeatseq.dup.vcf
    for vcf in ${vcfs}; do
        grep -v "^#" \$vcf >> ${prefix}_repeatseq.dup.vcf || [ \$? -eq 1 ]
    done

    cat ${calls} > ${prefix}_repeatseq.dup.calls

    cat ${repeatseqs} > ${prefix}_repeatseq.dup.repeatseq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MergeReps: "1.1.0"
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_repeatseq.dup.vcf
    touch ${prefix}_repeatseq.dup.calls
    touch ${prefix}_repeatseq.dup.repeatseq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MergeReps: "1.1.0"
    END_VERSIONS
    """
}

process MERGE_REPS_SIM {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta),  path(vcfs)
    tuple val(meta1), path(calls)
    tuple val(meta2), path(repeatseqs)

    output:
    tuple val(meta), path("${meta.id}_repeatseq.sim.vcf")      , emit: vcf
    tuple val(meta), path("${meta.id}_repeatseq.sim.calls")    , emit: calls
    tuple val(meta), path("${meta.id}_repeatseq.sim.repeatseq"), emit: repeatseq
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    head -n 8 ${vcfs[0]} > ${prefix}_repeatseq.sim.vcf
    for vcf in ${vcfs}; do
        grep -v "^#" \$vcf >> ${prefix}_repeatseq.sim.vcf || [ \$? -eq 1 ]
    done

    cat ${calls} > ${prefix}_repeatseq.sim.calls

    cat ${repeatseqs} > ${prefix}_repeatseq.sim.repeatseq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MergeReps: "1.1.0"
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_repeatseq.sim.vcf
    touch ${prefix}_repeatseq.sim.calls
    touch ${prefix}_repeatseq.sim.repeatseq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MergeReps: "1.1.0"
    END_VERSIONS
    """
}
