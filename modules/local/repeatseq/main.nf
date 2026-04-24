process REPEATSEQ {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/repeatseq:0.8.2' :
        'blancojmskcc/repeatseq:0.8.2' }"

    input:
    tuple val(meta), val(chunk), path(bam, stageAs: 'repeatseq_inputs/chunk.bam'), path(bai, stageAs: 'repeatseq_inputs/chunk.bam.bai')
    tuple val(meta1), path(fasta, stageAs: 'reference/genome.fasta')
    tuple val(meta2), path(fai, stageAs: 'reference/genome.fasta.fai')
    path(rep_regions, stageAs: 'reference/repeat_regions')

    output:
    tuple val(meta), path("*_repeatseq.vcf"), optional: true,       emit: vcf
    tuple val(meta), path("*_repeatseq.calls"), optional: true,     emit: calls
    tuple val(meta), path("*_repeatseq.repeatseq"), optional: true, emit: repeatseq
    path "versions.yml",                                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    set +e
    repeatseq \\
        -repeatseq \\
        -calls \\
        ${bam} \\
        ${fasta} \\
        ${rep_regions}
    repeatseq_exit=\$?
    set -e
    
    if [ \$repeatseq_exit -eq 0 ]; then
        mv *.vcf ${prefix}_${chunk}_repeatseq.vcf
        mv *.calls ${prefix}_${chunk}_repeatseq.calls
        mv *.repeatseq ${prefix}_${chunk}_repeatseq.repeatseq
    else
        echo "[WARN] repeatseq failed for ${prefix}_${chunk} (exit \$repeatseq_exit); continuing without this chunk." >&2
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatseq: \$(repeatseq -v 2>&1 | grep "RepeatSeq v" | sed 's/.*RepeatSeq v//')
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${chunk}_repeatseq.vcf
    touch ${prefix}_${chunk}_repeatseq.calls
    touch ${prefix}_${chunk}_repeatseq.repeatseq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatseq: \$(repeatseq -v 2>&1 | grep "RepeatSeq v" | sed 's/.*RepeatSeq v//')
    END_VERSIONS
    """
}

process REPEATSEQ_RAW {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/repeatseq:0.8.2' :
        'blancojmskcc/repeatseq:0.8.2' }"

    input:
    tuple val(meta), val(chunk), path(bam, stageAs: 'repeatseq_inputs/chunk.bam'), path(bai, stageAs: 'repeatseq_inputs/chunk.bam.bai')
    tuple val(meta1), path(fasta, stageAs: 'reference/genome.fasta')
    tuple val(meta2), path(fai, stageAs: 'reference/genome.fasta.fai')
    path(rep_regions, stageAs: 'reference/repeat_regions')

    output:
    tuple val(meta), path("*_repeatseq.raw.vcf")      , optional: true, emit: vcf
    tuple val(meta), path("*_repeatseq.raw.calls")    , optional: true, emit: calls
    tuple val(meta), path("*_repeatseq.raw.repeatseq"), optional: true, emit: repeatseq
    path "versions.yml"                                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    set +e
    repeatseq \\
        -repeatseq \\
        -calls \\
        ${bam} \\
        ${fasta} \\
        ${rep_regions}
    repeatseq_exit=\$?
    set -e
    
    if [ \$repeatseq_exit -eq 0 ]; then
        mv *.vcf ${prefix}_${chunk}_repeatseq.raw.vcf
        mv *.calls ${prefix}_${chunk}_repeatseq.raw.calls
        mv *.repeatseq ${prefix}_${chunk}_repeatseq.raw.repeatseq
    else
        echo "[WARN] repeatseq failed for ${prefix}_${chunk} (exit \$repeatseq_exit); continuing without this chunk." >&2
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatseq: \$(repeatseq -v 2>&1 | grep "RepeatSeq v" | sed 's/.*RepeatSeq v//')
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${chunk}_repeatseq.raw.vcf
    touch ${prefix}_${chunk}_repeatseq.raw.calls
    touch ${prefix}_${chunk}_repeatseq.raw.repeatseq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatseq: \$(repeatseq -v 2>&1 | grep "RepeatSeq v" | sed 's/.*RepeatSeq v//')
    END_VERSIONS
    """
}

process REPEATSEQ_CON {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/repeatseq:0.8.2' :
        'blancojmskcc/repeatseq:0.8.2' }"

    input:
    tuple val(meta), val(chunk), path(bam, stageAs: 'repeatseq_inputs/chunk.bam'), path(bai, stageAs: 'repeatseq_inputs/chunk.bam.bai')
    tuple val(meta1), path(fasta, stageAs: 'reference/genome.fasta')
    tuple val(meta2), path(fai, stageAs: 'reference/genome.fasta.fai')
    path(rep_regions, stageAs: 'reference/repeat_regions')

    output:
    tuple val(meta), path("*_repeatseq.con.vcf")      , optional: true, emit: vcf
    tuple val(meta), path("*_repeatseq.con.calls")    , optional: true, emit: calls
    tuple val(meta), path("*_repeatseq.con.repeatseq"), optional: true, emit: repeatseq
    path "versions.yml"                                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    set +e
    repeatseq \\
        -repeatseq \\
        -calls \\
        ${bam} \\
        ${fasta} \\
        ${rep_regions}
    repeatseq_exit=\$?
    set -e
    
    if [ \$repeatseq_exit -eq 0 ]; then
        mv *.vcf ${prefix}_${chunk}_repeatseq.con.vcf
        mv *.calls ${prefix}_${chunk}_repeatseq.con.calls
        mv *.repeatseq ${prefix}_${chunk}_repeatseq.con.repeatseq
    else
        echo "[WARN] repeatseq failed for ${prefix}_${chunk} (exit \$repeatseq_exit); continuing without this chunk." >&2
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatseq: \$(repeatseq -v 2>&1 | grep "RepeatSeq v" | sed 's/.*RepeatSeq v//')
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${chunk}_repeatseq.con.vcf
    touch ${prefix}_${chunk}_repeatseq.con.calls
    touch ${prefix}_${chunk}_repeatseq.con.repeatseq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatseq: \$(repeatseq -v 2>&1 | grep "RepeatSeq v" | sed 's/.*RepeatSeq v//')
    END_VERSIONS
    """
}

process REPEATSEQ_DUP {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/repeatseq:0.8.2' :
        'blancojmskcc/repeatseq:0.8.2' }"

    input:
    tuple val(meta), val(chunk), path(bam, stageAs: 'repeatseq_inputs/chunk.bam'), path(bai, stageAs: 'repeatseq_inputs/chunk.bam.bai')
    tuple val(meta1), path(fasta, stageAs: 'reference/genome.fasta')
    tuple val(meta2), path(fai, stageAs: 'reference/genome.fasta.fai')
    path(rep_regions, stageAs: 'reference/repeat_regions')

    output:
    tuple val(meta), path("*_repeatseq.dup.vcf")      , optional: true, emit: vcf
    tuple val(meta), path("*_repeatseq.dup.calls")    , optional: true, emit: calls
    tuple val(meta), path("*_repeatseq.dup.repeatseq"), optional: true, emit: repeatseq
    path "versions.yml"                                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    set +e
    repeatseq \\
        -repeatseq \\
        -calls \\
        ${bam} \\
        ${fasta} \\
        ${rep_regions}
    repeatseq_exit=\$?
    set -e
    
    if [ \$repeatseq_exit -eq 0 ]; then
        mv *.vcf ${prefix}_${chunk}_repeatseq.dup.vcf
        mv *.calls ${prefix}_${chunk}_repeatseq.dup.calls
        mv *.repeatseq ${prefix}_${chunk}_repeatseq.dup.repeatseq
    else
        echo "[WARN] repeatseq failed for ${prefix}_${chunk} (exit \$repeatseq_exit); continuing without this chunk." >&2
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatseq: \$(repeatseq -v 2>&1 | grep "RepeatSeq v" | sed 's/.*RepeatSeq v//')
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${chunk}_repeatseq.dup.vcf
    touch ${prefix}_${chunk}_repeatseq.dup.calls
    touch ${prefix}_${chunk}_repeatseq.dup.repeatseq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatseq: \$(repeatseq -v 2>&1 | grep "RepeatSeq v" | sed 's/.*RepeatSeq v//')
    END_VERSIONS
    """
}

process REPEATSEQ_SIM {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/repeatseq:0.8.2' :
        'blancojmskcc/repeatseq:0.8.2' }"

    input:
    tuple val(meta), val(chunk), path(bam, stageAs: 'repeatseq_inputs/chunk.bam'), path(bai, stageAs: 'repeatseq_inputs/chunk.bam.bai')
    tuple val(meta1), path(fasta, stageAs: 'reference/genome.fasta')
    tuple val(meta2), path(fai, stageAs: 'reference/genome.fasta.fai')
    path(rep_regions, stageAs: 'reference/repeat_regions')

    output:
    tuple val(meta), path("*_repeatseq.sim.vcf")      , optional: true, emit: vcf
    tuple val(meta), path("*_repeatseq.sim.calls")    , optional: true, emit: calls
    tuple val(meta), path("*_repeatseq.sim.repeatseq"), optional: true, emit: repeatseq
    path "versions.yml"                                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    set +e
    repeatseq \\
        -repeatseq \\
        -calls \\
        ${bam} \\
        ${fasta} \\
        ${rep_regions}
    repeatseq_exit=\$?
    set -e
    
    if [ \$repeatseq_exit -eq 0 ]; then
        mv *.vcf ${prefix}_${chunk}_repeatseq.sim.vcf
        mv *.calls ${prefix}_${chunk}_repeatseq.sim.calls
        mv *.repeatseq ${prefix}_${chunk}_repeatseq.sim.repeatseq
    else
        echo "[WARN] repeatseq failed for ${prefix}_${chunk} (exit \$repeatseq_exit); continuing without this chunk." >&2
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatseq: \$(repeatseq -v 2>&1 | grep "RepeatSeq v" | sed 's/.*RepeatSeq v//')
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${chunk}_repeatseq.sim.vcf
    touch ${prefix}_${chunk}_repeatseq.sim.calls
    touch ${prefix}_${chunk}_repeatseq.sim.repeatseq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatseq: \$(repeatseq -v 2>&1 | grep "RepeatSeq v" | sed 's/.*RepeatSeq v//')
    END_VERSIONS
    """
}
