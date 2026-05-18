process LINE_PROBE_CANDIDATES_RAW {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'docker://community.wave.seqera.io/library/bedtools_samtools:a2e619e7b9f7645b' :
        'community.wave.seqera.io/library/bedtools_samtools:a2e619e7b9f7645b' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path line_annot

    output:
    tuple val(meta), path("*.line_probe.candidates.raw.tsv"), emit: candidates
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    bash ${moduleDir}/resources/usr/bin/line_probe_candidates.sh \\
        ${bam} ${line_annot} ${prefix} Raw raw

    cat <<END_VERSIONS > versions.yml
"${task.process}":
    samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    bedtools: \$(bedtools --version | sed 's/bedtools v//')
END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.line_probe.candidates.raw.tsv

    cat <<END_VERSIONS > versions.yml
"${task.process}":
    samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    bedtools: \$(bedtools --version | sed 's/bedtools v//')
END_VERSIONS
    """
}

process LINE_PROBE_CANDIDATES_CON {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'docker://community.wave.seqera.io/library/bedtools_samtools:a2e619e7b9f7645b' :
        'community.wave.seqera.io/library/bedtools_samtools:a2e619e7b9f7645b' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path line_annot

    output:
    tuple val(meta), path("*.line_probe.candidates.con.tsv"), emit: candidates
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    bash ${moduleDir}/resources/usr/bin/line_probe_candidates.sh \\
        ${bam} ${line_annot} ${prefix} Consensus con

    cat <<END_VERSIONS > versions.yml
"${task.process}":
    samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    bedtools: \$(bedtools --version | sed 's/bedtools v//')
END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.line_probe.candidates.con.tsv

    cat <<END_VERSIONS > versions.yml
"${task.process}":
    samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    bedtools: \$(bedtools --version | sed 's/bedtools v//')
END_VERSIONS
    """
}

process LINE_PROBE_CANDIDATES_DUP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'docker://community.wave.seqera.io/library/bedtools_samtools:a2e619e7b9f7645b' :
        'community.wave.seqera.io/library/bedtools_samtools:a2e619e7b9f7645b' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path line_annot

    output:
    tuple val(meta), path("*.line_probe.candidates.dup.tsv"), emit: candidates
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    bash ${moduleDir}/resources/usr/bin/line_probe_candidates.sh \\
        ${bam} ${line_annot} ${prefix} Duplex dup

    cat <<END_VERSIONS > versions.yml
"${task.process}":
    samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    bedtools: \$(bedtools --version | sed 's/bedtools v//')
END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.line_probe.candidates.dup.tsv

    cat <<END_VERSIONS > versions.yml
"${task.process}":
    samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    bedtools: \$(bedtools --version | sed 's/bedtools v//')
END_VERSIONS
    """
}

process LINE_PROBE_CANDIDATES_SIM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'docker://community.wave.seqera.io/library/bedtools_samtools:a2e619e7b9f7645b' :
        'community.wave.seqera.io/library/bedtools_samtools:a2e619e7b9f7645b' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path line_annot

    output:
    tuple val(meta), path("*.line_probe.candidates.sim.tsv"), emit: candidates
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    bash ${moduleDir}/resources/usr/bin/line_probe_candidates.sh \\
        ${bam} ${line_annot} ${prefix} Simplex sim

    cat <<END_VERSIONS > versions.yml
"${task.process}":
    samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    bedtools: \$(bedtools --version | sed 's/bedtools v//')
END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.line_probe.candidates.sim.tsv

    cat <<END_VERSIONS > versions.yml
"${task.process}":
    samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    bedtools: \$(bedtools --version | sed 's/bedtools v//')
END_VERSIONS
    """
}
