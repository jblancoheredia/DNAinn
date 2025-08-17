process TIDDIT_SV {
    tag "$meta.patient"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/tiddit:3.9.1' :
        'blancojmskcc/tiddit:3.9.1' }"

    input:
    tuple val(meta), path(input), path(input_index), path(nbam), path(nbai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(bwa_index)

    output:
    tuple val(meta), path("*.vcf")         , emit: vcf
    tuple val(meta), path("*.ploidies.tab"), emit: ploidy
    path  "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.patient}"
    def bwa_command = bwa_index ? "[[ -d $bwa_index ]] && for i in $bwa_index/*; do [[ -f $fasta && ! \"\$i\" =~ .*\"$fasta.\".* ]] && ln -s \$i ${fasta}.\${i##*.} || ln -s \$i .; done" : ""

    """
    $bwa_command

    tiddit \\
        --sv \\
        $args \\
        --threads $task.cpus \\
        --bam $input \\
        --ref $fasta \\
        -o ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tiddit: \$(echo \$(tiddit 2>&1) | sed 's/^.*tiddit-//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix ="${meta.patient}"
    """
    touch ${prefix}.tiddit.unfiltered.vcf
    touch ${prefix}.ploidies.tab

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tiddit: \$(echo \$(tiddit 2>&1) | sed 's/^.*tiddit-//; s/ .*\$//')
    END_VERSIONS
    """
}
