process HLAIMGT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/imgt_hla:1.0.0' :
        'blancojmskcc/imgt_hla:1.0.0' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta), path(bai)
    path(bwa2)

    output:
    tuple val(meta), path("*.imgt.tsv"), emit: tsv
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools \\
        view \\
        -h ${bam} | \\
        grep -v "@" | \\
        awk '\$3 != \"*\" {count[\$3]++} END {for (hla in count) print hla \"\\t\" count[hla]}' | \\
        sort > HLA_COUNTS.tsv

    awk ' \$2 >= 2 { print \$1 }' HLA_COUNTS.tsv > HLA_ALLELE.tsv

    grep -f HLA_ALLELE.tsv ./imgt_hla/hla_gen.fa.hla > HLA_GREPNG.tsv

    cp "${workflow.projectDir}/bin/imgt_table.py" .

    python imgt_table.py

    mv imgt.tsv ${prefix}.imgt.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hlaimgt: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.imgt.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hlaimgt: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
