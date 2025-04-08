process SEQUENZA_SEQZ {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/sequenza:3.0.0' :
        'blancojmskcc/sequenza:3.0.0' }"

    input:
    tuple val(meta),  path(tbam), path(tbai)
    tuple val(meta1), path(fasta)
    tuple val(meta2), path(fastai)
    path (nbam)
    path (nbai)
    path (wigfile)

    output:
    tuple val(meta), path("*_bin50.seqz.gz"), emit: seqz
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """ 
    sequenza-utils bam2seqz --normal ${nbam} --tumor ${tbam} --fasta ${fasta} -gc ${wigfile} -o ${prefix}_1_DELETE.seqz.gz -C 1
    sequenza-utils bam2seqz --normal ${nbam} --tumor ${tbam} --fasta ${fasta} -gc ${wigfile} -o ${prefix}_2_DELETE.seqz.gz -C 2
    sequenza-utils bam2seqz --normal ${nbam} --tumor ${tbam} --fasta ${fasta} -gc ${wigfile} -o ${prefix}_3_DELETE.seqz.gz -C 3
    sequenza-utils bam2seqz --normal ${nbam} --tumor ${tbam} --fasta ${fasta} -gc ${wigfile} -o ${prefix}_4_DELETE.seqz.gz -C 4
    sequenza-utils bam2seqz --normal ${nbam} --tumor ${tbam} --fasta ${fasta} -gc ${wigfile} -o ${prefix}_5_DELETE.seqz.gz -C 5
    sequenza-utils bam2seqz --normal ${nbam} --tumor ${tbam} --fasta ${fasta} -gc ${wigfile} -o ${prefix}_6_DELETE.seqz.gz -C 6
    sequenza-utils bam2seqz --normal ${nbam} --tumor ${tbam} --fasta ${fasta} -gc ${wigfile} -o ${prefix}_7_DELETE.seqz.gz -C 7
    sequenza-utils bam2seqz --normal ${nbam} --tumor ${tbam} --fasta ${fasta} -gc ${wigfile} -o ${prefix}_8_DELETE.seqz.gz -C 8
    sequenza-utils bam2seqz --normal ${nbam} --tumor ${tbam} --fasta ${fasta} -gc ${wigfile} -o ${prefix}_9_DELETE.seqz.gz -C 9
    sequenza-utils bam2seqz --normal ${nbam} --tumor ${tbam} --fasta ${fasta} -gc ${wigfile} -o ${prefix}_10_DELETE.seqz.gz -C 10
    sequenza-utils bam2seqz --normal ${nbam} --tumor ${tbam} --fasta ${fasta} -gc ${wigfile} -o ${prefix}_11_DELETE.seqz.gz -C 11
    sequenza-utils bam2seqz --normal ${nbam} --tumor ${tbam} --fasta ${fasta} -gc ${wigfile} -o ${prefix}_12_DELETE.seqz.gz -C 12
    sequenza-utils bam2seqz --normal ${nbam} --tumor ${tbam} --fasta ${fasta} -gc ${wigfile} -o ${prefix}_13_DELETE.seqz.gz -C 13
    sequenza-utils bam2seqz --normal ${nbam} --tumor ${tbam} --fasta ${fasta} -gc ${wigfile} -o ${prefix}_14_DELETE.seqz.gz -C 14
    sequenza-utils bam2seqz --normal ${nbam} --tumor ${tbam} --fasta ${fasta} -gc ${wigfile} -o ${prefix}_15_DELETE.seqz.gz -C 15
    sequenza-utils bam2seqz --normal ${nbam} --tumor ${tbam} --fasta ${fasta} -gc ${wigfile} -o ${prefix}_16_DELETE.seqz.gz -C 16
    sequenza-utils bam2seqz --normal ${nbam} --tumor ${tbam} --fasta ${fasta} -gc ${wigfile} -o ${prefix}_17_DELETE.seqz.gz -C 17
    sequenza-utils bam2seqz --normal ${nbam} --tumor ${tbam} --fasta ${fasta} -gc ${wigfile} -o ${prefix}_18_DELETE.seqz.gz -C 18
    sequenza-utils bam2seqz --normal ${nbam} --tumor ${tbam} --fasta ${fasta} -gc ${wigfile} -o ${prefix}_19_DELETE.seqz.gz -C 19
    sequenza-utils bam2seqz --normal ${nbam} --tumor ${tbam} --fasta ${fasta} -gc ${wigfile} -o ${prefix}_20_DELETE.seqz.gz -C 20
    sequenza-utils bam2seqz --normal ${nbam} --tumor ${tbam} --fasta ${fasta} -gc ${wigfile} -o ${prefix}_21_DELETE.seqz.gz -C 21
    sequenza-utils bam2seqz --normal ${nbam} --tumor ${tbam} --fasta ${fasta} -gc ${wigfile} -o ${prefix}_22_DELETE.seqz.gz -C 22
    sequenza-utils bam2seqz --normal ${nbam} --tumor ${tbam} --fasta ${fasta} -gc ${wigfile} -o ${prefix}_X_DELETE.seqz.gz -C X
    sequenza-utils seqz_merge -o ${prefix}_a_DELETE.seqz.gz -1 ${prefix}_1_DELETE.seqz.gz -2 ${prefix}_2_DELETE.seqz.gz
    sequenza-utils seqz_merge -o ${prefix}_b_DELETE.seqz.gz -1 ${prefix}_a_DELETE.seqz.gz -2 ${prefix}_3_DELETE.seqz.gz
    sequenza-utils seqz_merge -o ${prefix}_c_DELETE.seqz.gz -1 ${prefix}_b_DELETE.seqz.gz -2 ${prefix}_4_DELETE.seqz.gz
    sequenza-utils seqz_merge -o ${prefix}_d_DELETE.seqz.gz -1 ${prefix}_c_DELETE.seqz.gz -2 ${prefix}_5_DELETE.seqz.gz
    sequenza-utils seqz_merge -o ${prefix}_e_DELETE.seqz.gz -1 ${prefix}_d_DELETE.seqz.gz -2 ${prefix}_6_DELETE.seqz.gz
    sequenza-utils seqz_merge -o ${prefix}_f_DELETE.seqz.gz -1 ${prefix}_e_DELETE.seqz.gz -2 ${prefix}_7_DELETE.seqz.gz
    sequenza-utils seqz_merge -o ${prefix}_g_DELETE.seqz.gz -1 ${prefix}_f_DELETE.seqz.gz -2 ${prefix}_8_DELETE.seqz.gz
    sequenza-utils seqz_merge -o ${prefix}_h_DELETE.seqz.gz -1 ${prefix}_g_DELETE.seqz.gz -2 ${prefix}_9_DELETE.seqz.gz
    sequenza-utils seqz_merge -o ${prefix}_i_DELETE.seqz.gz -1 ${prefix}_h_DELETE.seqz.gz -2 ${prefix}_10_DELETE.seqz.gz
    sequenza-utils seqz_merge -o ${prefix}_j_DELETE.seqz.gz -1 ${prefix}_i_DELETE.seqz.gz -2 ${prefix}_11_DELETE.seqz.gz
    sequenza-utils seqz_merge -o ${prefix}_k_DELETE.seqz.gz -1 ${prefix}_j_DELETE.seqz.gz -2 ${prefix}_12_DELETE.seqz.gz
    sequenza-utils seqz_merge -o ${prefix}_l_DELETE.seqz.gz -1 ${prefix}_k_DELETE.seqz.gz -2 ${prefix}_13_DELETE.seqz.gz
    sequenza-utils seqz_merge -o ${prefix}_m_DELETE.seqz.gz -1 ${prefix}_l_DELETE.seqz.gz -2 ${prefix}_14_DELETE.seqz.gz
    sequenza-utils seqz_merge -o ${prefix}_n_DELETE.seqz.gz -1 ${prefix}_m_DELETE.seqz.gz -2 ${prefix}_15_DELETE.seqz.gz
    sequenza-utils seqz_merge -o ${prefix}_o_DELETE.seqz.gz -1 ${prefix}_n_DELETE.seqz.gz -2 ${prefix}_16_DELETE.seqz.gz
    sequenza-utils seqz_merge -o ${prefix}_p_DELETE.seqz.gz -1 ${prefix}_o_DELETE.seqz.gz -2 ${prefix}_17_DELETE.seqz.gz
    sequenza-utils seqz_merge -o ${prefix}_q_DELETE.seqz.gz -1 ${prefix}_p_DELETE.seqz.gz -2 ${prefix}_18_DELETE.seqz.gz
    sequenza-utils seqz_merge -o ${prefix}_r_DELETE.seqz.gz -1 ${prefix}_q_DELETE.seqz.gz -2 ${prefix}_19_DELETE.seqz.gz
    sequenza-utils seqz_merge -o ${prefix}_s_DELETE.seqz.gz -1 ${prefix}_r_DELETE.seqz.gz -2 ${prefix}_20_DELETE.seqz.gz
    sequenza-utils seqz_merge -o ${prefix}_t_DELETE.seqz.gz -1 ${prefix}_s_DELETE.seqz.gz -2 ${prefix}_21_DELETE.seqz.gz
    sequenza-utils seqz_merge -o ${prefix}_u_DELETE.seqz.gz -1 ${prefix}_t_DELETE.seqz.gz -2 ${prefix}_22_DELETE.seqz.gz
    sequenza-utils seqz_merge -o ${prefix}.seqz.gz -1 ${prefix}_u_DELETE.seqz.gz -2 ${prefix}_X_DELETE.seqz.gz

    sequenza-utils \\
        seqz_binning \\
        --seqz ${prefix}.seqz.gz \\
        --window 50 \\
        -o ${prefix}_bin50.seqz.gz

    rm *DELETE*

    cat << EOF > run_sequenza.${meta.id}.R
    library(sequenza)

    ${meta.id} <- sequenza.extract("${meta.id}_bin50.seqz.gz", verbose = FALSE)

    CP <- sequenza.fit(${meta.id})

    sequenza.results(sequenza.extract = ${meta.id}, 
                    cp.table = CP, 
                    sample.id = "${meta.id}", 
                    out.dir = "${meta.id}.output")

    EOF

    Rscript run_sequenza.${meta.id}.R

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: \$(echo \$(sequenza-utils 2>&1) | sed 's/^.*is version //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_bin50.seqz.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: \$(echo \$(sequenza-utils 2>&1) | sed 's/^.*is version //; s/ .*\$//')
    END_VERSIONS
    """
}
