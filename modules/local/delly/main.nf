process DELLY {
    tag "$meta.patient"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/delly:1.3.3' :
        'blancojmskcc/delly:1.3.3' }"

    input:
    tuple val(meta),  path(tbam),  path(tbai), path(nbam), path(nbai)
    tuple val(meta1), path(fasta)
    tuple val(meta2), path(fai)
    path(exclude_bed)

    output:
    tuple val(meta), path("*.{csi,tbi}")            , emit: csi
    tuple val(meta), path("*.{bcf,delly.vcf.gz}")   , emit: bcf
    tuple val(meta), path("*.delly.unfiltered.vcf") , emit: vcf
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.patient}"
    def suffix = task.ext.suffix ?: "bcf"
    def exclude = exclude_bed ? "--exclude ${exclude_bed}" : ""
    def bcf_output = suffix == "bcf" ? "--outfile ${prefix}.bcf" : ""
    def bcf_filter = suffix == "bcf" ? "--outfile ${prefix}.filtered.bcf" : ""
    """
	cp ${tbam} tumor_modified.bam
	cp ${nbam} normal_modified.bam
	
	samtools view -H ${tbam} | sed 's/SM:[^[:space:]]*/SM:${meta.patient}/' > tumor_header.sam
	samtools view -H ${nbam} | sed 's/SM:[^[:space:]]*/SM:NORMAL/' > normal_header.sam
	
	samtools reheader tumor_header.sam tumor_modified.bam > tumor_final.bam
	samtools reheader normal_header.sam normal_modified.bam > normal_final.bam
	
	echo -e "${meta.patient}\\ttumor\\nNORMAL\\tcontrol" > sample_file.tsv

	samtools index tumor_final.bam
	samtools index normal_final.bam

	delly \\
	    call \\
	    ${args} \\
	    ${bcf_output} \\
	    --genome ${fasta} \\
	    ${exclude} \\
	    tumor_final.bam \\
	    normal_final.bam

	delly filter \\
	    -t \\
	    -v 3 \\
	    -m 50 \\
	    -a 0.01 \\
	    -f somatic \\
	    -s sample_file.tsv \\
	    ${bcf_filter} \\
	    ${prefix}.bcf

	bcftools convert -O v -o ${prefix}.delly.vcf ${prefix}.filtered.bcf

	awk 'BEGIN {FS=OFS=\"\\\\t\"}  /^#/ {print}' ${prefix}.delly.vcf > ${prefix}.delly.unfiltered.vcf
	awk 'BEGIN {FS=OFS=\"\\\\t\"}  \$1 ~ /^(1?[0-9]|2[0-2]|X|Y)\$/ {print}' ${prefix}.delly.vcf >> ${prefix}.delly.unfiltered.vcf

	bgzip ${prefix}.delly.vcf

	cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	    delly: \$( delly --version 2>&1 | head -1 | sed 's/Delly version: v//' )
	END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.patient}"
    """
	touch ${prefix}.delly.vcf.gz
	touch ${prefix}.delly.unfiltered.vcf

	cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	    delly: \$( delly --version 2>&1 | head -1 | sed 's/Delly version: v//' )
	END_VERSIONS
    """
}
