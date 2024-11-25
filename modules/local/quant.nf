process QUANT{
    tag "$meta.id"
    label 'process_medium'
    
    conda 'conda-forge::pandas==2.2.1 bioconda::pysam==0.22.1 conda-forge::scipy==1.14.0'
    container "qaqlans/sgr-accura-3"
    
    input:
    tuple val(meta), path(conv_sample), path(conv_bam,stageAs: "./well_file/*"), path(conv_bam_bai,stageAs: "./well_file/*"), path(conv_snp,stageAs: "./well_file/*"), path(raw_matrix)

    output:
    tuple val(meta), path("${meta.id}.matrix/"), emit:sample_matrix
    tuple val(meta), path("${meta.id}_raw.csv"), emit:sample_raw
    tuple val(meta), path("${meta.id}_filtered.csv"), emit:sample_filter
    tuple val(meta), path("${meta.id}.labeled_detail.txt")
    
    script:
    def prefix = "${meta.id}"
    def matrix_dir = "${prefix}.matrix"
    def conv_bam = conv_bam.join(",")
    def conv_snp = conv_snp.join(",")
    def args = task.ext.args ?: ''
    """
    mkdir ${matrix_dir}
    cp -r -L ${raw_matrix} ${matrix_dir}/raw
    gzip ${matrix_dir}/raw/*

    quant.py \\
        --sample $prefix \\
        --matrix_dir $matrix_dir \\
        --conv_sample $conv_sample \\
        --conv_bam $conv_bam \\
        --conv_snp $conv_snp \\
        --umi_cutoff ${params.umi_cutoff} \\
        --gene_cutoff ${params.gene_cutoff} \\
        --thread ${task.cpus} \\
        $args
    """
}