process SUBSTITUTION{
    tag "$meta.id"
    label 'process_medium'
    
    conda 'conda-forge::pandas==2.2.1 bioconda::pysam==0.22.1 bioconda::samtools==1.20'
    container "qaqlans/sgr-accura-2"
    
    input:
    tuple val(meta), path(conv_sample), path(conv_bam,stageAs: "?/*")

    output:
    tuple val(meta), path("*.substitution.csv"),emit:sample_subs

    script:
    def prefix = "${meta.id}"
    def conv_bam = conv_bam.join(",")
    """
    substitution.py \\
        --sample $prefix \\
        --conv_sample $conv_sample \\
        --conv_bam $conv_bam
    """
}