process CONVERSION{
    tag "$meta.id"
    label 'process_high'

    conda 'conda-forge::pandas==2.2.1 bioconda::pysam==0.22.1 bioconda::samtools==1.20'
    container "qaqlans/sgr-accura-2"

    input:
    tuple val(meta), path(well_bam)

    output:
    tuple val(meta), path("${meta.id}/*.PosTag.bam"), emit:conv_bam
    tuple val(meta), path("${meta.id}/*.PosTag.bam.bai"), emit:conv_bam_bai
    tuple val(meta), path("${meta.id}/*.PosTag.csv"), emit:conv_postag
    tuple val(meta), path("${meta.id}/*.snp.csv"), emit:conv_snp

    script:
    def prefix = "${meta.id}"
    """
    conversion.py \\
        --sample $prefix \\
        --wellBAM $well_bam \\
        --gtf ${params.gtf} \\
        --conversion_type ${params.conversion_type} \\
        --basequalilty ${params.basequalilty} \\
        --snp_threshold ${params.snp_threshold} \\
        --snp_min_depth ${params.snp_min_depth}
    """
}