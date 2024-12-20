process SPLIT_FASTQ {
    tag "$meta.id"
    label 'process_medium'

    conda 'conda-forge::pandas==2.2.1 bioconda::pyfastx==2.1.0'
    container "qaqlans/sgr-accura-1"
    
    input:
    tuple val(meta), path(reads)
    path assets_dir
    val protocol
    path split_inf

    output:
    tuple val(meta), path("${meta.id}/")
    tuple val(meta), path("${meta.id}.*.json")

    script:
    def prefix = "${meta.id}"
    def (forward, reverse) = reads.collate(2).transpose()
    def args = task.ext.args ?: ''
    """
    split_fastq.py \\
        --sample $prefix \\
        --fq1 ${forward.join( "," )} \\
        --fq2 ${reverse.join( "," )} \\
        --split_inf $split_inf \\
        --assets_dir $assets_dir \\
        --protocol $protocol \\
        --pattern ${params.pattern} \\
        $args
    """
}

workflow FASTQ_SPLIT {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:
    SPLIT_FASTQ (
        ch_samplesheet,
        "${projectDir}/assets/",
        params.protocol,
        params.split_inf,
    )
}