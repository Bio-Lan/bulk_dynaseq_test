process STARSOLO_SUMMARY{
    tag "$meta.id"
    label 'process_high'
    
    conda 'conda-forge::pandas==2.2.1'
    container "biocontainers/pandas:2.2.1"

    input:
    tuple val(meta), path(read_stats), path(summary)

    output:
    tuple val(meta), path("*.json"), emit:json

    script:
    def prefix = "${meta.id}"
    """
    starsolo_summary.py \\
        --sample $prefix \\
        --read_stats $read_stats \\
        --summary $summary
    """
}