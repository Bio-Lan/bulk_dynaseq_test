process REPORT_SUMMARY{
    tag "$meta.id"
    label 'process_single'
    
    conda 'conda-forge::pandas==2.2.1'
    container "biocontainers/pandas:2.2.1"

    input:
    tuple val(meta), path(read_stats), path(summary), path(substitution_stat), path(sample_raw), path(sample_filter)

    output:
    tuple val(meta), path("*.json"), emit:json

    script:
    def prefix = "${meta.id}"
    
    """
    report_summary.py \\
        --sample $prefix \\
        --read_stats $read_stats \\
        --summary $summary \\
        --sub_stats $substitution_stat \\
        --sample_raw $sample_raw \\
        --sample_filter $sample_filter
    """
}