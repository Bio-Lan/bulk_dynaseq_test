process RESULT_SUMMARY{
    tag "$meta.id"
    label 'process_high'
    
    conda 'conda-forge::pandas==2.2.1'
    container "biocontainers/pandas:2.2.1"

    input:
    tuple val(meta), path(read_stats), path(summary), path(sample_filter), path(sample_raw), path(sample_subs)

    output:
    tuple val(meta), path("*.json"), emit:json

    script:
    def prefix = "${meta.id}"
    """
    result_summary.py \\
        --sample $prefix \\
        --read_stats $read_stats \\
        --summary $summary \\
        --filter_csv $sample_filter \\
        --raw_csv $sample_raw \\
        --substitution $sample_subs
    """
}