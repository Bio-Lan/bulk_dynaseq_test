process CONVERSION_SUMMARY{
    tag "$meta.id"
    label 'process_single'

    conda 'conda-forge::pandas==2.2.1'
    container "biocontainers/pandas:2.2.1"

    input:
    tuple val(meta), path(postag,stageAs: "?/*"), path(snp,stageAs: "?/*")

    output:
    tuple val(meta), path("${meta.id}.conversion.csv"), emit:conv_sample

    script:
    def prefix = "${meta.id}"
    def postag = postag.join(",")
    def snp = snp.join(",")
    """
    conversion_summary.py \\
        --sample $prefix \\
        --postag $postag \\
        --snp $snp
    """
}