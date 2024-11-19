def multiqc_sgr_config():
    from multiqc import config

    """ Set up MultiQC config defaults for this package """
    sgr_search_patterns = {
        "bulk_dynaseq/stats": {
            "fn": "*bulk_dynaseq.*.stats.json"
        },
        "bulk_dynaseq/boxplot": {
            "fn": "*bulk_dynaseq.substitution.boxplot.json"
        },
        "bulk_dynaseq/barplot": {
            "fn": "*bulk_dynaseq.substitution.barplot.json"
        },
        "bulk_dynaseq/labeled_rate": {
            "fn": "*bulk_dynaseq.quant.labeled_rate.json"
        },
        "bulk_dynaseq/well_inf": {
            "fn": "*bulk_dynaseq.quant.well_inf.json"
        }
    }
    config.update_dict(config.sp, sgr_search_patterns)
