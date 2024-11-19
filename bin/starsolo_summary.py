#!/usr/bin/env python

import argparse
from collections import defaultdict

import pandas as pd

import utils
from __init__ import ASSAY

class Starsolo_summary:
    def __init__(self, args):
        self.args = args
        self.stats = {}

    def parse_read_stats(self):
        dtypes = defaultdict(lambda: "int")
        dtypes["CB"] = "object"
        df = pd.read_csv(
            self.args.read_stats, sep="\t", header=0, index_col=0, skiprows=[1], dtype=dtypes
        )  # skip first line cb not pass whitelist
        df = df.loc[
            :,
            [
                "cbMatch",
                "cbPerfect",
                "genomeU",
                "genomeM",
                "exonic",
                "intronic",
                "exonicAS",
                "intronicAS",
                "countedU",
                "nUMIunique",
                "nGenesUnique",
            ],
        ]
        s = df.sum()
        # json does not recognize NumPy data types. TypeError: Object of type int64 is not JSON serializable
        valid = int(s["cbMatch"])
        perfect = int(s["cbPerfect"])
        corrected = valid - perfect
        genome_uniq = int(s["genomeU"])
        genome_multi = int(s["genomeM"])
        mapped = genome_uniq + genome_multi
        exonic = int(s["exonic"])
        intronic = int(s["intronic"])
        antisense = int(s["exonicAS"] + s["intronicAS"])
        intergenic = mapped - exonic - intronic - antisense
        counted_uniq = int(s["countedU"])
        data_dict = {
            "Corrected Barcodes": corrected / valid,
            "Reads Mapped To Unique Loci": genome_uniq / valid,
            "Reads Mapped To Multiple Loci": genome_multi / valid,
            "Reads Mapped Uniquely To Transcriptome": counted_uniq / valid,
            "Mapped Reads Assigned To Exonic Regions": exonic / mapped,
            "Mapped Reads Assigned To Intronic Regions": intronic / mapped,
            "Mapped Reads Assigned To Intergenic Regions": intergenic / mapped,
            "Mapped Reads Assigned Antisense To Gene": antisense / mapped,
        }
        for k in data_dict:
            data_dict[k] = utils.get_frac(data_dict[k])
        self.stats.update(data_dict)

    def parse_summary(self):
        data = utils.csv2dict(self.args.summary)
        origin_new = {
            "Number of Reads": "Raw Reads",
            "Reads With Valid Barcodes": "Valid Reads",
            "Q30 Bases in CB+UMI": "Q30 Bases in CB+UMI",
            "Q30 Bases in RNA read": "Q30 Bases in RNA read"
        }
        parsed_data = {}
        for origin, new in origin_new.items():
            parsed_data[new] = data[origin]
        frac_names = {"Valid Reads","Q30 Bases in CB+UMI","Q30 Bases in RNA read"}
        for k in frac_names:
            parsed_data[k] = utils.get_frac(parsed_data[k])
        for k in set(origin_new.values()) - frac_names:
            parsed_data[k] = int(parsed_data[k])
        self.stats.update(parsed_data)

    def run(self):
        self.parse_read_stats()
        self.parse_summary()
        utils.write_multiqc(self.stats, self.args.sample, ASSAY, "starsolo_summary.stats")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Starsolo summary")
    parser.add_argument("--sample", help="sample name")
    parser.add_argument("--read_stats", help="cellReadsStats file")
    parser.add_argument("--summary", help="summary file")
    args = parser.parse_args()

    Starsolo_summary(args).run()