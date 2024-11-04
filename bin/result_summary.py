#!/usr/bin/env python

import argparse
from collections import defaultdict

import pandas as pd

import utils
from __init__ import ASSAY

class Result_Summary:
    def __init__(self, args):
        self.args = args
        self.stats = {}

    def parse_read_stats(self):
        dtypes = defaultdict(lambda: "int")
        dtypes["CB"] = "object"
        df = pd.read_csv(
            self.args.read_stats, sep="\t", header=0, index_col=0, skiprows=[1], dtype=dtypes
        )  # skip first line cb not pass whitelist
        rbs = list(df.index)
        umi_count = list(df["nUMIunique"])
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
    
    def parse_filterCSV(self):
        df = pd.read_csv(self.args.filter_csv,header=0,sep=",",index_col=0)
        if df.shape[0] == 0:
            df = pd.read_csv(self.args.raw_csv,header=0,sep=",",index_col=0)
        table_dict = df.to_dict(orient='index')
        utils.write_multiqc(table_dict, self.args.sample, ASSAY, "quant.table")

        
        labeled_df = df.dropna()
        labeled_df.loc[:,'rate'] = labeled_df.loc[:,'labeled_UMI'] / labeled_df.loc[:,'UMI']
        labeled_rate = round( labeled_df['rate'].median() * 100,2)
        labeled_rate2 = round( labeled_df['rate'].mean() * 100,2)

        stats = df.describe()
        stats.columns = ['UMI','Genes','labeled_UMI','labeled_gene']
        
        data_dict = {}
        for item in ['UMI', 'Genes']:
            temp = f'Median {item} across Well'
            data_dict[temp] = int(stats.loc['50%',item])
            temp = f'Mean {item} across Well'
            data_dict[temp] = int(stats.loc['mean',item])
        
        temp = 'Median labeled rate across wells'
        data_dict[temp] = labeled_rate
        temp = 'Mean labeled rate across wells'
        data_dict[temp] = labeled_rate2
        self.stats.update(data_dict)

        df['UMI_rate'] = df['labeled_UMI'] / df['UMI']
        df['gene_rate'] = df['labeled_gene'] / df['gene']
        label_df = df.loc[:,['UMI_rate','gene_rate']]
        label_df.columns = ['UMI','Gene']
        label_dict = label_df.to_dict(orient="list")
        utils.write_multiqc(label_dict, self.args.sample, ASSAY, "quant.labeled")

    def parse_substitution(self):
        file_path = self.args.substitution
        df = pd.read_csv(file_path,header=0,sep=",",index_col=0)
        box_dict = df.to_dict(orient='list')
        utils.write_multiqc(box_dict, self.args.sample, ASSAY, "substitution.boxplot")
        bar_dict = df.to_dict(orient="index")
        utils.write_multiqc(bar_dict, self.args.sample, ASSAY, "substitution.barplot")

    def run(self):
        self.parse_read_stats()
        self.parse_summary()
        self.parse_filterCSV()
        self.parse_substitution()
        utils.write_multiqc(self.stats, self.args.sample, ASSAY, "reads.stats")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Result summary")
    parser.add_argument("--sample", help="sample name")
    parser.add_argument("--read_stats", help="cellReadsStats file")
    parser.add_argument("--summary", help="summary file")
    parser.add_argument("--filter_csv", help="filtered_well information")
    parser.add_argument("--raw_csv", help="raw well information")
    parser.add_argument("--substitution", help="substitution result")
    args = parser.parse_args()

    Result_Summary(args).run()