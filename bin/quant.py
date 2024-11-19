#!/bin/env python
# coding=utf8

import os
import sys
import gzip
import argparse

import pysam
import scipy
import pandas as pd
from multiprocessing import Pool

import utils

from __init__ import ASSAY, BARCODE_FILE_NAME, FEATURE_FILE_NAME, MATRIX_FILE_NAME

BULK_DYNASEQ_MATRIX_DIR_SUFFIX = ["raw", "labeled", "unlabeled"]

class Quant:
    """
    Features
    - Quantify total RNA, unlabeled and labeled RNA.

    Output
    - `{sample}_{RNA}` The expression matrix of total/labeled/unlabeled in Matrix Market Exchange Formats.

    """

    def __init__(self, args):
        self.args = args
        #input
        self.sample = args.sample
        self.matrix_dir = args.matrix_dir
        bam = args.conv_bam.split(",")
        snp = args.conv_snp.split(",")
        self.wells = self.get_wells(args.conv_sample)

        #set
        barcodes_file = os.path.join(self.matrix_dir,BULK_DYNASEQ_MATRIX_DIR_SUFFIX[0], BARCODE_FILE_NAME)
        features_file = os.path.join(self.matrix_dir,BULK_DYNASEQ_MATRIX_DIR_SUFFIX[0], FEATURE_FILE_NAME)
        self.features = pd.read_csv(features_file, sep="\t", header=None, index_col=0)
        self.barcodes = utils.read_one_col(barcodes_file)
        self.totaldf = pd.DataFrame()

        self.bam_snp_dict = {}
        for i in self.wells:
            temp1 = [x for x in snp if i in x]
            temp2 = [y for y in bam if i in y]
            self.bam_snp_dict[i] = {"snp":temp1[0],"bam":temp2[0]}

        #output
        self.dir_labeled = f"{self.matrix_dir}/{BULK_DYNASEQ_MATRIX_DIR_SUFFIX[1]}"
        self.dir_unlabeled = f"{self.matrix_dir}/{BULK_DYNASEQ_MATRIX_DIR_SUFFIX[2]}"
        self.detail_txt = f"{self.sample}.labeled_detail.txt"
        self.rawcsv = f'{self.sample}_raw.csv'
        self.fltcsv = f'{self.sample}_filtered.csv'

    def get_wells(self,conv_file):
        df = pd.read_csv(conv_file, index_col=0)
        return df.index.to_list()
    
    def run(self):
        well_dfs = self.run_quant()
        for i in well_dfs:
            self.totaldf = pd.concat([self.totaldf,i], axis=0)
        self.totaldf.to_csv(self.detail_txt, sep="\t", index=False)

        self.split_matrix()
        self.write_csv_file()

    def run_quant(self):
        if self.args.snp_matchfile or self.args.snp_file:
            self.parse_snp_file()
        
        in_bam_list, snp_list = [], []
        for i in self.wells:
            in_bam_list.append(f"{self.bam_snp_dict[i]['bam']}")
            snp_list.append(f"{self.bam_snp_dict[i]['snp']}")
        
        mincpu = min(len(self.wells),self.args.thread)
        with Pool(mincpu) as pool:
            results = pool.starmap(Quant.quant,zip(in_bam_list,snp_list,self.wells)) 
        return results

    def parse_snp_file(self):
        if self.args.snp_matchfile:
            snp_matchfile = self.args.snp_matchfile
            if not os.path.exists(snp_matchfile):
                sys.exit('[Error] "snp_matchfile" not exists, please check the file.\n')
            else:
                with open(snp_matchfile) as f:
                    for i in f:
                        temp = i.strip().split(',')
                        if len(temp) == 1:
                            sys.exit(f"[Error] {temp[0]} does not have matched snp file. If not, please modify it in {snp_matchfile} !!\n")
                        if not os.path.exists(temp[1]):
                            sys.exit(f"[Error] snp file not exit, please check\n{temp[0]}")
                        if temp[0] in self.wells:
                            self.bam_snp_dict[temp[0]]['snp'] = self.bam_snp_dict[j]['snp'] + ',' + temp[1]
                        else:
                            sys.exit(f"[Error] {temp[0]} is not in well list, please check the well BC in {snp_matchfile} !!\nWell BC list can be found in {self.sample}.conversion.csv")
        else:
            snp_list = self.args.snp_file.strip().split(",")
            for i in snp_list:
                if not os.path.exists(i):
                    sys.exit(f"[Error] {i} not exists, please check the file.\n")
            snp_files = ",".join(snp_list)
            for j in self.wells:
                self.bam_snp_dict[j]['snp'] = self.bam_snp_dict[j]['snp'] + ',' + snp_files
    
    @staticmethod
    def quant(bam_file,snp_file,sample):
        # get backgroud snp
        bg = Quant.background_snp(snp_file)
        # get reads with TC
        well_df = Quant.extract_dem(bam_file, bg, sample)
        return well_df
    
    @staticmethod
    def extract_dem(bam, bg, sample):
        save = pysam.set_verbosity(0)
        bamfile = pysam.AlignmentFile(bam, 'rb')
        pysam.set_verbosity(save)
        readdict = {}

        for read in bamfile.fetch():
            cb = read.get_tag("CB")
            chro = read.reference_name
            ub = read.get_tag('UB')
            gene = read.get_tag('GX')
            tctag = 0
            true_tc = []

            if read.get_tag("ST") == "+":
                stag = read.get_tag("TL")
            else:
                stag = read.get_tag("AL")
            
            if stag == '-':
                tctag = 0
                true_tc = stag
            else:
                for si in range(0, len(stag)):
                    pos = chro + '_' + str(stag[si])
                    if pos not in bg.keys():
                        true_tc.append(int(stag[si]))
                tctag = len(true_tc)

            readid = ":".join([cb, ub, gene])
            if readid not in readdict:
                readdict[readid] = tctag
            else:
                if tctag > readdict[readid]:
                    readdict[readid] = tctag
        bamfile.close()
        
        tc_df = pd.DataFrame.from_dict(readdict, orient="index", columns=["TC"])
        tc_df = tc_df.reset_index()
        tc_df[['Barcode', 'UMI',"geneID"]] = tc_df['index'].str.split(':', expand=True)
        tc_df = tc_df[['Barcode', 'UMI',"geneID","TC"]]
        return tc_df

    @staticmethod
    def background_snp(bgfiles):
        ## dict.update()
        outdict = {}
        bgs=bgfiles.strip().split(',')

        for bgfile in bgs:
            if bgfile.endswith('.csv'):
                df = pd.read_csv(bgfile, dtype={"chrom":str})
                if df.shape[0] == 0:
                    continue
                df['chrpos'] = df['chrom'].astype(str) + '_' + df['pos'].astype(str)
                df1 = df[['chrpos','posratio']]
                df1.set_index('chrpos',inplace=True)
                if len(df1.index) > 0:
                    outdict.update(df1.to_dict(orient='index'))
            else:
                sys.exit('[Error] snp file must have chrom and pos columns and should be comma-separated.Only csv format is allowed.')
        return outdict

    def split_matrix(self):
        self.labeled = self.totaldf[self.totaldf["TC"] > 0]
        unlabeled = self.totaldf[self.totaldf["TC"] == 0]
        # matrix
        self.write_sparse_matrix(self.labeled.drop("TC", axis=1), self.dir_labeled)
        self.write_sparse_matrix(unlabeled.drop("TC", axis=1), self.dir_unlabeled)
    
    def write_sparse_matrix(self,df,matrix_dir):
        count_matrix = self.dataframe_to_matrix(df, features=self.features.index, barcodes=self.barcodes)
        self.to_matrix_dir(count_matrix, matrix_dir)
    
    def dataframe_to_matrix(self, df, features, barcodes, barcode_column="Barcode", feature_column="geneID", value="UMI"):
        if df.shape[0] > 0:
            series_grouped = df.groupby([barcode_column, feature_column], observed=True).size()
            series_grouped.name = value
            df_grouped = pd.DataFrame(series_grouped)
        else:
            empty_matrix = scipy.sparse.coo_matrix((len(features), len(barcodes)))
            return empty_matrix

        feature_index_dict = {}
        for index, gene_id in enumerate(features):
            feature_index_dict[gene_id] = index
        barcode_index_dict = {}
        for index, barcode in enumerate(barcodes):
            barcode_index_dict[barcode] = index

        # use all barcodes
        barcode_codes = [barcode_index_dict[barcode] for barcode in df_grouped.index.get_level_values(level=0)]
        # use all gene_id from features even if it is not in df
        gene_id_codes = [feature_index_dict[gene_id] for gene_id in df_grouped.index.get_level_values(level=1)]
        mtx = scipy.sparse.coo_matrix(
            (df_grouped[value], (gene_id_codes, barcode_codes)), shape=(len(features), len(barcodes))
        )

        return mtx
    
    def to_matrix_dir(self, count_matrix, matrix_dir):
        utils.check_mkdir(dir_name=matrix_dir)
        self.features.to_csv(f"{matrix_dir}/{FEATURE_FILE_NAME}", sep='\t', header=False)
        pd.Series(self.barcodes).to_csv(f"{matrix_dir}/{BARCODE_FILE_NAME}", index=False, sep='\t', header=False)
        matrix_path = f"{matrix_dir}/{MATRIX_FILE_NAME}"
        with gzip.open(matrix_path, 'wb') as f:
            scipy.io.mmwrite(f, count_matrix)
        
    def write_csv_file(self):
        df_sum = self.totaldf.groupby('Barcode').agg({
            'UMI': 'count',
            'geneID': 'nunique'
        })
        df_labeled_sum = self.labeled.groupby('Barcode').agg({
            'UMI': 'count',
            'geneID': 'nunique'
        })
        df = pd.concat([df_sum,df_labeled_sum], axis=1)
        df.columns=['UMI','gene','labeled_UMI','labeled_gene']
        df = df.fillna(0)
        df = df.astype(int)
        df.to_csv(self.rawcsv)

        df_filter = df[ (df['UMI']>=self.args.umi_cutoff) & (df['gene']>=self.args.gene_cutoff) ]
        df_filter.to_csv(self.fltcsv)

        json_df = df_filter if df_filter.shape[0] > 0 else df
        self.write_multiqc_json(json_df)
    
    def write_multiqc_json(self,df):
        json_dict = df.to_dict(orient="index")
        utils.write_multiqc(json_dict, self.sample, ASSAY, "quant.well_inf")
        
        stats = df.describe()
        stats.columns = ['UMI','Genes','labeled_UMI','labeled_gene']

        df['UMI_rate'] = df['labeled_UMI']/df['UMI']
        df['gene_rate'] = df['labeled_gene']/df['gene']
        labeled_rate = round( df['UMI_rate'].median() * 100,2)
        labeled_rate2 = round( df['UMI_rate'].mean() * 100,2)

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
        utils.write_multiqc(data_dict, self.sample, ASSAY, "quant.stats")
        
        df = df.loc[:,['UMI_rate','gene_rate']]
        df.columns = ['UMI','Gene']
        label_dict = df.to_dict(orient="list")
        utils.write_multiqc(label_dict, self.sample, ASSAY, "quant.labeled_rate")

if __name__ == "__main__":
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample",required=True)
    parser.add_argument("--matrix_dir",required=True)
    parser.add_argument("--conv_sample",required=True)
    parser.add_argument('--conv_bam',required=True)
    parser.add_argument('--conv_snp',required=True)
    parser.add_argument('--umi_cutoff', default=500, type=int,
        help='If the UMI number exceeds the threshold, it is considered a valid well and reported.'
    )
    parser.add_argument('--gene_cutoff', default=0, type=int,
        help='If the gene number exceeds the threshold, it is considered a valid well and reported.'
    )
    parser.add_argument('--thread', type=int)   
    parser.add_argument('--snp_file', type=str,
                        help="""backgroud snp file.
                            Can be multiple files, separated by commas. 
                            If this option is set, it is valid for all wells""")
    parser.add_argument('--snp_matchfile', type=str,
                        help="""backgroud snp file for each well, one well per line, the format is  \"well,snp_file\". 
                            This parameter takes precedence over snp_file. """)              
    args = parser.parse_args()

    runner = Quant(args)
    runner.run()