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

from __init__ import BARCODE_FILE_NAME, FEATURE_FILE_NAME, MATRIX_FILE_NAME

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
        self.rawcsv = f'{self.sample}_raw.csv'
        self.fltcsv = f'{self.sample}_filtered.csv'

    def get_wells(self,conv_file):
        df = pd.read_csv(conv_file, index_col=0)
        return df.index.to_list()
    
    def run(self):
        self.dfs = self.run_quant()
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
        df_total = Quant.extract_dem(bam_file, bg, sample)
        return df_total
    
    @staticmethod
    def extract_dem(bam, bg, sample):
        save = pysam.set_verbosity(0)
        bamfile = pysam.AlignmentFile(bam, 'rb')
        pysam.set_verbosity(save)
        countarr = []
        for read in bamfile.fetch():
            chro = read.reference_name
            ub = read.get_tag('UB')
            gene = read.get_tag('GX')

            if read.get_tag("ST") == "+":
                stag = read.get_tag("TL")
            else:
                stag = read.get_tag("AL")
            if stag == '-':
                tctag = 0
            else:
                fcount = 0
                for si in range(0, len(stag)):
                    pos = chro + '_' + str(stag[si])
                    if pos in bg.keys():
                        fcount += 1
                if fcount == len(stag):
                    tctag = 0
                else:
                    tctag = 1
            countarr.append([gene, sample, ub, tctag])
        bamfile.close()

        totaldf = pd.DataFrame(countarr)
        totaldf.columns = ['gene','Well','UMI','TC']
        df_unique = totaldf.groupby(['gene', 'Well', 'UMI'])[['gene', 'Well', 'UMI','TC']].apply(lambda x: x.loc[x['TC'].idxmax()]).reset_index(drop=True)
        return df_unique

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
        # merge matrix
        for i in self.dfs:
            self.totaldf = pd.concat([self.totaldf,i])
        self.newdf = self.totaldf[self.totaldf['TC']==1]
        self.olddf = self.totaldf[self.totaldf['TC']==0]
        outnewdf = self.newdf.groupby(['Well', 'gene']).agg({'UMI': 'count'})
        outolddf = self.olddf.groupby(['Well', 'gene']).agg({'UMI': 'count'})
        
        self.write_sparse_matrix(outnewdf, self.dir_labeled)
        self.write_sparse_matrix(outolddf, self.dir_unlabeled)
    
    def write_sparse_matrix(self,df,matrix_dir):
        count_matrix = self.dataframe_to_matrix(df, features=self.features.index, barcodes=self.barcodes)
        self.to_matrix_dir(count_matrix, matrix_dir)
    
    def dataframe_to_matrix(self,df,features,barcodes,value='UMI'):

        feature_index_dict = {}
        for index, gene_id in enumerate(features):
            feature_index_dict[gene_id] = index
        barcode_index_dict = {}
        for index, barcode in enumerate(barcodes):
            barcode_index_dict[barcode] = index
        
        # use all barcodes
        barcode_codes = [barcode_index_dict[barcode] for barcode in df.index.get_level_values(level=0)]
        # use all gene_id from features even if it is not in df
        gene_id_codes = [feature_index_dict[gene_id] for gene_id in df.index.get_level_values(level=1)]
        mtx = scipy.sparse.coo_matrix((df[value], (gene_id_codes, barcode_codes)),
            shape=(len(features), len(barcodes)))

        return  mtx
    
    def to_matrix_dir(self, count_matrix, matrix_dir):
        utils.check_mkdir(dir_name=matrix_dir)
        self.features.to_csv(f"{matrix_dir}/{FEATURE_FILE_NAME}", sep='\t', header=False)
        pd.Series(self.barcodes).to_csv(f"{matrix_dir}/{BARCODE_FILE_NAME}", index=False, sep='\t', header=False)
        matrix_path = f"{matrix_dir}/{MATRIX_FILE_NAME}"
        with gzip.open(matrix_path, 'wb') as f:
            scipy.io.mmwrite(f, count_matrix)
        
    def write_csv_file(self):
        df_sum = self.totaldf.groupby('Well').agg({
            'UMI': 'count',
            'gene': 'nunique'
        })
        df_new_sum = self.newdf.groupby('Well').agg({
            'UMI': 'count',
            'gene': 'nunique'
        })
        self.df = pd.concat([df_sum,df_new_sum], axis=1)
        self.df.columns=['UMI','gene','labeled_UMI','labeled_gene']
        self.df.to_csv(self.rawcsv)

        self.df = self.df[ (self.df['UMI']>=self.args.umi_cutoff) & (self.df['gene']>=self.args.gene_cutoff) ]
        self.df.to_csv(self.fltcsv)       

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