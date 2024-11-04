#!/usr/bin/env python

import argparse
import pandas as pd

if __name__ == "__main__":
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample",required=True)
    parser.add_argument('--postag',required=True)
    parser.add_argument("--snp",required=True)
    args = parser.parse_args()

    postag = args.postag.split(",")
    snp = args.snp.split(",")
    conversion_num = {}
    outfile_csv = args.sample+'.conversion.csv'
    
    for i in postag:
        df_postag = pd.read_csv(i,header=0,sep=",",dtype={"chrom":"str"})
        if df_postag.shape[0] == 0:
            continue

        wellBC = i.split("/")[1].split(".")[0]
        df_snp = [x for x in snp if wellBC in x]
        df_snp = pd.read_csv(df_snp[0],header=0,sep=",",dtype={"chrom":"str"})
        conversion_num[wellBC] = [df_postag.shape[0], df_snp.shape[0]]
    
    df = pd.DataFrame.from_dict(conversion_num,orient='index')
    df.columns = ['conv_num','snp_num']
    df.to_csv(outfile_csv)

