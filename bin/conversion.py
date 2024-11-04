#!/usr/bin/env python
import os
import argparse

import pysam
import pandas as pd

import utils
import filter_gtf
import func_conversion

class Conversion:
    """
    ## Features
    - Get conversion pos in each read.
        - Get snp info. 

    ## Output
    - `{sample}.PosTag.bam` Bam file with conversion info.
    - `{sample}.PosTag.csv` TC conversion sites info in csv format.
    - `{sample}.snp.csv` Candidated snp sites.
    """

    def __init__(self, args):
        #input
        self.args = args
        self.inbam = args.wellBAM
        gp = filter_gtf.GtfParser(self.args.gtf)
        self.strand = gp.get_id_strand()

        #output
        os.makedirs(args.sample, exist_ok=True)
        self.wellBC = os.path.splitext(os.path.basename(self.inbam))[0]
        output_prefix = args.sample + "/"+ self.wellBC
        self.outbam = output_prefix + '.PosTag.bam'
        self.outPostag = output_prefix + '.PosTag.csv'
        self.outsnp = output_prefix + '.snp.csv'
    
    def run(self):
        # Adding tags and parse snps
        self.depth_dict = func_conversion.addTags(self.inbam,self.outbam,self.args.conversion_type,self.args.basequalilty,self.strand)
        self.loci_process()

    def loci_process(self):
        utils.index_bam(self.outbam)
        if len(self.depth_dict) == 0:
            df = pd.DataFrame(columns=['chrom', 'pos', 'convs','covers','posratio'])
            df.to_csv(self.outPostag)
            df.to_csv(self.outsnp)
        else:
            loci_dict = {}

            save = pysam.set_verbosity(0)
            bam = pysam.AlignmentFile(self.outbam, 'rb')
            pysam.set_verbosity(save)
            for key in self.depth_dict.keys():
                chrom_pos = key.split("+")                
                ConvsDepth = self.depth_dict[key] #conversion pos depth

                # read with conversion loci count
                CoverofPosWithConvs = bam.count(str(chrom_pos[0]), int(chrom_pos[1]), int(chrom_pos[1])+1)
                posratio = int(ConvsDepth)/int(CoverofPosWithConvs)
                temp_dict = {"chrom":str(chrom_pos[0]),"pos":int(chrom_pos[1]),
                             "convs":ConvsDepth,"covers":CoverofPosWithConvs,
                             "posratio":posratio}
                loci_dict[key] = temp_dict
            bam.close()
            df = pd.DataFrame.from_dict(loci_dict, orient='index')
            snp_df = df[ df['convs'] >= self.args.snp_min_depth ]
            snp_df = snp_df[ snp_df['posratio'] >= self.args.snp_threshold]
            
            df.to_csv(self.outPostag,index=False)
            snp_df.to_csv(self.outsnp,index=False)

if __name__ == "__main__":
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--wellBAM", required=True,
                        help='STARsolo well BAM(sortedByCoord),must have "MD" tag')
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--conversion_type", type=str,default="TC",
                        help='conversion type, TC for dynaseq', required=False)
    parser.add_argument("--basequalilty", default=20, type=int,
                        help='min base quality of the read sequence', required=False)
    parser.add_argument("--snp_threshold", type=float, default=0.5,
                        help='snp threshold filter, greater than snp_threshold will be recognized as snp')
    parser.add_argument("--snp_min_depth",default=10, type=int,
                        help='Minimum depth to call a variant')
    # add version
    parser.add_argument("--version", action="version", version="1.0")

    args = parser.parse_args()

    runner = Conversion(args)
    runner.run()



