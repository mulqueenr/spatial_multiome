"""
singularity shell \
--bind /volumes/USR2/Ryan/projects/spatial_wgs/tools/spatial_multiome/src:/src/,/volumes/USR2/Ryan/tools/cellranger-atac-2.1.0/:/cellranger/,/volumes/USR2/Ryan/projects/spatial_wgs/data/250129_First_Experiment/DNA_SampleSheet.csv:/samplesheet.tsv \
    /volumes/USR2/Ryan/projects/10x_MET/src/amethyst.sif
"""

import os
from Bio.Seq import Seq
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--whitelist',default="/cellranger/lib/python/atac/barcodes/737K-cratac-v1.txt.gz")
parser.add_argument('--samplesheet',default="/samplesheet.tsv")
parser.add_argument('--gem_idx',default="initial_gem_idx.txt")
parser.add_argument('--prefix',default="spatialdna")
parser.add_argument('--gem_cutoff',default=5000)
parser.add_argument('--outdir',default="./")
parser.add_argument('--sequencing_cycles',default="Y151;I8N2;I24;Y151")
args = parser.parse_args()

os.system("mkdir -p "+ args.outdir)  

#read gem counts
gem_idx=pd.read_csv(args.gem_idx,names=['counts',"idx"],sep="\t",dtype={'counts':'Int64','idx':str})
gem_idx['corrected_idx'] = [i[9:] for i in gem_idx['idx'].to_list()]

# read whitelist
whitelist=pd.read_csv(args.whitelist,names=["idx"]) #Read whitelist
whitelist['revcomp']=[str(Seq(i).reverse_complement()) for i in whitelist['idx']]

# read in samplesheet to extract i7 indexes
i7_idx=pd.read_csv(args.samplesheet,skiprows=4)

#subset gem counts by whitelist matches (assuming highest numbers have less misscalls)
gem_idx=gem_idx[gem_idx['corrected_idx'].isin(whitelist['revcomp'])]
gem_idx = gem_idx.sort_values(by='counts',ascending=False)
gem_idx_pass=pd.DataFrame()
gem_idx_pass['index2']=gem_idx.iloc[0:int(args.gem_cutoff)+1]['corrected_idx']
gem_idx_pass['index2']=[str(Seq(i).reverse_complement()) for i in gem_idx_pass['index2']]
gem_idx_pass['Sample_ID']=[args.prefix+"_"+i for i in gem_idx_pass['index2']]
gem_idx_pass=pd.concat([gem_idx_pass] * len(i7_idx['index']))

index_list=[[i7] * gem_idx_pass.shape[0] for i7 in i7_idx['index']]
gem_idx_pass['index'] = [element for innerList in index_list for element in innerList]

pd.concat(gem_idx_pass,gem_idx_pass)
gem_i7_df=[gem_idx_pass['index']=[i7 for i7 in i7_idx["index"].to_list()] #replicate column per i7 index

df = pd.concat(gem_i7_df) #concat list
df.columns=['Sample_ID','index2','index'] #set names
gem_idx_pass= df.col[['Sample_ID','index','index2']]

with open("samplesheet_gemidx.csv", "w") as f:
    f.write("""[Settings],
CreateFastqForIndexReads,0
OverrideCycles,"""+args.sequencing_cycles+"""
BarcodeMismatchesIndex1,2
BarcodeMismatchesIndex2,0
[Data],
Sample_ID,index,index2
""")

gem_idx_pass.to_csv('samplesheet_gemidx.csv',mode='a',sep=",",header=False,index=False)