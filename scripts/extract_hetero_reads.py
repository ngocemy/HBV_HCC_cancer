import pandas as pd
import numpy as np

def extract_hetero_reads(pair,out):
    """Extract hetero reads comes from both human and HBV and write to file. 
    Reads are not ensure to be chimeric, but rather HiC contacts.
    
    Parameters
    ----------
    pair: pair file has been extracted.
    out: output file with read names, eg. hetero_read.tsv
    @author: ngocusth
    """
    
    for fpair in pd.read_csv(pair,sep='\t',comment='#',encoding='utf-8',header=None, chunksize=100000):
    
        #Initiate array of 0
        # Column 2,4 are chromosomes
        col2=np.zeros(fpair.shape[0])
        col4=np.zeros(fpair.shape[0])
    
        col4[fpair.index[fpair.iloc[:,3] == "HBVayw"]] = 1
    
        col2[fpair.index[fpair.iloc[:,1] == "HBVayw"]] = 1
    
        # Keep only index with (0,1) or (1,0) in pair (col2,col4) 
        het_col = np.logical_xor(col2,col4)
    
        read_het = fpair.iloc[het_col]
    
        try:
            read_het=read_het.iloc[:,0]   
            read_het.to_csv(out, mode='a', sep='\t')
        except IndexingError:
            pass
extract_hetero_reads('hbv_data/tmp/pairs/PLC_1_captn.pairs.sorted','hetero_read.tsv')

#extract_hetero_reads(snakemake.input[0],snakemake.output[0])
