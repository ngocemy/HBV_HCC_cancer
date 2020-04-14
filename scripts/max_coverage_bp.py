!# usr/bin/env python3
import pysam as ps
import pandas as pd
import numpy as np

bamfile = ps.AlignmentFile(snakemake.input["bam"], "rb")

hetero_read = set(open(snakemake.input["hetero"].read().split('\n'))
bed = pd.read_csv(snakemake.input["bed"],sep="\t",header=None,names=["chr","start","end"])
with open(snakemake.output[0]) as bedout:
    for index,line in bed.iterrows():
        cov_array = np.zeros(10000)
        for read in bamfile.fetch(line["chr"],line["start"],line["end"]):
            if read.query_name in hetero_reads:
                ref_arr_mapped = read.get_reference_position() # Give an array of positions in reference genome with alignment
                cov_array[ref_arr_mapped] += 1
        site_inte = np.amax(cov_array)
        bedout.write('\t'.join(line["chr"],site_inte)) 
bamfile.close()
