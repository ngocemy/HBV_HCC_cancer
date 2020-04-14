!# usr/bin/env python3
import pysam as ps
import pandas as pd
import numpy as np

bamfile = ps.AlignmentFile(snakemake.input["bam"], "rb")

hetero_read = set(open(snakemake.input["hetero"].read().split('\n'))

with open(snakemake.input["bed"],"r") as bed: 
    for line in bed:
        cov_array = np.zeros(10000)
        for read in bamfile.fetch(line[0],line[1],line[2]):
            if read.query_name in hetero_reads:
                ref_arr_mapped = read.get_reference_position() # Give an array of positions in reference genome with alignment
                cov_array[ref_arr_mapped] += 1
        site_inte = np.amax(cov_array)
        
bamfile.close()
