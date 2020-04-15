import pysam as ps
import pandas as pd
import numpy as np
import pdb
bamfile = ps.AlignmentFile(snakemake.input["bam"], "rb")

hetero_reads = set(open(snakemake.input["hetero"]).read().split('\n'))
bed = pd.read_csv(snakemake.input["bed"],sep="\t",header=None,names=["chr","start","end"])
with open(snakemake.output[0],"w") as bedout:
    i=0
    for index,line in bed.iterrows():
        cov_array = np.zeros(10000)
        for read in bamfile.fetch(line["chr"],line["start"],line["end"]):
            if read.query_name in hetero_reads:
                ref_arr_mapped = np.array(read.get_reference_positions()) # Give an array of positions in reference genome with alignment
                ref_arr_mapped = ref_arr_mapped[ref_arr_mapped < (line["start"] + 10000)]
                cov_array[ref_arr_mapped - line["start"]] += 1
        site_inte = np.argmax(cov_array)
        if i < 10:
            plt.plot(cov_array)
            plt.savefig("{read} array.png")
        else:
            system.exit(0)
        bedout.write(''.join([line["chr"],"\t", str(line["start"] +site_inte),"\n"])) 
        i+=1
bamfile.close()
bedout.close()
