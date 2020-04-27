import matplotlib as mpl
mpl.use("Agg")
import pysam as ps
import pandas as pd
import numpy as np
import pdb
import matplotlib.pyplot as plt
#Read in a list of BAM paths with both end1 and end2
bamfiles = [ps.AlignmentFile(b, "rb") for b in snakemake.input["bams"]]
#Read in list of read names that have human-virus contacts
hetero_reads = set(open(snakemake.input["hetero"]).read().split('\n'))
# Read in BED file contain the regions of 10KB that containing integration sites
bed = pd.read_csv(snakemake.input["bed"],sep="\t",header=None,names=["chr","start","end"])
# Give output as file with integration sites and count, plot coverage for those regions 
with open(snakemake.output["txt_out"],"w") as bedout:
    plt.figure()
    fig, axs = plt.subplots(bed.shape[0])
    i = 0
    for index,line in bed.iterrows():
        cov_array = np.zeros(10000)
        for bamfile in bamfiles:
            for read in bamfile.fetch(line["chr"],line["start"],line["end"]):
                if read.query_name in hetero_reads:
                    ref_arr_mapped = np.array(read.get_reference_positions()) # Give an array of positions in reference genome with alignment
                    ref_arr_mapped = ref_arr_mapped[ref_arr_mapped < (line["start"] + 10000)]
                    cov_array[ref_arr_mapped - line["start"]] += 1
        site_inte = np.argmax(cov_array)
        value_max = np.amax(cov_array)
        bedout.write(''.join([line["chr"],"\t", str(line["start"] +site_inte),"\t",str(value_max),"\n"])) 
        #Plot cov_array to see noise for each region - i:

        axs[i].plot(cov_array)
        axs[i].scatter(site_inte, cov_array[site_inte])
        axs[i].set_title(f"{line['chr']}: {line['start']} - {line['end']}")
        i+=1
   # Hide x labels and tick labels for all but bottom plot.
    for ax in axs:
        ax.label_outer() 
        ax.set(ylabel='Number of reads mapped per bp')
    fig.suptitle(f"Coverage of integration events in regions of 10KB: {snakemake.wildcards['sample']}")
    plt.savefig(snakemake.output["svg"])
