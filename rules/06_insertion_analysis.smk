#!/bin/env snakemake -s

polyidus_path = 'polyidus/src/polyidus.py'

#rule run_polyidus:
#    input: 
#        e1 = join(TMP, '{libtype}_merge','{sample}_{libtype}.end1.fq.gz'),
#        e2 = join(TMP, '{libtype}_merge','{sample}_{libtype}.end2.fq.gz')
#    output: directory(join(OUT,'polyidus','{sample}_{libtype}'))
#    params:
#        polyidus_path = polyidus_path,
#        human = join(GENOME_HUMAN_INDEX,'hg19'), 
#        hbv = join(GENOME_HBV_INDEX,'hbv')
#    shell:
#        """
#        mkdir {output}
#        echo "This is '{output}'"
#        python {params.polyidus_path} {params.human} {params.hbv} \
#            --fastq {input.e1} {input.e2} \
#            --outdir {output}
#        """
 #Automatically get coordinates of HBV insertion sites using Hi-C contacts

rule get_insertion_coords:
	input: join(OUT, 'cov_virus_{sample}.bedgraph')
	output: join(OUT, 'insertions', 'insertions_{sample}.bed')
	params:
		height_threshold = 200,
		distance_threshold = 22 # To get rid of the same integration event in 2 consecutive bins (1 bin= 10000 bp)
	run:
		from scipy.signal import find_peaks
		cov = pd.read_csv(
			input[0],
			sep='\t',
		   	header=None,
		  	names=['chrom', 'start', 'end', 'virus'],
			)
		peaks_arr = []
		for (i,chrom) in enumerate(np.unique(cov.chrom)):
			print(i, chrom)
			cov_sub = cov.loc[cov.chrom == chrom]
			peaks_arr.append(find_peaks(cov_sub["virus"], height = params["height_threshold"], distance = params["distance_threshold"]))
			print(peaks_arr[i][0])
			cov_peak = cov_sub.iloc[peaks_arr[i][0],0:3]
## Do not consider peaks on the virus itself
			cov_peak = cov_peak.loc[~cov_peak.chrom.str.contains("HBV"), :]
			cov_peak.to_csv(output[0], mode='a',sep="\t",header=None,index=False)
 
 
 ## Retrieve Hi-C windows and annotations around insertions
 #rule get_insertion_regions:
 #    input: join(OUT, 'insertions', 'insertions_{sample}.bed')
 #    output: join()
 #
 #
 ## Compare differential expression at HBV insertion sites versus rest
 ## of the genome (Ideally Transloc + insertion vs transloc)
 #rule compare_expr_insertion:
 #    input: join(OUT, 'insertions', 'insertions_{sample}.bed')
 #    output: join()
