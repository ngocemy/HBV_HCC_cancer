# Snakemake workflow to handle the border calling via cooltools.
# @ngocemy, 22052020; adapted from loop calling @cmdoret
conda: "envs/hic_processing.yaml"
L_TMP = join(TMP, "borders")
L_OUT = join(OUT, "borders")

# Merge all Hi-C matrices to obtain union of all border coordinates
rule merge_matrices:
    input: 
        expand(
            join(OUT, 'cool', '{sample}_hic_' + f'{MAX_RES_STR}.cool'),
            sample=hic_libs
        )
    output: join(OUT, 'cool', 'all_merged.cool')
    run: 
        mats = ' '.join([cool for cool in input[:]])
        shell("cooler merge {output} " + f"{mats}")


# Perform border detection on the merged matrices
# chromosight version: 
rule loop_calling:
    input: 
        cool = join(OUT, 'cool', 'all_merged.cool')
    output:
        coords = join(L_OUT, 'merged_contacts', 'borders_out.txt'),
        wins = join(L_OUT, 'merged_contacts', 'borders_out.npy')
    threads: 1
    shell:
        """ 
        mkdir -p $(dirname {output.coords})
        chromosight detect \
            --pattern borders \
            --threads {threads} \
            --max-dist 20000000 \
            --min-dist 20000 \
            --min-separation 40000 \
            --perc-undetected 90 \
            --iterations 1 \
            --win-fmt npy \
            {input.cool}
            $(dirname {output.coords})
        """

# Matrices must be subsampled to the same coverage for border scores to be 
# comparable. Find the lowest number of contacts among all libraries
rule find_subsampling_value:
    input:
        expand(
            join(OUT, 'cool', '{sample}_hic_' + f'{MAX_RES_STR}.cool'),
            sample=hic_libs
        )
    output: join(L_OUT, 'target_contacts.txt')
    run:
        min_contacts = np.inf
        for inp in input[:]:
            c = cooler.Cooler(inp)
            c_sum = c.info['sum']
            if c_sum < min_contacts:
                min_contacts = c_sum
            print(f"{inp}: {c_sum}")
        print(f"lowest : {min_contacts}")
        with open(output[0], 'w') as out:
            out.write(str(min_contacts))


# Quantify border scores for all detected positions on individual matrices.
rule quantify_border_scores:
    input:
        cool = join(OUT, 'cool', '{sample}_hic_' + f'{MAX_RES_STR}.cool'),
        coords = join(L_OUT, 'merged_contacts', 'borders_out.txt'),
        subsample = join(L_OUT, 'target_contacts.txt')
    output:
        coords = join(L_OUT, '{sample}', 'borders_quant.txt'),
        wins = join(L_OUT, '{sample}', 'borders_quant.npy')
    shell:
        """
        chromosight quantify {input.coords} \
                             {input.cool} \
                             --win-fmt npy \
                             --subsample $(cat {input.subsample}) \
                             $(dirname {output.coords})
        """


# 02 Aggregate borders into pileups for different bins of genomic distances
rule pileup_borders:
    input:
        cool = join(OUT, 'cool', '{sample}_hic_' + f'{MAX_RES_STR}.cool'),
        border = join(L_OUT, '{sample}', 'borders_out.txt')
    output:
        pileup = join(L_OUT, 'pileup', '{sample}.txt'),
        bedpe = join(L_TMP, '{sample}_borders.bed2d')
    threads: 6
    run:
        from os.path import dirname, basename, join
        # Format bedpe to have only required colummns
        shell("tail -n +2 {input.border} | cut -f1-6 > {output.bedpe}")
        # Make a pileup for each sample - distance bin combination
        shell(
            "coolpup.py --n_proc {threads} \
                --mindist 40000 \
                --maxdist 2000000 \
                --outname {file} \
                --outdir {dir} \
                {cool} \
                {bedpe}".format(
                    threads=threads,
                    dir=dirname(output['pileup']),
                    file=basename(output['pileup']),
                    cool=input['cool'],
                    bedpe=output['bedpe']
                )
        )

# Generate one figure per distance, with one panel per sample on each figure.
rule pileup_figure:
    input: expand(join(L_OUT, 'pileup', '{sample}.txt'), sample=np.unique(units.loc[units.libtype == "hic", "sample"]))
    output: join(OUT, 'figures', 'borders.pdf')
    run:
        from matplotlib import pyplot as plt
        from os.path import basename
        import numpy as np
        import matplotlib
        matplotlib.use('Agg')
        fig, axi = plt.subplots(len(input)//2, 2)
        for i, ax in enumerate(axi.reshape(-1)):
            a = np.loadtxt(input[i])
            ax.imshow(np.log2(a), vmin=-2, vmax=2)
            fname = basename(input[i])
            ax.set_title(fname.split('_')[0])
        fig.savefig(output[0])


rule concat_border_stats:
    input: expand(join(L_OUT, '{sample}', 'borders_quant.txt'), sample=hic_libs)
    output: join(L_OUT, 'all_borders.tsv')
    run:
        from os.path import dirname, basename
        libs = []
        for lib in input:
            lib_df = pd.read_csv(lib, sep='\t')
            lib_df['sample'] = basename(dirname(lib))
            libs.append(lib_df)

        borders = pd.concat(libs, axis=0, ignore_index=True)
        borders.to_csv(output[0], sep='\t', index=False)
        

rule plot_border_stats:
    input: join(L_OUT, 'all_borders.tsv')
    output: join(OUT, 'figures', 'border_stats.pdf')
    shell: "touch {output}"
