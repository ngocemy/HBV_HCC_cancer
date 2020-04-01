# Comparative analysis of A/B compartments

# Make a combined table with hbv cancer and healthy lines
lib = dual_samples + comp_lines
tissue = dual_cells + comp_lines
pop = ["hbv_cancer"] * len(dual_samples) + ["normal"] * len(comp_lines)
merged_samples = pd.DataFrame({"library": lib, "population": pop, 'tissue': tissue})

rule call_compartments:
    input: join(OUT, 'cool', '{sample}_hic_' + f'{COMP_RES_STR}.cool')
    output: join(OUT, 'compartments', 'compartments_{sample}.bedgraph')
    params:
        comp_res = COMP_RES
    run:
        compartment_to_bedgraph(
            input[0], 
            output[0]
        )

# 03b: Use cooltools to compute eigenvectors on cool files and then rank and orient
# them based on gene count in hg19 (fetched online).
rule call_compartments_comp_lines:
  input: join(IN, 'comp_lines', '{tissue}_' + f'{COMP_RES_STR}.cool')
  output: join(OUT, 'compartments', 'compartments_{tissue}.bedgraph')
  run:
    compartment_to_bedgraph(input[0], output[0])



# Combine all compartment files for comparative analysis
rule merge_compartments:
  input : expand(join(OUT, 'compartments', 'compartments_{sample}.bedgraph'), sample=dual_samples + comp_lines)
  output: join(OUT, 'compartments', 'merged_compartments.tsv')
  params:  
    libs = dual_libs + comp_lines
  run:
    # Run bin coordinates from ref file
    bedgraph = pd.read_csv(
      input[0], 
      sep='\t', 
      usecols=[0, 1, 2], 
      names=['chrom', 'start', 'end'], 
      header=None
    )
    # Append eigenvector of libraries iteratively
    for i, file in enumerate(input[:]):
      bedgraph[params['libs'][i]] = pd.read_csv(file, sep='\t', usecols=[3])
    bedgraph.to_csv(output[0], index=False, sep='\t')


rule plot_eigenvectors:
    input: join(OUT, 'compartments', 'merged_compartments.tsv')
    output: join(OUT, 'plots', 'compartments', 'eigens', '{chrom}_eigen.pdf')
    params:
      samples_df = merged_samples,
      chrom = lambda w: f"{w.chrom}"
    run:
      import matplotlib
      matplotlib.use('Agg') 
      comp_df = pd.read_csv(input[0], sep='\t')
      plot_eigens(comp_df, params['samples_df'], chrom=params['chrom'], out=output[0])

# Measure compartment similarity on whole genome using PCA
rule plot_compartment_pca:
  input: join(OUT, 'compartments', 'merged_compartments.tsv')
  output: join(OUT, 'plots', 'compartments', 'compartments_pca.pdf')
  params:
    samples_df = merged_samples
  run:
    import matplotlib
    matplotlib.use('Agg')
    comp_df = pd.read_csv(input[0], sep='\t')
    pca_compartments(comp_df, params['samples_df'], out=output[0])
  

rule compute_compartments_divergence:
  input: join(OUT, 'compartments', 'merged_compartments.tsv')
  output: join(OUT, 'compartments', 'compartments_divergence.tsv')
  params:
    merged_samples = merged_samples,
    winsize = 20,
    max_boots=np.inf
  threads : 12
  run:
    comp_df = pd.read_csv(input[0], sep='\t')
    div_df = compute_compartments_div(comp_df, params['merged_samples'], win_size=params['winsize'], max_boots=params['max_boots'], cpus=threads)
    div_df.to_csv(output[0], sep='\t', index=False)

rule plot_compartments_divergence:
  input: join(OUT, 'compartments', 'compartments_divergence.tsv')
  output: join(OUT, 'plots', 'compartments', 'compartments_divergence.pdf')
  params:
    winsize = 20
  run:
    import matplotlib
    matplotlib.use('Agg')
    div_df = pd.read_csv(input[0], sep='\t')
    plot_compartments_div(div_df, win_size=params['winsize'], out=output[0])
