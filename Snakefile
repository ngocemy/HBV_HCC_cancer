#!/bin/env snakemake -s
# This file can be run using snakemake. It runs all different HBV cancer 
# analyses on the input pairs files
# cmdoret, 20190501

#from snakemake.utils import validate
from os.path import join
import numpy as np
import pandas as pd
import sys
from glob import glob
import seaborn

### USER DEFINED VARIABLES
# =============================================================================

# Define samples on which the code should run
configfile: "config.yaml"
#validate(config, schema="schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep='\t', dtype=str, comment='#').set_index(["sample"], drop=False)
#validate(samples, schema="schemas/samples.schema.yaml")

units = pd.read_csv(config["units"], sep='\t', dtype=str, comment='#')
# Enforce str in index
#units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])


def bp_to_suffix(size):
    """
    Given a number of basepairs, returns the notation with suffix
    
    Examples
    --------
    >>> get_bp_scale(10000)
    "10kb"
    """
    # Defind mapping between powers of 10 and suffixes
    pow_to_suffix = {0: "bp", 3: "kb", 6: "Mb", 9: "Gb", 12: "Tb"}
    sorted_pows = sorted(pow_to_suffix.keys())
    # Find order of magnitude of input size
    input_len_pow = int(np.log10(size))
    # Find which power matches order of magnitude
    valid_pow_idx = max(0, np.searchsorted(sorted_pows, input_len_pow, side='right') - 1)
    input_valid_pow = sorted_pows[valid_pow_idx]
    # Get corresponding suffix
    suffix = pow_to_suffix[input_valid_pow]
    scale = 10 ** input_valid_pow
    str_bp = f"{int(size // scale)}{suffix}"
    return str_bp

## Function to build remote path to input / output files
def access_remote(local_path):
    """
    Given the local path to a file, return the full path
    for remote access by accessing snakemake config. When using
    SFTP (ssh access to remote files), the bucket parameter is used
    as a prefix to the path.

    Parameters
    ----------
    local_path : str
        The local path within filesystem that need to be wrapped.

    Returns
    -------
        The path with remote information (provider, host, bucket name) prepended.

    """
    rc = config['remote']
    provider = rc['provider']
    bucket = rc['bucket']
    host = rc['host']
    username = rc['username']
    ssh_key = rc['ssh_key']
    password = rc['password']
    # Failsafe, expansion must be done after calling this function
    if not isinstance(local_path, str):
        print("Remote information cannot be added to expanded path. Expand after access_remote.")
        sys.exit(1)
    if provider == 'GS':
        GS = GSRemoteProvider()
        remote_path = GS.remote(join(bucket, local_path))
    elif provider == 'SFTP':
        if password == "":
            SFTP = SFTPRemoteProvider(username=username, private_key=ssh_key)
        else:
            SFTP = SFTPRemoteProvider(username=username, password=password)
        remote_path = SFTP.remote(host + ":22" + join(bucket, local_path))
    elif provider in ('', "local"):
        remote_path = local_path
    else: 
        print('Pipeline not configured for remote provider: {}'.format(provider))
        print('Please edit the config.yaml file to select a valid provider.')
        sys.exit(1)
    return remote_path

# Set input / output paths
DATA_DIR = 'human_hbv_cancer'
IN = join(DATA_DIR, 'input')
OUT = join(DATA_DIR, 'output')
TMP = join(DATA_DIR, 'tmp')
GENOME_HUMAN_INDEX = config['reference']['human']
GENOME_HBV_INDEX = config['reference']['hbv']
MAX_RES = config['contact_maps']['max_res']
MAX_RES_STR = bp_to_suffix(MAX_RES)
COMP_RES = config['contact_maps']['comp_res']
COMP_RES_STR = bp_to_suffix(COMP_RES)

# Libraries that are available in Hi-C
hic_libs = np.unique(units.loc[units.libtype == "hic", "sample"])
captn_libs= np.unique(units.loc[units.libtype == "captn", "sample"])
# Libraries that are available both in RNAseq and Hi-C and their cell types
dual_libs = list(set(hic_libs).intersection(
            set(units.loc[units.libtype == "rnaseq", "sample"])))
# Get sample names and cell types in the cell order
dual_samples = samples['sample'][np.isin(samples['sample'], dual_libs)].tolist() 
dual_cells = samples.loc[np.isin(samples['sample'], dual_libs), 'cell_type'].reset_index(drop=True).tolist() 

# # Load lines from other studies to be used for comparison
# comp_lines = open(join(IN, 'comp_lines', 'tissues.tsv')).read().split('\n')
# comp_lines = [l for l in comp_lines if l != ""]

# Set constraints on wildcards based on data
wildcard_constraints:
  cell_type="|".join(np.unique(samples.cell_type)),
  sample="|".join(samples['sample']),
  libtype="|".join(np.unique(units.libtype)),
  # tissue="|".join(comp_lines)




### ANALYSIS
# =============================================================================

# Final files to generate
rule all:
  input:
    # expand(join(OUT, 'all_signals_{sample}.bedgraph'), sample=dual_libs),
    # join(OUT, 'rnaseq', 'diff_expr', 'integration_vs_control.tsv'),
    # join(OUT, 'figures', 'loops.pdf'),
    # join(OUT, 'figures', 'loop_stats.pdf'),
    # join(OUT, 'compartments', 'merged_compartments.tsv'),
    # join(OUT, 'plots', 'compartments', 'compartments_divergence.pdf'),
    # join(OUT, 'plots', 'compartments', 'compartments_pca.pdf'),
    # expand(join(OUT, 'plots', 'compartments', 'eigens', '{chrom}_eigen.pdf'), chrom=[f"chr{i}" for i in range(1, 23)] + ['chrX'])
    # #expand(join(OUT, 'hint', '{sample}'), sample=dual_libs)
#    expand(join(OUT,'polyidus','{sample}_hic'),sample=hic_libs) 
	expand(join(TMP, 'hetero_reads_{sample}_{libtype}.txt'),sample=dual_libs,libtype=['captn','hic'])
# Python helper functions
include: "scripts/pairs_utils.py"
include: "scripts/mat_utils.py"
include: "scripts/compartments_utils.py"

# Pipeline sub-workflows
include: 'rules/01_common.smk'
# include: 'rules/02_hic_processing.smk'
# include: 'rules/03_compartment_analysis.smk'
# include: 'rules/04_loop_calling.smk'
# include: 'rules/05_rna_seq.smk'
include: 'rules/06_insertion_analysis.smk'

# 04:Combine signals at different resolutions into a single bedgraph
rule aggregate_signals:
  input: 
    comp = join(OUT, 'compartments', 'compartments_{sample}.bedgraph'),
    insul = join(OUT, 'insulation_{sample}.bedgraph'),
    virus = join(OUT, 'cov_virus_{sample}.bedgraph'),
    chroms = join(IN, 'chrom.sizes'),
    rna = join(OUT, 'rnaseq', 'rna_coverage_{sample}.bedgraph')
  output: join(OUT, 'all_signals_{sample}.bedgraph')
  shell:
    """
    mkdir -p {TMP}/signal_tracks

    # Generate output bins as a bed file
    bedtools makewindows -w {MAX_RES} \
                         -g {input.chroms} \
                         > {OUT}/bins.bed
    
    # Intersect each signal with fixed output bins
    # Group overlapping bins and average signal for each group
    for in_file in {input.comp} {input.virus} {input.insul}; do
      tmp_file={TMP}/signal_tracks/$(basename $in_file)
      bedtools intersect -a $in_file -b {OUT}/bins.bed -wb |
        bedtools groupby -g 1,2,3 -c 4 -o mean 2> /dev/null |
        awk '{{print $4}}' > $tmp_file
    done
    
    # RNAseq file have smaller bins than MAX_RES: need to groupby different col
    bedtools intersect -a {input.rna} -b {OUT}/bins.bed -wb |
      sort -k5,5 -k6,6n |
      bedtools groupby -g 5,6,7 -c 4 -o mean 2> /dev/null |
      awk '{{print $4}}' > {TMP}/signal_tracks/$(basename {input.rna})
      
    echo -e "chrom\tstart\tend\tcompartment_{COMP_RES}\tvirus_{MAX_RES}\tlog2_insulation{COMP_RES}\trna_coverage" > {output}
    paste {OUT}/bins.bed {TMP}/signal_tracks/*{wildcards.sample}* >> {output}

    """
