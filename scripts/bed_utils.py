# Functions for operating on bed / bedgraph files
# cmdoret, 20190507
import csv
import numpy as np

def merge_bedgraphs(chrom_sizes, output, resolution=5000,  mode='avg',*args):
    """
    Combined an arbitrary number of bedgraph files into a single one.
    The bedgraph files can have different resolutions and will be forced into
    the selected resolution. Assumes the files have no header.

    NOTE : This is currently achieved by bedtools, this code might be removed eventually
    Parameters
    ----------
    chrom_sizes : list of tuples
        A list of (chrom_name: chromsize) tuples. [(str, int), ...]
    output : str
        Path to the output merged bedgraph file
    resolution : int
        The bin size in the output BEDgraph file.
    mode : str
        How to transform signal when splitting / merging bins. Can be 'avg' to
        average the signal of merged bins and keep original signal when
        splitting bins. Can also be set to 'sum' to sum signals when merging and
        divide signal when splitting.
    args* : str
        Filenames of the BED files to combine
    """
    # Define functions to apply on signal when merging or splitting bins
    try:
        split_fun = {'avg': lambda x: x, 'sum': np.divide}
        merge_fun = {'avg': np.mean, 'sum': np.sum}
        split_fun, merge_fun = split_fun[mode], merge_fun[mode]
    except IndexError as err:
        print("Invalid value for 'mode': {0}".format(mode))
        raise err
    
    # Open all input files and handles into a list
    input_handles = {in_bg: open(in_bg,  'r') for in_bg in *args}

    # Make a reader object with each bedgraph to reference columns by name
    bgcols = ['chrom', 'start', 'end', 'signal']
    for in_bg in input_handles:
        input_handles[in_bg] = csv.DictReader(input_handles[in_bg], delimiter='\t', fieldnames=bgcols)
            
    curr_bin = ["", 0, 0] + [0] * len(input_handles)
    # Loop over chromosomes
    for chrom, size in chrom_sizes:
        # Split each chromosome into equally sized bins
        for bed_bin in range(0, size, resolution):
            for in_bg in input_handles:
                while 
    