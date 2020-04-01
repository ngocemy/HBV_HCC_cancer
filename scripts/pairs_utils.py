#!/bin/env python3
# Utility script to compute common things on .pairs files
# 20190430, cmdoret
from hicstuff import io as hio
from hicstuff import hicstuff as hcs
from collections import OrderedDict
import pandas as pd
import numpy as np
import csv


def chrom_contact_to_bedgraph(
    pairs_file, bg_file, focus_chrom, binsize=5000, threads=1, buffer="2G"
):
    """
    Compute 1D coverage track of contacts with a given chromosome from a pairs
    file and writes it to a bed file.

    Parameters
    ----------
    pairs_file : str
        Path to the input pairs file. It is tab separated with header lines
        starting by "#" and contains 7 columns: readID,chr1,pos1,chr2,pos2,
        strand1,strand2. Assumes the pairs is an upper triangle (chr1,pos1 is
        always smaller than chr2,pos2).
    bg_file : str
        Path to the output bedgraph file. It is tab separated and contains 4
        columns: chrom,start,end,value.
    focus_chrom : str
        Name of the chromosome with which contacts should be computed.
    bin_size : int
        Size of bins in which to compute contacts, in basepairs.
    
    Examples
    --------
    chrom_contact_to_bedgraph('libxyz.pairs', 'libxyz.bedgraph' 'HBV')
    """
    # Sort pairs file by coordinate and check how many lines are in header row
    hio.sort_pairs(
        pairs_file,
        pairs_file + ".sorted",
        ["chr1", "pos1", "chr2", "pos2"],
        threads=threads,
        buffer=buffer,
    )
    header = hio.get_pairs_header(pairs_file + ".sorted")
    # Create bedgraph bins for each chromosome
    bg_bins = OrderedDict()
    for line in header:
        if line.startswith("#chromsize"):
            chrom = line.split()[1:]
            # Initialize empty bins for each chromosome (at least 1)
            bg_bins[chrom[0]] = np.zeros(max(int(chrom[1]) // binsize, 1))

    # bg_bins should look like that: {'chr1': [0, 0, 0, ...], 'chr2':...}
    # Where each array has as many 0's as there are bins in the chrom.

    with open(pairs_file + ".sorted", "r") as pairs:
        # Ignore header lines
        for _ in range(len(header)):
            next(pairs)
        # Define reader object for input pairs file.
        pairs_cols = ["readID", "chr1", "pos1", "chr2", "pos2", "strand1", "strand2"]
        pairs_reader = csv.DictReader(pairs, delimiter="\t", fieldnames=pairs_cols)

        # Read and record pair contacts with focus chromosome
        for pair in pairs_reader:
            try:
                # If chromosome of interest is first mate, increment second mate
                # bin coverage by 1
                if pair["chr1"] == focus_chrom:
                    bg_bins[pair["chr2"]][int(pair["pos2"]) // binsize] += 1
                # If it is second mate, increment first mate bin coverage
                elif pair["chr2"] == focus_chrom:
                    bg_bins[pair["chr1"]][int(pair["pos1"]) // binsize] += 1
            except IndexError:
                print(pair)

    # Write contact value as a 1D bedgraph file
    with open(bg_file, "w") as bg:
        for chrom in bg_bins:
            for bin_idx, bin_contacts in enumerate(bg_bins[chrom]):
                # Write single line as "chrom start end contacts"
                bg.write(
                    "{0}\t{1}\t{2}\t{3}\n".format(
                        chrom, bin_idx * binsize, (bin_idx + 1) * binsize, bin_contacts
                    )
                )
