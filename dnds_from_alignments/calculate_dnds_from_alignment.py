# SJS
# This script computes a site-wise dN/dS, based on Spielman and Wilke (2015), MBE for each alignment based on amino-acid frequencies. Mutation rates are assumed equal.
# Dependencies: numpy, biopython, pyvolve, compute_dnds_from_mutsel module (packaged w/ pyvolve)

from compute_dnds_from_mutsel import *
from pyvolve import Genetics
from pyvolve import state_freqs
g = Genetics()
from Bio import AlignIO
import os
import re
from copy import deepcopy
pitype = "unequalpi"
blank_codons = {'ACC': 0, 'ATG': 0, 'AAG': 0, 'AAA': 0, 'ATC': 0, 'AAC': 0, 'ATA': 0, 'AGG': 0, 'CCT': 0, 'CTC': 0, 'AGC': 0, 'ACA': 0,
                'AGA': 0, 'CAT': 0, 'AAT': 0, 'ATT': 0, 'CTG': 0, 'CTA': 0, 'ACT': 0, 'CAC': 0, 'ACG': 0, 'CAA': 0, 'AGT': 0, 'CCA': 0,
                'CCG': 0, 'CCC': 0, 'TAT': 0, 'GGT': 0, 'TGT': 0, 'CGA': 0, 'CAG': 0, 'CGC': 0, 'GAT': 0, 'CGG': 0, 'CTT': 0, 'TGC': 0,
                'GGG': 0, 'GGA': 0, 'TGG': 0, 'GGC': 0, 'TAC': 0, 'GAG': 0, 'TCG': 0, 'TTA': 0, 'GAC': 0, 'CGT': 0, 'GAA': 0, 'TCA': 0,
                'GCA': 0, 'GTA': 0, 'GCC': 0, 'GTC': 0, 'GCG': 0, 'GTG': 0, 'TTC': 0, 'GTT': 0, 'GCT': 0, 'TTG': 0, 'TCC': 0, 'TCT': 0, 'TTT':0}

pi_a = 0.32; pi_t = 0.28; pi_c = 0.18; pi_g = 0.22;
kappa = 4.0
mu_dict = {'AT': pi_t, 'TA':pi_a, 'CG': pi_g, 'GC':pi_c, 'AC': pi_c, 'CA':pi_a, 'GT':pi_t, 'TG':pi_g, 'AG': kappa*pi_g, 'GA':kappa*pi_a, 'CT':kappa*pi_t, 'TC':kappa*pi_c}


def compute_dnds(d, mu):
    c = dNdS_from_MutSel(d, mu)
    dn = None
    ds = None
    # Assertion error will be raised in dnds calculation if no evolution has occurred. These are uninformative sites.
    try:
        dn = c.compute_dn()
        if np.isnan(dn) or np.isinf(dn):
            dn = "NA"
    except AssertionError:
        dnds = "NA"
    try:
        ds = c.compute_ds()
        if np.isnan(ds) or np.isinf(ds):
            ds = "NA"
    except AssertionError:
        ds = "NA"

    assert(dn != None), "dn not computed."
    assert(ds != None), "ds not computed."
    return dn, ds



def extract_dnds(dir, file):
    with open(dir + file, "r") as f:
        aln = AlignIO.read(dir + file, "fasta")
    get_info = re.search("rep(\d+)_n(\d+)_bl(0\.\d+)_unequalpi_(\w+)\.fasta", file)
    if get_info:
        rep = get_info.group(1)
        ntaxa = get_info.group(2)
        bl = get_info.group(3)
        bias = get_info.group(4)
    outlines = []
    for x in range(0, len(aln[0]), 3):
        # Create dictionary of codon frequencies
        column = list(aln[:, x:x+3])
        col_codons = deepcopy(blank_codons)
        for c in column:
            c2 = str(c.seq)
            col_codons[c2] += 1./len(aln)
        assert( sum(col_codons.values()) == 1.), "\nError in frequency calculations from an alignment column."

        # Calculate dn,ds for column
        dn, ds = compute_dnds(col_codons, mu_dict)
        outline = str(x/3 + 1) + "," + ntaxa + "," + bl + "," + rep + "," + bias + "," + str(dn) + "," + str(ds) + ",empirical"
        outlines.append(outline)
    return outlines

def main():

    alignment_directory = "/Users/sjspielman/Dropbox/dnds1rate2rate_data_results/alignments/balanced_trees/" #sys.argv[1]
    outfile = "dnds_from_balanced_alignments.csv" # unequalpi only
    with open(outfile, "w") as f:
        f.write("site,ntaxa,bl,rep,type,dn,ds,method\n") #method is "empirical"
    pitype = "unequalpi"
    files = os.listdir(alignment_directory)
    for file in files:
        if pitype in file and file.endswith(".fasta"):
            print file
            outlines = extract_dnds(alignment_directory, file)
            with open(outfile, "a") as f:
                f.write("\n".join(outlines) + "\n")
main()
