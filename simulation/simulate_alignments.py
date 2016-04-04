# SJS
# Simulate alignments along to specified tree according to MutSel model. 

import numpy as np
import os
import sys
import re
from Bio import AlignIO
from pyvolve import *
codons = ["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]

assert(len(sys.argv) == 5), "Usage: python simulate_alignments.py <treefile> <type> <seqoutfile> <seqanoutfile>"
treefile        = sys.argv[1]
type            = sys.argv[2]
seq_outfile     = sys.argv[3]
seq_anc_outfile = sys.argv[4]

allowed_types = ["equalpi_nobias", "equalpi_bias", "unequalpi_nobias", "unequalpi_bias"]
assert(type in allowed_types), "\nBad simulation type specified."

if type.startswith("equalpi"):
    pi_a = 0.25; pi_t = 0.25; pi_c = 0.25; pi_g = 0.25;
elif type.startswith("unequalpi"):
    pi_a = 0.32; pi_t = 0.28; pi_c = 0.18; pi_g = 0.22;
kappa = 4.0
mu_dict = {'AT': pi_t, 'TA':pi_a, 'CG': pi_g, 'GC':pi_c, 'AC': pi_c, 'CA':pi_a, 'GT':pi_t, 'TG':pi_g, 'AG': kappa*pi_g, 'GA':kappa*pi_a, 'CT':kappa*pi_t, 'TC':kappa*pi_c}

freqfile = "codon_freq_lib_" + type + ".txt"
codon_freqs = np.loadtxt(freqfile)
all_partitions = []

# Set up partitions
for i in range(len(codon_freqs)):
    m = Model("mutsel", {'state_freqs':codon_freqs[i], 'mu':mu_dict})
    p = Partition(models = [m], size = 1)
    all_partitions.append( p )
     
     
# Evolve sequences and save
tree = read_tree(file = treefile)
evolve = Evolver(partitions = all_partitions, tree = tree)
evolve(seqfile = seq_anc_outfile, infofile = False, ratefile = False, write_anc = True)
tipseqs = evolve.get_sequences()
with open(seq_outfile, "w") as outf:
    for record in tipseqs:
        outf.write(">" + record + "\n" + tipseqs[record] + "\n")
        
