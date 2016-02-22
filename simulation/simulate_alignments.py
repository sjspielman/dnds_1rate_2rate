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


# Setup mutation rates and state frequencies
if type in ["bias", "nobias"]:
    kappa = 4.0; mu = 1.
    mu_dict = {'AT': mu, 'TA':mu, 'CG': mu, 'GC':mu, 'AC': mu, 'CA':mu, 'GT':mu, 'TG':mu, 'AG': kappa*mu, 'GA':kappa*mu, 'CT':kappa*mu, 'TC':kappa*mu}

elif type == "gtr":
    pi_a = 0.32; pi_t = 0.34; pi_c = 0.16; pi_g = 0.18
    gtr_rates = [ 1.64390601,  1.27668478,  0.795571,  0.44377381,  0.32759197,  0.25651819]
    mu_dict = {'AG':gtr_rates[0]*pi_g, 'GA':gtr_rates[0]*pi_a, 'CT':gtr_rates[1]*pi_t, 'TC':gtr_rates[1]*pi_c, 'AC':gtr_rates[2]*pi_c, 'CA':gtr_rates[2]*pi_a, 'TG':gtr_rates[3]*pi_g, 'GT':gtr_rates[3]*pi_t, 'AT':gtr_rates[4]*pi_t, 'TA':gtr_rates[4]*pi_a, 'GC':gtr_rates[5]*pi_c, 'CG':gtr_rates[5]*pi_g}

else:
    raise AssertionError("\n\nWrong simulation type specified.")

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
        
