# SJS
# Simulate according to an actual dN/dS model (MG94)

import numpy as np
import os
import sys
import re
from Bio import AlignIO
from pyvolve import *
codons = ["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]

treefile        = sys.argv[1]
pitype = "unequalpi"
length = 100
if pitype == "equalpi":
    pi_a = 0.25; pi_t = 0.25; pi_c = 0.25; pi_g = 0.25;
elif pitype == "unequalpi":
    pi_a = 0.32; pi_t = 0.28; pi_c = 0.18; pi_g = 0.22;
kappa = 4.0
mu_dict = {'AT': pi_t, 'TA':pi_a, 'CG': pi_g, 'GC':pi_c, 'AC': pi_c, 'CA':pi_a, 'GT':pi_t, 'TG':pi_g, 'AG': kappa*pi_g, 'GA':kappa*pi_a, 'CT':kappa*pi_t, 'TC':kappa*pi_c}

dn = np.random.uniform(low=0.01, high = 1.5, size = length)
ds = np.random.uniform(low=0.01, high = 1.5, size = length)
m = Model("MG", {'mu': mu_dict, 'beta': dn, 'alpha': ds})
p = Partition(models = [m], size = length)


# Evolve sequences and save
tree = read_tree(file = treefile)
evolve = Evolver(partitions = p, tree = tree)
evolve(seqfile = "dnds_sim1.fasta", infofile = "dnds_sim1_info.txt", ratefile = "dnds_sim1_rates.txt")
#tipseqs = evolve.get_sequences()
#with open(seq_outfile, "w") as outf:
#    for record in tipseqs:
#        outf.write(">" + record + "\n" + tipseqs[record] + "\n")
