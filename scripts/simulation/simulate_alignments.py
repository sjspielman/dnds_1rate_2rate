# SJS
# Simulate alignments along to specified tree according to MutSel model. 

import numpy as np
import os
import sys
import re
from Bio import AlignIO
from pyvolve import *
codons = ["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]

assert(len(sys.argv) == 3), "Usage: python simulate_alignments.py <treefile> <type>"
treefile        = sys.argv[1]
type            = sys.argv[2]
seq_anc_outfile = sys.argv[3]
seq_outfile     = sys.argv[4]


# Setup mutation rates and state frequencies
if type in ["bias", "nobias"]:
    kappa = 4.0; mu = 1.
    mu_dict = {'AT': mu, 'TA':mu, 'CG': mu, 'GC':mu, 'AC': mu, 'CA':mu, 'GT':mu, 'TG':mu, 'AG': kappa*mu, 'GA':kappa*mu, 'CT':kappa*mu, 'TC':kappa*mu}

elif type == "asym":
    mu = 1.67e-10 # this is the mean per generation per nucleotide mutation rate. Ultimately, as this is relative, it doesn't matter.
    mu_dict = {'AG':0.144/2*mu, 'TC':0.144/2*mu, 'GA':0.349/2*mu, 'CT':0.349/2*mu, 'AC':0.11/2*mu, 'TG':0.11/2*mu, 'CA':0.182/2*mu, 'GT':0.182/2*mu, 'AT':0.063/2*mu, 'TA':0.063/2*mu, 'GC':0.152/2*mu, 'CG':0.152/2*mu}

elif type == "gtr":
    mu_dict = {'AC': 0.12806341921812744, 'GT': 0.036330424021680262, 'AG': 0.11300569340489693, 'CA': 0.19209512882719118, 'CG': 0.023165377637234918, 'GC': 0.023165377637234918, 'AT': 0.17689851535472259, 'GA': 0.16950854010734537, 'CT': 0.51939282525298114, 'TG': 0.024220282681120177, 'TC': 0.3462618835019875, 'TA': 0.17689851535472259}

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
        
