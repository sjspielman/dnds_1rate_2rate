# SJS 
# Script to simulate 100 sets of stationary codon frequencies with codon bias. 
# Uses the same amino-acid frequencies as the nobias set, but values divvied up differently among codons rather than evenly
# Codon bias is implemented by randomly selecting a preferred codon for each amino acid and incrementing it by a certain value, while the remaining unpreferred are decreased by a certain value.


import sys
from numpy import *
from random import shuffle
from dnds_mutsel_functions import *
        

# For computing dN/dS
mu = 1.
kappa = 4.0
mu_dict = {'AT': mu, 'TA':mu, 'CG': mu, 'GC':mu, 'AC': mu, 'CA':mu, 'GT':mu, 'TG':mu, 'AG': kappa*mu, 'GA':kappa*mu, 'CT':kappa*mu, 'TC':kappa*mu}

freqs   = []
bias    = []

unbiased_freqs = np.loadtxt("../../results/codon_freq_lib.txt")
freqfile = "../../results/codon_freq_lib_bias.txt"
biasfile = "../../results/codonbias_term.txt" # To save the bias "factors" 

rep = 0
for rawfreqs in unbiased_freqs:
    
    w = 1.5
    while w > 1.:
        b = np.random.uniform(0.6, 0.9)
        cf, cf_dict = add_bias(rawfreqs, b)
    
    freqs.append(cf)
    bias.append(b)

    
 
# Save all
np.savetxt(freqfile, np.array(freqs))
with open(biasfile, "w") as outf:
    outf.write("site\tbias_term\n")
    for i in range(len(bias)):
        outf.write(str(i + 1) + "\t" + str(bias[i]) + "\n")

