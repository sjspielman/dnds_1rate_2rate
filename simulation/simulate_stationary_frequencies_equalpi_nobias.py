# SJS 
# Script to simulate 100 sets of stationary codon frequencies. Frequencies follow an exponential distribution, ala Ramsey paper.
# dN/dS values are all in range 0.05 - 0.95, evenly spaced

import sys
from numpy import *
from random import shuffle
from function_library import *

def sim_aa_freqs(lambda_):
    aa_freqs = np.zeros(20)
    
    amino_acids  = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
    shuffle(amino_acids)
    
    for i in range(20):
        aa_freqs[i] = np.exp(-1. * lambda_ * (i+1))
    aa_freqs /= np.sum(aa_freqs)
    return aa_freqs 
        

inc = 0.1
min = 0.025
threshold = min + inc
w_per_cat = 11
n = 9

freqs  = []
lambdas = []

# Simulate
freqfile   = "codon_freq_lib_equalpi_nobias.txt"
lambdafile = "lambda_term.csv" # To save the exponents used to compute frequency vectors

total = 0.
for i in range(n):
    print i, threshold
    
    # We need 100, so just give this category another dN/dS
    if threshold == 1.:
        w_per_cat /= 2
    j = 0
    while j < w_per_cat:
        lambda_ = np.random.uniform(low=ZERO, high=3.)
        aa_freqs = sim_aa_freqs(lambda_)
        cf, cf_dict = aa_to_codon_freqs(aa_freqs)
        w = derive_dnds(cf_dict)
        if (w >= threshold - inc and w < threshold and w > min):
            j += 1
            freqs.append(cf)
            lambdas.append(lambda_) 
    threshold += inc
 
# Finally, add an amino-acid neutral
lambda_ = 0.
aa_freqs = sim_aa_freqs(lambda_)
cf, cf_dict = aa_to_codon_freqs(aa_freqs)
freqs.append(cf)
lambdas.append(lambda_) 

# Save all
np.savetxt(freqfile, np.array(freqs))
with open(lambdafile, "w") as outf:
    outf.write("site,lambda\n")
    for i in range(len(dnds)):
        outf.write(str(i + 1) + "," + str(lambdas[i]) + "\n")

