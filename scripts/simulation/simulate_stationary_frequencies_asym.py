# SJS 
# Script to simulate 100 sets of stationary codon frequencies with mutational asymmetry, specifically with yeast mutation rates from Zhu et al. (2014) PNAS
# Uses the same codon fitness values as the nobias set

import numpy as np
import csv
from dnds_mutsel_functions import *
ZERO = 1e-8

# Import regular frequencies and lambda values
raw_freqfile = "codon_freq_lib_nobias.txt"
raw_codon_freqs = np.loadtxt(raw_freqfile)

# Output file names
asym_freqfile = "codon_freq_lib_asym.txt"

# Yeast asymmetric mutation rates, from Zhu et al. 2014 (PNAS)
mu = 1.67e-10 # this is the mean per generation per nucleotide mutation rate. Ultimately, as this is relative, it doesn't matter.
mu_dict = {'AG':0.144/2*mu, 'TC':0.144/2*mu, 'GA':0.349/2*mu, 'CT':0.349/2*mu, 'AC':0.11/2*mu, 'TG':0.11/2*mu, 'CA':0.182/2*mu, 'GT':0.182/2*mu, 'AT':0.063/2*mu, 'TA':0.063/2*mu, 'GC':0.152/2*mu, 'CG':0.152/2*mu}


asym_codon_freqs = []
for raw_freqs in raw_codon_freqs:
    
    # Determine codon fitness values
    codon_fit = np.log(raw_freqs)
    
    # Extract equilibrium frequencies
    matrix = build_mutsel_matrix(mu_dict, codon_fit)
    cf = get_eq_from_eig(matrix) 
    assert( abs(np.sum(cf) - 1.) < ZERO ), "codon frequencies do not sum to 1" 
    asym_codon_freqs.append(cf)
    
np.savetxt(asym_freqfile, np.array(asym_codon_freqs))

    
    
    
    
    
    
    
    
    
    
    




