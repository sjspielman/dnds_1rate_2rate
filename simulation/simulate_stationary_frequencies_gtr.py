# SJS 
# Script to simulate 100 sets of stationary codon frequencies with mutational gtrmetry, specifically with yeast mutation rates from Zhu et al. (2014) PNAS
# Uses the same codon fitness values as the nobias set

import numpy as np
import csv
from dnds_mutsel_functions import *
ZERO = 1e-8

# Import regular frequencies and lambda values
raw_freqfile = "codon_freq_lib_nobias.txt"
raw_codon_freqs = np.loadtxt(raw_freqfile)

# Output file names
gtr_freqfile = "codon_freq_lib_gtr.txt"

# gtr_rates = np.random.gamma(2.,0.5, 6)
# pi_a = 0.3
# pi_t = 0.3
# pi_c = 0.2
# pi_g = 0.2
# mu_dict = {'AG':gtr_rates[0]*pi_g, 'GA':gtr_rates[0]*pi_a, 'CT':gtr_rates[1]*pi_t, 'TC':gtr_rates[1]*pi_c, 'AC':gtr_rates[2]*pi_c, 'CA':gtr_rates[2]*pi_a, 'TG':gtr_rates[3]*pi_g, 'GT':gtr_rates[3]*pi_t, 'AT':gtr_rates[4]*pi_t, 'TA':gtr_rates[4]*pi_a, 'GC':gtr_rates[5]*pi_c, 'CG':gtr_rates[5]*pi_g}
mu_dict = {'AC': 0.12806341921812744, 'GT': 0.036330424021680262, 'AG': 0.11300569340489693, 'CA': 0.19209512882719118, 'CG': 0.023165377637234918, 'GC': 0.023165377637234918, 'AT': 0.17689851535472259, 'GA': 0.16950854010734537, 'CT': 0.51939282525298114, 'TG': 0.024220282681120177, 'TC': 0.3462618835019875, 'TA': 0.17689851535472259}



gtr_codon_freqs = []
for raw_freqs in raw_codon_freqs:
    
    # Determine codon fitness values
    codon_fit = np.log(raw_freqs)
    
    # Extract equilibrium frequencies
    matrix = build_mutsel_matrix(mu_dict, codon_fit)
    cf = get_eq_from_eig(matrix) 
    assert( abs(np.sum(cf) - 1.) < ZERO ), "codon frequencies do not sum to 1" 
    gtr_codon_freqs.append(cf)
    
np.savetxt(gtr_freqfile, np.array(gtr_codon_freqs))

    
    
    
    
    
    
    
    
    
    
    




