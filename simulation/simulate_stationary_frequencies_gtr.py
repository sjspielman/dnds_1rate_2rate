# SJS 
# Script to simulate 100 sets of stationary codon frequencies with under a GTR mutation model
# Uses the same codon fitness values as the original nobias dataset

import numpy as np
import csv
from dnds_mutsel_functions import *
ZERO = 1e-8

# Import regular frequencies and lambda values
raw_freqfile = "codon_freq_lib_nobias.txt"
raw_codon_freqs = np.loadtxt(raw_freqfile)

# Output file names
gtr_freqfile = "codon_freq_lib_gtr.txt"

# The GTR rates were derived with the following code:
#### gtr_rates =  np.sort( np.random.gamma(1., 1.5, 6))
#### gtr_rates = gtr_rates[::-1] # This way largest are first, as transitions are first in the dictionary.
# These pi values were simply chosen to create some amount of A/T bias
pi_a = 0.32
pi_t = 0.34
pi_c = 0.16
pi_g = 0.18
gtr_rates = [ 1.64390601,  1.27668478,  0.795571  ,  0.44377381,  0.32759197,  0.25651819]
mu_dict = {'AG':gtr_rates[0]*pi_g, 'GA':gtr_rates[0]*pi_a, 'CT':gtr_rates[1]*pi_t, 'TC':gtr_rates[1]*pi_c, 'AC':gtr_rates[2]*pi_c, 'CA':gtr_rates[2]*pi_a, 'TG':gtr_rates[3]*pi_g, 'GT':gtr_rates[3]*pi_t, 'AT':gtr_rates[4]*pi_t, 'TA':gtr_rates[4]*pi_a, 'GC':gtr_rates[5]*pi_c, 'CG':gtr_rates[5]*pi_g}



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

    
    
    
    
    
    
    
    
    
    
    



