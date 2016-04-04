# SJS 
# Script to simulate 100 sets of stationary codon frequencies without and with codon bias, under an HKY model with unequal nucleotide frequencies.


import sys
from numpy import *
from function_library import *

# For computing new frequencies
pi_a = 0.32
pi_t = 0.28
pi_c = 0.18
pi_g = 0.22
kappa = 4.0
mu_dict = {'AT': pi_t, 'TA':pi_a, 'CG': pi_g, 'GC':pi_c, 'AC': pi_c, 'CA':pi_a, 'GT':pi_t, 'TG':pi_g, 'AG': kappa*pi_g, 'GA':kappa*pi_a, 'CT':kappa*pi_t, 'TC':kappa*pi_c}




for type in ["bias", "nobias"]:
    rawfreqs = np.loadtxt("codon_freq_lib_equalpi_" + type + ".txt")
    freqfile = "codon_freq_lib_unequalpi_" + type + ".txt"
    new_freqs = []
    
    for f in rawfreqs:
    
        fitness = np.log(f)
        
        matrix = build_mutsel_matrix(mu_dict, fitness)
        cf = get_eq_from_eig(matrix) 
        assert( abs(np.sum(cf) - 1.) < ZERO ), "new codon frequencies do not sum to 1" 
        new_freqs.append( cf )
      
    np.savetxt(freqfile, np.array(new_freqs))
