# SJS 
# Script to simulate 100 sets of stationary codon frequencies with codon bias. 
# Uses the same amino-acid frequencies as the nobias set, but values divvied up differently among codons rather than evenly
# Codon bias is implemented by randomly selecting a preferred codon for each amino acid and incrementing it by a certain value, while the remaining unpreferred are decreased by a certain value.


import sys
from numpy import *
from random import shuffle
from dnds_mutsel_functions import *
        

# For computing dN/dS
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



freqs   = []
bias    = []
dnds    = []
unbiased_freqs = np.loadtxt("codon_freq_lib_nobias.txt")
freqfile = "codon_freq_lib_bias_gtr.txt"
biasfile = "codonbias_gtr_term.txt" # To save the bias "factors" 

rep = 0
for rawfreqs in unbiased_freqs:
    
    w = 1.5
    while w > 1.:
        b = np.random.uniform(0.6, 0.9)
        cf, cf_dict = add_bias(rawfreqs, b)
        w = derive_dnds(cf_dict, mu_dict)
    dnds.append(w)
    freqs.append(cf)
    bias.append(b)

print dnds
 
# Save all
np.savetxt(freqfile, np.array(freqs))
with open(biasfile, "w") as outf:
    outf.write("site\tbias_term\n")
    for i in range(len(bias)):
        outf.write(str(i + 1) + "\t" + str(bias[i]) + "\n")

