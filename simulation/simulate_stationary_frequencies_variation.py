# SJS 

# IMPORTANT THOUGHT:
## From where might dS variation arise? If it comes from mutation rate differences across sites, then fitting a global mutation model is already the problem. The assumption is that dS would capture the variation of mutation rate differences, but if mutation rate parameters are already set as identical across sites, then it is unclear exactly what information the dS parameter is capturing.



import numpy as np
import pickle
from dnds_mutsel_functions import *
sys.path.append("../")
from compute_dnds_from_mutsel import *
ZERO = 1e-8

# Import regular frequencies and lambda values
raw_freqfile = "codon_freq_lib_nobias.txt"
raw_codon_freqs = np.loadtxt(raw_freqfile)

# Output file names
freqfile = "codon_freq_lib_variation.txt"
mufile   = "variation_mutation_dictionaries.pkl" # pickle the dictionaries. This is a list of mu_dicts, where list index corresponds to site.
dndsfile = "codon_freq_lib_variation_info.txt" # dnds

def sample_gtr_mutation(nsites, outfile):
    mutation_rates = []
    for x in range(nsites):
        gtr_rates = np.random.uniform(0.1, 10000., 6)
        pi = np.random.uniform(0.1, 5., size=4)
        pi/= np.sum(pi)
        pi_a = 0.01#pi[0]
        pi_t = 0.01#pi[1]
        pi_g = 0.49#pi[2]
        pi_c = 0.49#pi[3]
        print pi
        assert(1. - (pi_a + pi_t + pi_c + pi_g) <= ZERO), "Bad nucleotide frequencies."
        mu_dict = {'AG':gtr_rates[0]*pi_g, 'GA':gtr_rates[0]*pi_a, 'CT':gtr_rates[1]*pi_t, 'TC':gtr_rates[1]*pi_c, 'AC':gtr_rates[2]*pi_c, 'CA':gtr_rates[2]*pi_a, 'TG':gtr_rates[3]*pi_g, 'GT':gtr_rates[3]*pi_t, 'AT':gtr_rates[4]*pi_t, 'TA':gtr_rates[4]*pi_a, 'GC':gtr_rates[5]*pi_c, 'CG':gtr_rates[5]*pi_g}
        print mu_dict
        mutation_rates.append(mu_dict)
    with open(outfile, "w") as f:
        pickle.dump( mutation_rates, f )
    return mutation_rates
    

mutation_rates = sample_gtr_mutation(len(raw_codon_freqs), mufile)
new_codon_freqs = []
dn = []
ds = []

i = 0
for raw_freqs in raw_codon_freqs:
    
    # Determine codon fitness values
    codon_fit = np.log(raw_freqs)
    
    # Extract equilibrium frequencies
    matrix = build_mutsel_matrix(mutation_rates[i], codon_fit)
    cf = get_eq_from_eig(matrix) 
    assert( abs(np.sum(cf) - 1.) < ZERO ), "codon frequencies do not sum to 1" 
    new_codon_freqs.append(cf)
    
    c = dNdS_from_MutSel(cf, mutation_dictionary = mutation_rates[i])
    dn.append( c.compute_dn() )
    ds.append( c.compute_ds() )
    
    i+=1
    
np.savetxt(freqfile, np.array(new_codon_freqs))

outf = open(dndsfile, "w")
outf.write("site\tdn\tds\n")
for i in range(len(dn)):
    outf.write(str(i+1)+ "\t" + str(dn[i]) + "\t" + str(ds[i]) + "\n")
