from count_simulated_dnds import *
import re
import os
import sys

# For computing new frequencies
pi_a = 0.32
pi_t = 0.28
pi_c = 0.18
pi_g = 0.22
kappa = 4.0
mu_dict = {'AT': pi_t, 'TA':pi_a, 'CG': pi_g, 'GC':pi_c, 'AC': pi_c, 'CA':pi_a, 'GT':pi_t, 'TG':pi_g, 'AG': kappa*pi_g, 'GA':kappa*pi_a, 'CT':kappa*pi_t, 'TC':kappa*pi_c}

alnfile = sys.argv[1]
treefile = sys.argv[2]
outfile = sys.argv[3]

c = dNdS_Counter(ALNDIR + alnfile, treefile, mu_dict)
c.calculate_dnds(savefile = outfile)
