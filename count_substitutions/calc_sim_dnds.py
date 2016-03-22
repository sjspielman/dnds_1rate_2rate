from count_simulated_dnds import *
import re
import os

pi_a = 0.32; pi_t = 0.34; pi_c = 0.16; pi_g = 0.18
gtr_rates = [ 1.64390601,  1.27668478,  0.795571,  0.44377381,  0.32759197,  0.25651819]
mu_dict = {'AG':gtr_rates[0]*pi_g, 'GA':gtr_rates[0]*pi_a, 'CT':gtr_rates[1]*pi_t, 'TC':gtr_rates[1]*pi_c, 'AC':gtr_rates[2]*pi_c, 'CA':gtr_rates[2]*pi_a, 'TG':gtr_rates[3]*pi_g, 'GT':gtr_rates[3]*pi_t, 'AT':gtr_rates[4]*pi_t, 'TA':gtr_rates[4]*pi_a, 'GC':gtr_rates[5]*pi_c, 'CG':gtr_rates[5]*pi_g}

ALNDIR  = "../data/alignments_with_ancestors/"
TREEDIR = "../data/trees/"
OUTDIR  = "counted/"
files = [x for x in os.listdir(ALNDIR) if x.endswith("withanc.fasta")]

for alnfile in files: 
    print alnfile  
    find = re.search("(rep\d+_)(n\d+_bl0\.\d+)_gtr_withanc.fasta", alnfile)
    treefile = TREEDIR + find.group(2) + ".tre"
    outfile = "counted/" + find.group(1) + "_" + find.group(2) + "_gtr_counted.txt"
    if not os.path.exists(outfile):
        c = dNdS_Counter(ALNDIR + alnfile, treefile, mu_dict)
        c.calculate_dnds(savefile = outfile)  