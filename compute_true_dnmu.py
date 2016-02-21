import sys
from compute_dnds_from_mutsel import *

results_directory = "../results/"

types = ["nobias", "bias", "asym"]
for type in types:

    print "Calculuting true dN/dS for ", type, "simulations"
    
    freqfile = results_directory + "codon_freq_lib_" + type + ".txt"
    infofile = results_directory + "codon_freq_lib_" + type + "_info.txt"

    true_frequencies = np.loadtxt(freqfile)

    dn = np.zeros(len(true_frequencies))
    mu = np.zeros(len(true_frequencies))


    if type == "asym":
        m = 1.67e-10 # this is the mean per generation per nucleotide mutation rate. Ultimately, as this is relative, it doesn't matter.
        mu_dict = {'AG':0.144/2*m, 'TC':0.144/2*m, 'GA':0.349/2*m, 'CT':0.349/2*m, 'AC':0.11/2*m, 'TG':0.11/2*m, 'CA':0.182/2*m, 'GT':0.182/2*m, 'AT':0.063/2*m, 'TA':0.063/2*m, 'GC':0.152/2*m, 'CG':0.152/2*m}

    else:
        m = 1.
        kappa = 4.0
        mu_dict = {'AT': m, 'TA':m, 'CG': m, 'GC':m, 'AC': m, 'CA':m, 'GT':m, 'TG':m, 'AG': kappa*m, 'GA':kappa*m, 'CT':kappa*m, 'TC':kappa*m}




    for i in range(len(true_frequencies)):
        freqs = true_frequencies[i]
        m = dNdS_from_MutSel(freqs, mu_dict)
        dn[i] = m.compute_dn()
        mu[i] = m.compute_mu()

    mu_global = np.mean(mu)
    dnmu2 = dn / mu
    dnmu1 = dn / mu_global


    with open(infofile, "w") as f:
        f.write("site\tdnmu2\tdnmu1\tdn\tmu\n")
        for i in range(len(dnmu2)):
            f.write(str(i+1) + "\t" + str(dnmu2[i]) + "\t" + str(dnmu1[i]) + "\t" + str(dn[i]) + "\t" + str(mu[i]) + "\n")


