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
    ds = np.zeros(len(true_frequencies))


    if type == "asym":
        mu = 1.67e-10 # this is the mean per generation per nucleotide mutation rate. Ultimately, as this is relative, it doesn't matter.
        mu_dict = {'AG':0.144/2*mu, 'TC':0.144/2*mu, 'GA':0.349/2*mu, 'CT':0.349/2*mu, 'AC':0.11/2*mu, 'TG':0.11/2*mu, 'CA':0.182/2*mu, 'GT':0.182/2*mu, 'AT':0.063/2*mu, 'TA':0.063/2*mu, 'GC':0.152/2*mu, 'CG':0.152/2*mu}

    else:
        mu = 1.
        kappa = 4.0
        mu_dict = {'AT': mu, 'TA':mu, 'CG': mu, 'GC':mu, 'AC': mu, 'CA':mu, 'GT':mu, 'TG':mu, 'AG': kappa*mu, 'GA':kappa*mu, 'CT':kappa*mu, 'TC':kappa*mu}




    for i in range(len(true_frequencies)):
        freqs = true_frequencies[i]
        m = dNdS_from_MutSel(freqs, mu_dict)
        dn[i] = m.compute_dn()
        ds[i] = m.compute_ds()

    ds_global = np.mean(ds)
    dnds_2rate = dn / ds
    dnds_1rate = dn / ds_global


    with open(infofile, "w") as f:
        f.write("site\tdnds2\tdnds1\tdn\tds_site\tds_global\n")
        for i in range(len(dnds_2rate)):
            f.write(str(i+1) + "\t" + str(dnds_2rate[i]) + "\t" + str(dnds_1rate[i]) + "\t" + str(dn[i]) + "\t" + str(ds[i]) + "\t" + str(ds_global) + "\n")


