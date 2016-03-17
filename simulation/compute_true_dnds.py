import sys
from compute_dnds_from_mutsel import *


types = ["gtr", "bias_gtr"]
for type in types:

    print "Calculuting true dN/dS for ", type, "simulations"
    
    freqfile = "codon_freq_lib_" + type + ".txt"
    dndsfile = "codon_freq_lib_" + type + "_true_dnds.txt"

    true_frequencies = np.loadtxt(freqfile)

    dn = np.zeros(len(true_frequencies))
    ds = np.zeros(len(true_frequencies))


    if "gtr" in type:
        pi_a = 0.32; pi_t = 0.34; pi_c = 0.16; pi_g = 0.18
        gtr_rates = [ 1.64390601,  1.27668478,  0.795571,  0.44377381,  0.32759197,  0.25651819]
        mu_dict = {'AG':gtr_rates[0]*pi_g, 'GA':gtr_rates[0]*pi_a, 'CT':gtr_rates[1]*pi_t, 'TC':gtr_rates[1]*pi_c, 'AC':gtr_rates[2]*pi_c, 'CA':gtr_rates[2]*pi_a, 'TG':gtr_rates[3]*pi_g, 'GT':gtr_rates[3]*pi_t, 'AT':gtr_rates[4]*pi_t, 'TA':gtr_rates[4]*pi_a, 'GC':gtr_rates[5]*pi_c, 'CG':gtr_rates[5]*pi_g}
    else:
        mu = 1.
        kappa = 4.0
        mu_dict = {'AT': mu, 'TA':mu, 'CG': mu, 'GC':mu, 'AC': mu, 'CA':mu, 'GT':mu, 'TG':mu, 'AG': kappa*mu, 'GA':kappa*mu, 'CT':kappa*mu, 'TC':kappa*mu}

    for i in range(len(true_frequencies)):
        freqs = true_frequencies[i]
        m = dNdS_from_MutSel(freqs, mu_dict)
        dn[i] = m.compute_dn()
        ds[i] = m.compute_ds()

    dnds_2rate = dn / ds
    dnds_1rate = dn / np.mean(ds)


    with open(dndsfile, "w") as f:
        f.write("site\tdnds2\tdnds1\tdn\tds\n")
        for i in range(len(dnds_2rate)):
            f.write(str(i+1) + "\t" + str(dnds_2rate[i]) + "\t" + str(dnds_1rate[i]) + "\t" + str(dn[i]) + "\t" + str(ds[i]) + "\n")


