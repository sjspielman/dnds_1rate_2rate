import sys
from function_library import *


types = ["equalpi_nobias", "equalpi_bias", "unequalpi_nobias", "unequalpi_bias"]
for type in types:

    print "Calculuting true dN/dS for ", type, "simulations"
    
    freqfile = "codon_freq_lib_" + type + ".txt"
    dndsfile = "truednds_" + type + ".csv"

    true_frequencies = np.loadtxt(freqfile)

    dn = np.zeros(len(true_frequencies))
    ds = np.zeros(len(true_frequencies))


    if type.startswith("equalpi"):
        pi_a = 0.25; pi_t = 0.25; pi_c = 0.25; pi_g = 0.25;
    elif type.startswith("unequalpi"):
        pi_a = 0.32; pi_t = 0.28; pi_c = 0.18; pi_g = 0.22;
    kappa = 4.0
    mu_dict = {'AT': pi_t, 'TA':pi_a, 'CG': pi_g, 'GC':pi_c, 'AC': pi_c, 'CA':pi_a, 'GT':pi_t, 'TG':pi_g, 'AG': kappa*pi_g, 'GA':kappa*pi_a, 'CT':kappa*pi_t, 'TC':kappa*pi_c}
    print mu_dict

    for i in range(len(true_frequencies)):
        freqs = true_frequencies[i]
        m = dNdS_from_MutSel(freqs, mu_dict)
        dn[i] = m.compute_dn()
        ds[i] = m.compute_ds()

    dnds = dn / ds


    with open(dndsfile, "w") as f:
        f.write("site,dnds,dn,ds\n")
        for i in range(len(dnds)):
            f.write(str(i+1) + "," + str(dnds[i]) + "," + str(dn[i]) + "," + str(ds[i]) + "\n")


