# Quick script to extract mean optimized branch length for each simulation condition, specifically \Pi_unequal
# Note that all methods use ~same nucleotide fit, so these SLAC should extend out to all other methds.

import re
import numpy as np

def extract_bl(rep, n, bl, pitype, bias, path):
    file = path + "rep" + str(rep) + "_n" + str(n) + "_bl" + str(bl) + "_" + pitype + "_" + bias + "_SLAC_nucfit.txt"
    branches = []
    with open(file, "rU") as f:
        for line in f:
            find_bl = re.search(r"givenTree\.\w+\.t=(\d+\.\d+);$", line)
            if find_bl:
                branches.append( float(find_bl.group(1)))

    branches = np.array(branches)
    return str(rep)+","+str(2**n)+","+str(bl)+","+biastype + "," + str(np.mean(branches)) + "\n"

outfile = "dataframes/optimized_bl.csv"
with open(outfile, "w") as f:
    f.write("rep,ntaxa,bl,biastype,meanbl\n")
directory = "../results/unequalpi_nucfits/"
for biastype in ["nobias", "bias"]:
    for rep in range(1,51):
        for bl in [0.0025, 0.01, 0.04, 0.16, 0.64]:
            for n in range(7,12):
                s = extract_bl(rep,n,bl,"unequalpi", biastype, directory)
                with open(outfile, "a") as f:
                    f.write(s)
