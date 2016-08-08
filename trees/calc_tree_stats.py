## Calculate tree statistics, shown in Table 1 in manuscript

from dendropy import *
import numpy as np
import os

trees = ["amine", "vertrho", "camelid", "h3", "hivrt"]

for tree in trees:
    t = Tree.get_from_path('../phylogenies/' + dir + file, 'newick')

            # tree length
            tree_length = t.length()

            # root-to-tip
            treetips = t.leaf_nodes()
            rtt = []
            for tip in treetips:
                rtt.append( tip.distance_from_root() )

            # patristic
            pd = []
            dist = treecalc.PatristicDistanceMatrix(tree=t)
            for i, t1 in enumerate(t.taxon_set):
                for t2 in t.taxon_set[i:]:
                        d = dist(t1,t2)
                        pd.append( float(d) )


            outf.write(type + '\t' + name + '\t' + str(round(tree_length, 5)) + '\t' + str(round(np.mean(rtt), 5)) + '\t' + str(round(np.mean(pd), 5)) + '\n')
outf.close()
