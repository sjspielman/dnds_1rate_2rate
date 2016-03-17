## Balanced trees
Balanced trees created with the following R code:
Trees are named as "n<exponent>_bl<branch_length>.tre", where 2^exponent = number of taxa and bl = branch length.
```r
require(ape)
branch_lengths <- c( 0.0025, 0.01, 0.04, 0.16, 0.64 )
for (exp in 7:11){
    t <- stree(2^exp, type = "balanced")
    for (bl in branch_lengths){
        t <- compute.brlen(t, bl)
        name <- paste("n", exp, "_bl", bl, ".tre", sep="")
        write.tree(t, name)
    }
}
```


## Real data trees
1. `camelid.tre`, `hivrt.tre`, and `vertrho.tre` are from Murrell et al. 2012, 2013. Note that rhodospin tree taken directly but trees for other two datasets made here with `FastTree -nt -gtr -nosupport`, since BL were not provided.  
2. `amine.tre` is from Spielman et al. 2015
3. `h3.tre` is the H3N2 hemagglutinin tree from Meyer Wilke 2015 (PLoS Path) 


Trees represent..
1. Vertebrate rhodopsin: small and high divergence
2. Camelid             : intermediate and high divergence
3. HIV RT              : intermediate and low divergence
4. H3N2 HA             : large and low divergence 
5. amine               : large and high divergence

Note that simulations for (4) and (5) will be analyzed only with FUBAR and SLAC for runtime considerations.

