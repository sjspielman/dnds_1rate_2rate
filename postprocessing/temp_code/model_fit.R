# SJS
# Script to merge dN/dS inferences for each dataset into a single data frame.

require(dplyr)
require(readr)
require(cowplot)

numcol <- 100
ntaxa <- 7:11
branch_lengths <- c(0.0025, 0.01, 0.04, 0.16, 0.64)
nreps <- 50


# Initialize results data frame
dat <- data.frame("type"     = factor(),
                  "ntaxa"     = numeric(),
                  "bl"    = numeric(),
                  "site"  = numeric(),
                  "rep"  = numeric(),
                  "fel1_logl"= numeric(),
                  "fel2_logl"= numeric())

levels(dat$type) <- as.factor(c("gtr", "bias_gtr"))

for (type in c("gtr", "bias_gtr")){
for (n in ntaxa){
    print(n)
    for (bl in branch_lengths){
        print(bl)
        for (repl in 1:nreps){

            fel1_raw <- paste0("../results/raw_results/rep", repl, "_n", n, "_bl", bl, "_", type, "_FEL1_GTR.txt")
            fel2_raw <- paste0("../results/raw_results/rep", repl, "_n", n, "_bl", bl, "_", type, "_FEL2_GTR.txt")
            fel1 <- read.csv(fel1_raw)
            fel2 <- read.csv(fel2_raw)        
            fel1_logl <- fel1$Log.L.
            fel2_logl <- fel2$Full.Log.L.
            #if (fel1_logl == 0 || fel2_logl == 0){ next }
                
            # create data frame per method  
            temp  <- data.frame("type" = type, "ntaxa" = rep(2^n, numcol), "bl" = rep(bl, numcol), "site" = 1:numcol, "rep" = repl, "fel1_logl" = fel1_logl, "fel2_logl" = fel2_logl )    
 
            # Merge
            dat <- rbind(dat, temp)
        }
    }
}}

dat %>% group_by(rep, type, ntaxa, bl) %>% 
        summarize(total_fel1_logl = sum(fel1_logl), total_fel2_logl = sum(fel2_logl) ) %>%
        mutate(logl_diff = 2*(total_fel2_logl - total_fel1_logl), 
               p         = (1 - pchisq(logl_diff, 100)),   # Right df? 
               sig       = (p < 0.01),                     # Correct p-value? 
               twobetter = (total_fel2_logl > total_fel1_logl)) %>%
               select(-total_fel2_logl, -total_fel2_logl, -p, -logl_diff) %>% 
        ungroup() %>% group_by(type, ntaxa, bl) %>% tally(sig) %>% 
        spread(type, n) %>% print.data.frame()
# the columns here are the number of replicates (of 50!) which are significantly in favor of two-rate model via LRT
#    ntaxa     bl gtr bias_gtr
# 1    128 0.0025   0        0
# 2    128 0.0100   0        0
# 3    128 0.0400   2       43
# 4    128 0.1600  14       50
# 5    128 0.6400  44       50
# 6    256 0.0025   0        0
# 7    256 0.0100   0        3
# 8    256 0.0400   8       50
# 9    256 0.1600  26       50
# 10   256 0.6400  50       50
# 11   512 0.0025   0        0
# 12   512 0.0100   3       49
# 13   512 0.0400  23       50
# 14   512 0.1600  46       50
# 15   512 0.6400  50       50
# 16  1024 0.0025   0        3
# 17  1024 0.0100  10       50
# 18  1024 0.0400  41       50
# 19  1024 0.1600  49       50
# 20  1024 0.6400  50       50
# 21  2048 0.0025   2       42
# 22  2048 0.0100  28       50
# 23  2048 0.0400  50       50
# 24  2048 0.1600  50       50
# 25  2048 0.6400  50       50        
#         
        
        
# There are 100 tests per dataset, so with Bonferroni we have p < 1e-4
alpha = 1e-4
dat %>% mutate(logl_diff = 2*(fel2_logl - fel1_logl), p = (1 - pchisq(logl_diff, 1)), sig = (p < alpha), twobetter = (fel2_logl > fel1_logl)) %>% select(-logl_diff, -p, -fel1_logl, -fel2_logl) %>% group_by(type, ntaxa, bl, site) %>% summarize(sumsig = sum(sig)) -> site.sig
# All significant differences between sites are in favor of a two-rate model


# 
# 
# 
# 
# 
# dat.lrt$better <- "onerate"
# dat.lrt$better[dat.lrt$fel2_logl > dat.lrt$fel1_logl] <- "tworate"
# dat.lrt %>% filter(sig == TRUE, better == "tworate", type == "gtr") %>% group_by(ntaxa,bl, site) %>% tally() -> tally.tworate.better.nobias
# 
# dat.lrt %>% filter(sig == TRUE, better == "tworate", type == "bias_gtr") %>% group_by(ntaxa,bl, site) %>% tally() -> tally.tworate.better.bias
# 
# 
# gamma <- read.table("../simulation/codonbias_gtr_term.txt", header=T)
# true_bias <- read.table("../simulation/codon_freq_lib_bias_gtr_true_dnds.txt", header=T)
# 
# dat.lrt %>% group_by(type, ntaxa, bl, site) %>% tally(sig) %>% filter(type == "bias_gtr") %>% left_join(true_bias) -> bias.lrt
# 
# 
# 
