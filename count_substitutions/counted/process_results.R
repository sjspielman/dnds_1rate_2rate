# SJS
# Script to merge dN/dS inferences for each dataset into a single data frame.

require(dplyr)
require(readr)

max_threshold = 1e2 # if dN/dS greater than this, clearly we're in an uninformative situation. Also, 5e-9 is HyPhy's bailed 0.

# Clean dN/dS values for those methods which return separate dN and dS (hence divide)
# Methods include SLAC, FUBAR1, FUBAR2
clean_dnds_divide <- function(w, numcol)
{
    for (i in 1:numcol)
    {
        if (is.na(w[i]))
        {
            w[i] <- NA
        }
        else if (w[i] >= max_threshold | w[i] == Inf)
        {
            w[i] <- NA
        }
    }
    w
}

# Clean dN/dS values for fel1 and fel2 (detect NAs)
clean_dnds_fel <- function(df.fel, numcol)
{   
    w <- df.fel$dN.dS
    for (i in 1:numcol)
    {
        if (df.fel$LRT[i] == 0 & df.fel$p.value[i] == 1 & (w[i] == 0 | w[i] == 1)) #uninformative
        {
            w[i] <- NA
        }
        else if (w[i] >= max_threshold | w[i] == Inf)
        {
            w[i] <- NA
        }
    }
    w
} 
numcol <- 100
ntaxa <- 7:11
branch_lengths <- c(0.0025, 0.01, 0.04, 0.16, 0.64)
nreps <- 50


# Initialize results data frame
dat <- data.frame("ntaxa"     = numeric(),
                         "bl"    = numeric(),
                         "site"  = numeric(),
                         "rep"  = numeric(),
                         "ncount"= numeric(),
                         "scount"   = numeric())
i <- 1
for (n in ntaxa){
    print(n)
    for (bl in branch_lengths){
        print(bl)
        for (repl in 1:nreps){

            name <- paste0("rep", repl, "_n", n, "_bl", bl, "_bias_gtr_counted.txt")
            counted <- read.table(name, header=T)        
            ns = counted$ns_changes
            ss = counted$s_changes

            # create data frame per method  
            temp  <- data.frame("ntaxa" = rep(2^n, numcol), "bl" = rep(bl, numcol), "site" = 1:numcol, "rep" = repl, "ncount" = ns, "scount" = ss )    
 
            # Merge
            dat <- rbind(dat, temp)
        }
    }
}

dat %>% group_by(ntaxa, bl, rep) %>% mutate(ratio_ns = ncount/scount) %>% na.omit() %>% summarize(mean_ratio = mean(ratio_ns)) %>% filter(!is.infinite(mean_ratio)) -> ns_count_ratio

dat %>% group_by(ntaxa, bl, rep) %>% mutate(ratio_ns = ncount/scount) %>% na.omit() %>% filter(!is.infinite(ratio_ns)) %>% summarize(mean_ratio = mean(ratio_ns))  -> ns_count_ratio

ggplot(ns_count_ratio, aes(x = as.factor(bl), y = mean_ratio, fill = as.factor(ntaxa))) + geom_boxplot() + geom_hline(yintercept=1)
