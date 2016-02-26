# SJS
# Script to merge dN/dS inference, simulation root rank results into a single data frame. 
# Saves one file per data type -  default (assumptions all met), codon bias, and mutational asymmetry.
# Additionally creates and saves a summary data frame with correlations and estimator bias for each simulation.
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


# Function to create summary data frame with information about mean correlation and estimator bias
summarize_results <- function(dat, type)
{

    # Summarize, for true1
    dat %>% select(-true2) %>% na.omit() %>% group_by(ntaxa, bl, method, rep) %>% 
    do( bias.raw = glm(dnds ~ offset(true1), dat = .),  
        cor.raw  = cor(.$true1, .$dnds),
        rmsd.raw = sqrt(mean((.$true1 - .$dnds)^2))) %>% 
    mutate(estbias1 = summary(bias.raw)$coeff[1], r1 = cor.raw[1], rmsd1 = rmsd.raw[[1]]) %>% 
    select(-bias.raw, -cor.raw, -rmsd.raw) %>% na.omit() %>% print.data.frame()-> summ1
    
    # Summarize, for true2
    dat %>% select(-true1) %>% na.omit() %>% group_by(ntaxa, bl, method, rep) %>% 
    do( bias.raw = glm(dnds ~ offset(true2), dat = .),  
        cor.raw  = cor(.$true2, .$dnds),
        rmsd.raw = sqrt(mean((.$true2 - .$dnds)^2))) %>% 
    mutate(estbias2 = summary(bias.raw)$coeff[1], r2 = cor.raw[1], rmsd2 = rmsd.raw[[1]]) %>% 
    select(-bias.raw, -cor.raw, -rmsd.raw) %>% na.omit() %>% print.data.frame()-> summ2
    
    summ.full <- left_join(summ1, summ2)
    summ.full$type <- type
    summ.full
    
}

# Initialize summary data frame. Won't contain too many rows, so just rbind it
df.sum <- data.frame("ntaxa"    = factor(),   # number of taxa
                     "bl"       = factor(),   # branch lengths
                     "method"   = factor(),   # inference method
                     "type"     = factor(),   # default, bias, asym
                     "r1"        = numeric(), # correlation with true1
                     "estbias1"  = numeric(), # estimator bias with true1
                     "rmsd1"     = numeric(), # RMSD with true1
                     "r2"        = numeric(), # correlation with true2
                     "estbias2"  = numeric(), # estimator bias with true2                  
                     "rmsd2"     = numeric()) # RMSD with true2    


DATADIR <- "../results/"
numcol <- 100
ntaxa <- 7:11
types <- c("nobias" , "bias", "asym")
branch_lengths <- c(0.0025, 0.01, 0.04, 0.16, 0.64)
nreps <- 50
numrow <- 1e6  # Rows per dataframe, calc'd as 100*5*5*6*50 = alnlen * ntaxa * bl * methods * reps



for (type in types){
    
    # Load true simulated dN/dS values
    true <- read.table(paste0(DATADIR,"codon_freq_lib_", type, "_info.txt"), header=T)
    true_dnds1 <- true$dnds1
    true_dnds2 <- true$dnds2
    
    # Initialize results data frame
    df.results <- data.frame("ntaxa"     = rep(0, numrow), # number of taxa
                             "bl"        = rep(0, numrow),  # branch lengths
                             "site"      = rep(0, numrow),  # site index (1:numcol)
                             "rep"       = rep(0, numrow),
                             "true1"      = rep(0, numrow),  # True dN/dS, calculated in 1-rate framework (note that, for no bias, this is same as true2)
                             "true2"      = rep(0, numrow),  # True dN/dS, calculated in 2-rate framework 
                             "type"      = rep(as.factor(type), numrow),   # nobias, bias, asym
                             "dnds"      = rep(0, numrow),  # inferred dN/dS
                             "dn"        = rep(0, numrow),   # inferred dN (If 1-rate method, then same as dN/dS). NA for FEL2_1, FUBAR2_1
                             "ds"        = rep(0, numrow),   # inferred dS (If 1-rate method, then 1, or mean(dS) if SLAC). NA for FEL2_1, FUBAR2_1
                             "method"    = rep("slac", numrow))   # pvalue or pp to indicate which type of value is in the signif column
    levels(df.results$type) <- as.factor(type)
    levels(df.results$method) <- c("FEL1", "FEL2", "FEL2_1", "FUBAR1", "FUBAR2", "FUBAR2_1","SLAC1", "SLAC2")

    i <- 1
    for (n in ntaxa){
        print(n)
        for (bl in branch_lengths){
            print(bl)
            for (repl in 1:nreps){

                # Read in raw results
                fel1   <- read.csv(paste(DATADIR, "raw_results/rep", repl, "_n", n, "_bl", bl, "_", type, "_FEL1.txt", sep=""))        
                fel2   <- read.csv(paste(DATADIR, "raw_results/rep", repl, "_n", n, "_bl", bl, "_", type, "_FEL2.txt", sep=""))        
                slac   <- read.table(paste(DATADIR, "raw_results/rep", repl, "_n", n, "_bl", bl, "_", type, "_SLAC.txt", sep=""), header=T)
                fubar1 <- read.csv(paste(DATADIR, "raw_results/rep", repl, "_n", n, "_bl", bl, "_", type, "_FUBAR1.txt", sep=""))
                fubar2 <- read.csv(paste(DATADIR, "raw_results/rep", repl, "_n", n, "_bl", bl, "_", type, "_FUBAR2.txt", sep=""))
                
                # Clean up dN/dS values to replace uninformative with NA
                fel1_w = clean_dnds_fel(fel1, numcol)
                fel2_w = clean_dnds_fel(fel2, numcol)
                slac1_w = slac$dN/(mean(slac$dS))
                raw_slac2_w  = slac$dN/slac$dS
                slac2_w = clean_dnds_divide( raw_slac2_w, numcol )
                
                raw_fubar1_w = fubar1$beta / fubar1$alpha
                raw_fubar2_w = fubar2$beta / fubar2$alpha               
                fubar1_w = clean_dnds_divide( raw_fubar1_w, numcol )
                fubar2_w = clean_dnds_divide( raw_fubar2_w, numcol )

                fel21_w = fel2$dN / mean(fel2$dS)
                fubar21_w = fubar2$beta / mean(fubar2$alpha)
            

                # create data frame per method  
                firstpart <-   data.frame("ntaxa" = rep(2^n, numcol), "bl" = rep(bl, numcol), "site" = 1:numcol, "rep" = repl, "true1" = true_dnds1, "true2" = true_dnds2, "type" = rep(type, numcol) )    
                temp_fel1 <-   cbind(firstpart, data.frame("dnds" = fel1_w, "dn" = fel1_w, "ds" = 1., "method" = rep("FEL1", numcol)))
                temp_fel2 <-   cbind(firstpart, data.frame("dnds" = fel2_w, "dn" = fel2$dN, "ds" = fel2$dS, "method" = rep("FEL2", numcol)))
                temp_slac1 <- cbind(firstpart, data.frame("dnds" = slac1_w, "dn" = slac$dN, "ds" = slac$dS, "method" = rep("SLAC1", numcol)))
                temp_slac2 <- cbind(firstpart, data.frame("dnds" = slac2_w, "dn" = slac$dN, "ds" = mean(slac$dS), "method" = rep("SLAC2", numcol)))
                temp_fubar1 <- cbind(firstpart, data.frame("dnds" = fubar1_w, "dn" = fubar1$beta, "ds" = fubar1$alpha, "method" = rep("FUBAR1", numcol)))
                temp_fubar2 <- cbind(firstpart, data.frame("dnds" = fubar2_w, "dn" = fubar2$beta, "ds" = fubar2$alpha, "method" = rep("FUBAR2", numcol)))
  
                temp_fel21 <-   cbind(firstpart, data.frame("dnds" = fel21_w, "dn" = NA, "ds" = NA, "method" = rep("FEL2_1", numcol)))
                temp_fubar21 <-   cbind(firstpart, data.frame("dnds" = fubar21_w, "dn" = NA, "ds" = NA, "method" = rep("FUBAR2_1", numcol)))


                # Merge
                temp <- rbind(temp_fel1, temp_fel2, temp_fel21, temp_slac1, temp_slac2, temp_fubar1, temp_fubar2, temp_fubar21)
                df.results[i:(i+nrow(temp)-1),] <- temp
                i <- i + nrow(temp)
                
            }
        }
    }

    # Save full dataset
    write_csv(df.results, paste0(DATADIR, "processed_results/full_results_", type, "_dataset.csv"))
    
    # Build up summary dataframe with correlations, estimator bias, percent sites estimated
    df.sum <- rbind(df.sum, summarize_results(df.results, type) )
}

write_csv(df.sum, paste0(DATADIR, "processed_results/summarized_results.csv"))





