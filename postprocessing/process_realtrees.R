# SJS
# Script to merge all dN/dS inferences into a single data frame, one for each simulation set, for simulations using real trees.

require(dplyr)
require(readr)

max_threshold = 9999.99 # Hyphy assigns this value (or greater, some decimal threshold lots of points out) to parameter upon failure to converge


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


# Clean dN/dS values for fel1, fubar1
clean_dnds_fel1 <- function(df.fel, numcol)
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

# Clean dN/dS values for fel2
clean_dnds_fel2 <- function(df.fel, numcol)
{   
    dn <- df.fel$dN
    ds <- df.fel$dS
    w <- c()
    for (i in 1:numcol)
    {
        k <- TRUE
        if (df.fel$LRT[i] == 0 & df.fel$p.value[i] == 1) #uninformative
        {
            k <- FALSE
        }
        else if (dn >= max_threshold | ds >= max_threshold |  dn == Inf | ds ==  Inf)
        {
            k <- FALSE
        }
        if (k){ 
            w <- c(w, dn[i]/ds[i]) 
        }
        else{
            w <- c(w, NA)
        }
    }
    w
} 


RESULTDIR <- "../results/results_realtrees/"
TRUEDIR <- "../simulation/"
numcol <- 100
datasets <- c("amine", "h3", "camelid", "vertrho", "hivrt")
types <- c("gtr", "bias_gtr")
nreps <- 50

# Initialize results data frame
realdat <- data.frame("dataset"   = character(),
                         "site"      = numeric(),  # site index (1:numcol)
                         "rep"       = numeric(),
                         "true"     = numeric(),  # True dN/dS
                         "truedn"    = numeric(), # True dN
                         "trueds"    = numeric(), # True dS
                         "type"      = factor(),   # nobias, bias
                         "dnds"      = numeric(),  # inferred dN/dS
                         "dn"        = numeric(),   # inferred dN (If 1-rate method, then same as dN/dS). NA for FEL2_1, FUBAR2_1
                         "ds"        = numeric(),   # inferred dS (If 1-rate method, then 1, or mean(dS) if SLAC). NA for FEL2_1, FUBAR2_1
                         "method"    = factor())   # pvalue or pp to indicate which type of value is in the signif column
levels(realdat$type) <- as.factor(types)
levels(realdat$method) <- c("SLAC1", "SLAC2")  #"FUBAR1", "FUBAR2", 
 

for (type in types){
    print(type)
    # Load true simulated dN/dS values
    true <- read.table(paste0(TRUEDIR,"codon_freq_lib_", type, "_true_dnds.txt"), header=T)
    true_dnds1 <- true$dnds1
    true_dnds2 <- true$dnds2
    true_dn    <- true$dn
    true_ds    <- true$ds

    for (d in datasets){
        print(d)
        for (repl in 1:nreps){
            print(repl)
            # Read in raw results
            slac   <- read.table(paste(RESULTDIR, "rep", repl, "_", d, "_", type, "_SLAC_GTR.txt", sep=""), header=T)
            #fubar1 <- read.csv(paste(RESULTDIR, "rep", repl, "_", d, "_", type, "_FUBAR1.txt", sep=""))
            #fubar2 <- read.csv(paste(RESULTDIR, "rep", repl, "_", d, "_", type, "_FUBAR2.txt", sep=""))
            
            # Clean up dN/dS values to replace uninformative with NA
            slac1_w = slac$dN/(mean(slac$dS))
            raw_slac2_w  = slac$dN/slac$dS
            slac2_w = clean_dnds_divide( raw_slac2_w, numcol )
            
            #raw_fubar1_w = fubar1$beta / fubar1$alpha
            #raw_fubar2_w = fubar2$beta / fubar2$alpha               
            #fubar1_w = clean_dnds_divide( raw_fubar1_w, numcol )
            #fubar2_w = clean_dnds_divide( raw_fubar2_w, numcol )            

            # create data frame per method  
            firstpart   <- data.frame("dataset" = rep(d, numcol), "site" = 1:numcol, "rep" = repl, "true1" = true_dnds1, "true2" = true_dnds2, "truedn" = true_dn, "trueds" = true_ds, "type" = rep(type, numcol) )    
            temp_slac1  <- cbind(firstpart, data.frame("dnds" = slac1_w, "dn" = slac$dN, "ds" = mean(slac$dS), "method" = rep("SLAC1", numcol)))
            temp_slac2  <- cbind(firstpart, data.frame("dnds" = slac2_w, "dn" = slac$dN, "ds" = slac$dS, "method" = rep("SLAC2", numcol)))
            #temp_fubar1 <- cbind(firstpart, data.frame("dnds" = fubar1_w, "dn" = fubar1$beta, "ds" = fubar1$alpha, "method" = rep("FUBAR1", numcol)))
            #temp_fubar2 <- cbind(firstpart, data.frame("dnds" = fubar2_w, "dn" = fubar2$beta, "ds" = fubar2$alpha, "method" = rep("FUBAR2", numcol)))
            temp <- rbind(temp_slac1, temp_slac2)  #, temp_fubar1, temp_fubar2)
            realdat <- rbind(realdat, temp)
        }
    }
}

# Save full dataset
write.csv(realdat, "full_results_realtrees.csv", quote=F, row.names=F)




