# SJS
# Script to merge all dN/dS inferences into a single data frame, one for each simulation set.

args<-commandArgs(TRUE)
if (length(args) != 1)
{
     stop("Supply the directory where results are stored as a cmd line argument.")
}


require(dplyr)
require(readr)
require(stringr)
max_threshold = 9999     # Hyphy assigns this value (or greater, some decimal threshold lots of points out) to parameter upon failure to converge
zero_threshold = 1e-15   # My zero threshold
OUTDIR="dataframes/"

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
        else if (dn >= max_threshold | ds >= max_threshold |  dn == Inf | ds ==  Inf | ds <= zero_threshold) #no convergence
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


RESULTDIR <- args[1]
if (str_sub(RESULTDIR, start=-1) != "/"){ RESULTDIR <- paste0(RESULTDIR, "/") }
TRUEDIR <- "../simulation/"
numcol <- 100
ntaxa <- 7:11
pi.types <- c( "unequalpi")
types <- c("bias")
branch_lengths <- c(0.0025, 0.01, 0.04, 0.16, 0.64)
nreps <- 50
numrow <- 750000  # Rows per dataframe, calc'd as 5*5*100*50*6 = ntaxa * bl * alnlen * reps * methods



for (pi in pi.types){

    for (type in types)
    {
    
        true <- read.csv(paste0(TRUEDIR, "truednds_", pi, "_", type, ".csv"))
        true_dnds  <- true$dnds
        true_dn    <- true$dn
        true_ds    <- true$ds
    
        # Initialize results data frame
        df.results <- data.frame("ntaxa"     = rep(0, numrow),        # number of taxa
                                 "bl"        = rep(0, numrow),        # branch lengths
                                 "site"      = rep(0, numrow),        # site index (1:numcol)
                                 "rep"       = rep(0, numrow),        # replicate (1-50)
                                 "true"      = rep(0, numrow),        # True dN/ddS
                                 "truedn"    = rep(0, numrow),        # True dN
                                 "trueds"    = rep(0, numrow),        # True dS
                                 "type"      = rep(type, numrow),     # nobias, bias
                                 "dnds"      = rep(0, numrow),        # inferred dN/dS
                                 "dn"        = rep(0, numrow),        # inferred dN (If 1-rate method, then same as dN/dS).
                                 "ds"        = rep(0, numrow),        # inferred dS (If 1-rate method, then 1, or mean(dS) if SLAC).
                                 "method"    = rep("SLAC1", numrow))  # Inference method
        levels(df.results$type) <- as.factor(type)
        levels(df.results$method) <- c("FEL1", "FEL2", "FUBAR1", "FUBAR2", "SLAC1", "SLAC2")

        i <- 1
        for (n in ntaxa){
            print(n)
            for (bl in branch_lengths){
                print(bl)
                for (repl in 1:nreps){

                    # Read in raw results
                    fel1   <- read.csv(paste(RESULTDIR, "rep", repl, "_n", n, "_bl", bl, "_", pi, "_", type, "_FEL1.txt", sep=""))        
                    fel2   <- read.csv(paste(RESULTDIR, "rep", repl, "_n", n, "_bl", bl, "_", pi, "_", type, "_FEL2.txt", sep=""))        
                    slac   <- read.table(paste(RESULTDIR, "rep", repl, "_n", n, "_bl", bl, "_", pi, "_", type, "_SLAC.txt", sep=""), header=T)
                    fubar1 <- read.csv(paste(RESULTDIR, "rep", repl, "_n", n, "_bl", bl, "_", pi, "_", type, "_FUBAR1.txt", sep=""))
                    fubar2 <- read.csv(paste(RESULTDIR, "rep", repl, "_n", n, "_bl", bl, "_", pi, "_", type, "_FUBAR2.txt", sep=""))
                
                    # Clean up dN/dS values to replace uninformative with NA
                    fel1_w = clean_dnds_fel1(fel1, numcol)
                    fel2_w = clean_dnds_fel2(fel2, numcol)
                    slac1_w = slac$dN/(mean(slac$dS))
                    raw_slac2_w  = slac$dN/slac$dS
                    slac2_w = clean_dnds_divide( raw_slac2_w, numcol )
                
                    raw_fubar1_w = fubar1$beta / fubar1$alpha
                    raw_fubar2_w = fubar2$beta / fubar2$alpha               
                    fubar1_w = clean_dnds_divide( raw_fubar1_w, numcol )
                    fubar2_w = clean_dnds_divide( raw_fubar2_w, numcol )            

                    # create data frame per method  
                    firstpart   <- data.frame("ntaxa" = rep(2^n, numcol), "bl" = rep(bl, numcol), "site" = 1:numcol, "rep" = repl, "true" = true_dnds, "truedn" = true_dn, "trueds" = true_ds, "type" = rep(type, numcol) )    
                    temp_fel1   <- cbind(firstpart, data.frame("dnds" = fel1_w, "dn" = fel1_w, "ds" = 1., "method" = rep("FEL1", numcol)))
                    temp_fel2   <- cbind(firstpart, data.frame("dnds" = fel2_w, "dn" = fel2$dN, "ds" = fel2$dS, "method" = rep("FEL2", numcol)))
                    temp_slac1  <- cbind(firstpart, data.frame("dnds" = slac1_w, "dn" = slac$dN, "ds" = mean(slac$dS), "method" = rep("SLAC1", numcol)))
                    temp_slac2  <- cbind(firstpart, data.frame("dnds" = slac2_w, "dn" = slac$dN, "ds" = slac$dS, "method" = rep("SLAC2", numcol)))
                    temp_fubar1 <- cbind(firstpart, data.frame("dnds" = fubar1_w, "dn" = fubar1$beta, "ds" = fubar1$alpha, "method" = rep("FUBAR1", numcol)))
                    temp_fubar2 <- cbind(firstpart, data.frame("dnds" = fubar2_w, "dn" = fubar2$beta, "ds" = fubar2$alpha, "method" = rep("FUBAR2", numcol)))

                    # Merge
                    temp <- rbind(temp_fel1, temp_fel2, temp_slac1, temp_slac2, temp_fubar1, temp_fubar2)
                    df.results[i:(i+nrow(temp)-1),] <- temp
                    i <- i + nrow(temp)
                }
            }
        }
        # Save full dataset
        write_csv(df.results, paste0(OUTDIR, "full_results_", pi, "_", type, ".csv"))
    }
}


