# SJS
# Script to merge all dN/dS inferences into a single data frame, one for each simulation set, for simulations using real trees.
#args<-commandArgs(TRUE)
#if (length(args) != 1)
#{
#    stop("Supply the directory where results are stored as a cmd line argument.")
#}

require(dplyr)
require(readr)
require(stringr)
max_threshold = 9999.99 # Hyphy assigns this value (or greater, some decimal threshold lots of points out) to parameter upon failure to converge
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

#RESULTDIR <- args[1]
#if (str_sub(RESULTDIR, start=-1) != "/"){ RESULTDIR <- paste0(RESULTDIR, "/") }
RESULTDIR <- "~/Dropbox/dnds1rate2rate_data_results/results/realtrees_results/"
TRUEDIR <- "../simulation/"
numcol <- 100
datasets <- c("amine", "h3", "camelid", "vertrho", "hivrt")
pi <- "unequalpi"
types <- c("nobias", "bias")
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
                         "method"    = factor())   
levels(realdat$type) <- as.factor(types)
levels(realdat$method) <- c("FUBAR1", "FUBAR2")
 

for (type in types){
    print(type)
    # Load true simulated dN/dS values
    true <- read.csv(paste0(TRUEDIR, "truednds_", pi, "_", type, ".csv"))
    true_dnds <- true$dnds
    true_dn   <- true$dn
    true_ds   <- true$ds

    for (d in datasets){
        print(d)
        for (repl in 1:nreps){
            print(repl)
            # Read in raw results
            fubar1 <- read.csv(paste0(RESULTDIR, "rep", repl, "_", d, "_", pi, "_", type, "_FUBAR1.txt"))
            fubar2 <- read.csv(paste0(RESULTDIR, "rep", repl, "_", d, "_", pi, "_", type, "_FUBAR2.txt"))
            
            raw_fubar1_w = fubar1$beta / fubar1$alpha
            raw_fubar2_w = fubar2$beta / fubar2$alpha               
            fubar1_w = clean_dnds_divide( raw_fubar1_w, numcol )
            fubar2_w = clean_dnds_divide( raw_fubar2_w, numcol )            

            # create data frame per method  
            firstpart   <- data.frame("dataset" = rep(d, numcol), "site" = 1:numcol, "rep" = repl, "true" = true_dnds, "truedn" = true_dn, "trueds" = true_ds, "type" = rep(type, numcol) )    
            temp_fubar1 <- cbind(firstpart, data.frame("dnds" = fubar1_w, "dn" = fubar1$beta, "ds" = fubar1$alpha, "method" = rep("FUBAR1", numcol)))
            temp_fubar2 <- cbind(firstpart, data.frame("dnds" = fubar2_w, "dn" = fubar2$beta, "ds" = fubar2$alpha, "method" = rep("FUBAR2", numcol)))
            temp <- rbind(temp_fubar1, temp_fubar2)
            realdat <- rbind(realdat, temp)
        }
    }
}

# Save full dataset
write.csv(realdat, paste0(OUTDIR, "full_results_realtrees.csv"), quote=F, row.names=F)




