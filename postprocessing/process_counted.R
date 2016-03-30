# SJS
args<-commandArgs(TRUE)
if (length(args) != 1)
{
    stop("Supply the directory where counted results are stored as a cmd line argument.")
}

require(dplyr)
require(readr)
require(stringr)


RESULTDIR <- args[1]
if (str_sub(RESULTDIR, start=-1) != "/"){ RESULTDIR <- paste0(RESULTDIR, "/") }

ntaxa <- 7:11
branch_lengths <- c(0.0025, 0.01, 0.04, 0.16, 0.64)
nreps <- 50
mutypes <- c("hky", "gtr")
types <- c("nobias", "bias")


for (mutype in mutypes){

    for (type in types){
        # Initialize results data frame
        dat <- data.frame("ntaxa"  = numeric(),
                          "bl"     = numeric(),
                          "site"   = numeric(),
                          "rep"    = numeric(),
                          "type"   = character(),
                          "ncount" = numeric(),
                          "scount" = numeric())
        for (n in ntaxa){
            print(n)
            for (bl in branch_lengths){
                print(bl)
                for (repl in 1:nreps){

                    name <- paste0(RESULTDIR,"rep", repl, "_n", n, "_bl", bl, "_", mutype, "_", type, "_counted.txt")
                    counted <- read.table(name, header=T)        
                    ns = counted$ns_changes
                    ss = counted$s_changes

                    # create data frame per method  
                    temp  <- data.frame("ntaxa" = rep(2^n, 100), "bl" = rep(bl, 100), "site" = 1:100, "rep" = repl, "type" = type, "ncount" = ns, "scount" = ss )    
 
                    # Merge
                    dat <- rbind(dat, temp)
                }
            }
        }
        write_csv(dat, paste0("substitution_counts_", mutype, "_", type, ".csv"))
    }
}

