# SJS
# Post-processing full pipeline.
# NOTE: If you only want to re-create plots, you might want just run the "plot_results.R" script to save some time.

require(dplyr)
require(tidyr)
require(readr)
require(cowplot)


source("summary_functions.R")        # Load functions required for summarizing
source("process_raw_results.R")      # Extracts dN/dS inferences from raw FEL, SLAC, and FUBAR inference files
source("process_realtrees.R")        # Extracts dN/dS inferences from raw SLAC, and FUBAR inference files for simulations performed along real trees
source("summarize_results.R")        # Summarize dN/dS, dN, and dS estimates to compute, for each quantity, correlation, estimator bias, and RMSD with true values. The resulting statistics are averaged across replicates. 
source("process_gtr_bias_counted.R") # Process counted number of changes for bias simulations
#source("plot_results.R")