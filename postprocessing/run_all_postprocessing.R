# SJS
# Post-processing full pipeline.
# NOTE: If you only want to re-create plots, you might want to comment out the first two lines (those calling scripts process_raw_results.R and summarize_results.R) to avoid unnecessary stuff.

source("process_raw_results.R") # Extracts dN/dS inferences from raw FEL, SLAC, and FUBAR inference files
source("summarize_results.R")   # Summarize dN/dS, dN, and dS estimates to compute, for each quantity, correlation, estimator bias, and RMSD with true values. The resulting statistics are averaged across replicates. 
source("plot_results.R")