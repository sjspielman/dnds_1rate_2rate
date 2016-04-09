# SJS
# Build and save summary statistic dataframes for dN/dS inference. For each parameter dN/dS, dN, and dS, this script computes the Pearson correlation, estimator bias, and RMSD for all inferences. Note that RMSD and estimator bias for SLAC dN and dS parameter estimates are meaningless due to scaling, but they are processed anyways.

require(dplyr)
require(tidyr)
require(readr)
source("summary_functions.R")

for (pi in c("equalpi", "unequalpi"))
{
    for (type in c("nobias", "bias"))
    {
        full <- read_csv(paste0("full_results_", pi, "_", type, ".csv"))

        dnds.sum   <- summarize_dnds(full, type)
        dn.sum     <- summarize_dn(full, type)
        ds.sum     <- summarize_ds(full, type)
        dn.ds.sum <- rbind(dn.sum, ds.sum)
    
        write_csv(dnds.sum, paste0("dnds_summary_", pi, "_", type, ".csv"))
        write_csv(dn.ds.sum, paste0("dn_ds_summary_", pi, "_", type, ".csv"))
    }
}
