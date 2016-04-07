# SJS
# Build and save summary statistic dataframes for dN/dS inference. For each parameter dN/dS, dN, and dS, this script computes the Pearson correlation, estimator bias, and RMSD for all inferences. Note that RMSD and estimator bias for SLAC dN and dS parameter estimates are meaningless due to scaling, but they are processed anyways.

require(dplyr)
require(tidyr)
require(readr)
source("summary_functions.R")

for (pi in c("equalpi", "unequalpi"))
{
    nobias <- read_csv(paste0("full_results_", pi, "_nobias.csv"))
    bias   <- read_csv(paste0("full_results_", pi, "_bias.csv"))

    nobias.dnds.sum   <- summarize_dnds(nobias, "nobias")
    bias.dnds.sum     <- summarize_dnds(bias, "bias")
    dnds.sum          <- rbind(nobias.dnds.sum, bias.dnds.sum)
    dnds.sum$type     <- factor(dnds.sum$type, levels=c("nobias", "bias"))


    nobias.dn.sum  <- summarize_dn(nobias, "nobias")
    nobias.ds.sum  <- summarize_ds(nobias, "nobias")
    bias.dn.sum    <- summarize_dn(bias, "bias")
    bias.ds.sum    <- summarize_ds(bias, "bias")
    dn.ds.sum      <- rbind(nobias.ds.sum, nobias.dn.sum, bias.ds.sum, bias.dn.sum)
    dn.ds.sum$type <- factor(dn.ds.sum$type, levels=c("nobias", "bias"))

    write_csv(dnds.sum, paste0("dnds_summary_", pi, ".csv"))
    write_csv(dn.ds.sum, paste0("dn_ds_summary_", pi, ".csv"))

}
