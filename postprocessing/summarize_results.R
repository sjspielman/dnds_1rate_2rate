# SJS
# Build and save summary statistic dataframes for dN/dS inference. For each parameter dN/dS, dN, and dS, this script computes the Pearson correlation, estimator bias, and RMSD for all inferences. Note that RMSD and estimator bias for SLAC dN and dS parameter estimates are meaningless due to scaling, but they are processed anyways.

require(dplyr)
require(tidyr)
require(readr)


RESULTDIR <- "../results/processed_results/"   
nobias <- read_csv(paste0(RESULTDIR, "full_results_gtr_dataset.csv"))
bias   <- read_csv(paste0(RESULTDIR, "full_results_bias_gtr_dataset.csv"))

nobias.dnds.sum   <- summarize_dnds(nobias, "nobias")
bias.dnds.sum     <- summarize_dnds(bias, "bias")
dnds.sum          <- rbind(nobias.dnds.sum, bias.dnds.sum)
dnds.sum$truetype <- factor(dnds.sum$truetype, levels=c("true1", "true2"))
dnds.sum$type     <- factor(dnds.sum$type, levels=c("nobias", "bias"))


nobias.dn.sum  <- summarize_dn(nobias, "nobias")
nobias.ds.sum  <- summarize_ds(nobias, "nobias")
bias.dn.sum    <- summarize_dn(bias, "bias")
bias.ds.sum    <- summarize_ds(bias, "bias")
dn.ds.sum      <- rbind(nobias.ds.sum, nobias.dn.sum, bias.ds.sum, bias.dn.sum)
dn.ds.sum$type <- factor(dn.ds.sum$type, levels=c("nobias", "bias"))

write_csv(dnds.sum, "dnds_summary.csv")
write_csv(dn.ds.sum, "dn_ds_summary.csv")

