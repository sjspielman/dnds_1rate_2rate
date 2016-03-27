# SJS
# Build and save summary statistic dataframes for dN/dS inference. For each parameter dN/dS, dN, and dS, this script computes the Pearson correlation, estimator bias, and RMSD for all inferences. Note that RMSD and estimator bias for SLAC dN and dS parameter estimates are meaningless due to scaling, but they are processed anyways.

require(dplyr)
require(tidyr)
require(readr)


nobias <- read_csv("full_results_gtr_dataset.csv")
bias   <- read_csv("full_results_bias_gtr_dataset.csv")

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

write_csv(dnds.sum, "dnds_summary.csv")
write_csv(dn.ds.sum, "dn_ds_summary.csv")


#  Random code dump!
# require(dplyr)
# require(lme4)
# require(lmerTest)
# nobias <- read.csv("substitution_counts_gtr.csv")
# bias <- read.csv("substitution_counts_bias_gtr.csv")
# all %>% filter(bl >= 0.01) %>% mutate(ratio = ncount/scount) %>% na.omit() %>% filter(!is.infinite(ratio))->all.ratio
# summary(lmer(ratio ~ type + (1|ntaxa:bl) + (1|rep), data = all.ratio))
# # Fixed effects:
# #               Estimate Std. Error        df t value Pr(>|t|)    
# # (Intercept)  1.830e+00  1.315e-01 2.400e+01  13.916 5.60e-13 ***
# # typebias_gtr 4.720e-02  1.169e-02 2.026e+05   4.038 5.39e-05 ***
# # Bias simulations have a slightly higher ratio, meaning while both simulation sets show that nonsynonymous changes occur more frequently than do synonymous changes, t
