# SJS
# Build and save summary statistic dataframes for dN/dS inference. For each parameter dN/dS, dN, and dS, this script computes the Pearson correlation, estimator bias, and RMSD for all inferences. Note that RMSD and estimator bias for SLAC dN and dS parameter estimates are meaningless due to scaling, but they are processed anyways.

require(dplyr)
require(tidyr)
require(readr)

# Function to create summary data frame with information about mean correlation and estimator bias
summarize_dnds <- function(dat, type)
{

    # Summarize, for true1
    dat %>% select(-true2) %>% na.omit() %>% filter(!is.infinite(dnds)) %>% group_by(ntaxa, bl, method, rep) %>% 
    do( bias.raw = glm(dnds ~ offset(true1), dat = .),  
        cor.raw  = cor(.$true1, .$dnds),
        rmsd.raw = sqrt(mean((.$true1 - .$dnds)^2))) %>% 
    mutate(estbias = summary(bias.raw)$coeff[1], r = cor.raw[1], rmsd = rmsd.raw[[1]]) %>% 
    select(-bias.raw, -cor.raw, -rmsd.raw) %>% na.omit() %>% mutate(truetype = "true1") -> summ1
    
    # Summarize, for true2
    dat %>% select(-true1) %>% na.omit() %>% filter(!is.infinite(dnds)) %>% group_by(ntaxa, bl, method, rep) %>% 
    do( bias.raw = glm(dnds ~ offset(true2), dat = .),  
        cor.raw  = cor(.$true2, .$dnds),
        rmsd.raw = sqrt(mean((.$true2 - .$dnds)^2))) %>% 
    mutate(estbias = summary(bias.raw)$coeff[1], r = cor.raw[1], rmsd = rmsd.raw[[1]]) %>% 
    select(-bias.raw, -cor.raw, -rmsd.raw) %>% na.omit() %>% mutate(truetype = "true2") -> summ2
    
    summ.full <- rbind(summ1, summ2)
    summ.full$type <- type
    summ.full
    
}

summarize_dn <- function(dat, type)
{

    dat %>% select(-true1, -true2, -ds) %>% na.omit() %>% filter(!is.infinite(dn)) %>% group_by(ntaxa, bl, method, rep) %>% 
    do( bias.raw = glm(dn ~ offset(truedn), dat = .),  
        cor.raw  = cor(.$dn, .$truedn),
        rmsd.raw = sqrt(mean((.$dn - .$truedn)^2))) %>% 
    mutate(estbias = summary(bias.raw)$coeff[1], r = cor.raw[1], rmsd = rmsd.raw[[1]]) %>% 
    select(-bias.raw, -cor.raw, -rmsd.raw) %>% na.omit() -> summ
    
    summ$type <- type
    summ$parameter <- "dn"
    summ
    
}

summarize_ds <- function(dat, type)
{

    dat %>% select(-true1, -true2, -dn) %>% na.omit() %>% filter(!is.infinite(ds)) %>% group_by(ntaxa, bl, method, rep) %>% 
    do( bias.raw = glm(ds ~ offset(trueds), dat = .),  
        cor.raw  = cor(.$ds, .$trueds),
        rmsd.raw = sqrt(mean((.$ds - .$trueds)^2))) %>% 
    mutate(estbias = summary(bias.raw)$coeff[1], r = cor.raw[1], rmsd = rmsd.raw[[1]]) %>% 
    select(-bias.raw, -cor.raw, -rmsd.raw) %>% na.omit() -> summ
    
    summ$type <- type
    summ$parameter <- "ds"
    summ
}


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

