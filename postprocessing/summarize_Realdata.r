
# Function to create summary data frame with information about mean correlation and estimator bias
summarize_dnds <- function(dat)
{

    # Summarize, for true1
    dat %>% select(-true2) %>% na.omit() %>% filter(!is.infinite(dnds)) %>% group_by(dataset, method, type, rep) %>% 
    do( bias.raw = glm(dnds ~ offset(true1), dat = .),  
        cor.raw  = cor(.$true1, .$dnds),
        rmsd.raw = sqrt(mean((.$true1 - .$dnds)^2))) %>% 
    mutate(estbias = summary(bias.raw)$coeff[1], r = cor.raw[1], rmsd = rmsd.raw[[1]]) %>% 
    select(-bias.raw, -cor.raw, -rmsd.raw) %>% na.omit() %>% mutate(truetype = "true1") -> summ1
    
    # Summarize, for true2
    dat %>% select(-true1) %>% na.omit() %>% filter(!is.infinite(dnds)) %>% group_by(dataset, method, type, rep) %>% 
    do( bias.raw = glm(dnds ~ offset(true2), dat = .),  
        cor.raw  = cor(.$true2, .$dnds),
        rmsd.raw = sqrt(mean((.$true2 - .$dnds)^2))) %>% 
    mutate(estbias = summary(bias.raw)$coeff[1], r = cor.raw[1], rmsd = rmsd.raw[[1]]) %>% 
    select(-bias.raw, -cor.raw, -rmsd.raw) %>% na.omit() %>% mutate(truetype = "true2") -> summ2
    
    summ.full <- rbind(summ1, summ2)
    summ.full
    
}


summarize_dn <- function(dat)
{

    dat %>% select(-true1, -true2, -ds) %>% na.omit() %>% filter(!is.infinite(dn)) %>% group_by(dataset, type, method, rep) %>% 
    do( bias.raw = glm(dn ~ offset(truedn), dat = .),  
        cor.raw  = cor(.$dn, .$truedn),
        rmsd.raw = sqrt(mean((.$dn - .$truedn)^2))) %>% 
    mutate(estbias = summary(bias.raw)$coeff[1], r = cor.raw[1], rmsd = rmsd.raw[[1]]) %>% 
    select(-bias.raw, -cor.raw, -rmsd.raw) %>% na.omit() -> summ
    
    summ$parameter <- "dn"
    summ
    
}

summarize_ds <- function(dat)
{

    dat %>% select(-true1, -true2, -dn) %>% na.omit() %>% filter(!is.infinite(ds)) %>% group_by(dataset, type, method, rep) %>% 
    do( bias.raw = glm(ds ~ offset(trueds), dat = .),  
        cor.raw  = cor(.$ds, .$trueds),
        rmsd.raw = sqrt(mean((.$ds - .$trueds)^2))) %>% 
    mutate(estbias = summary(bias.raw)$coeff[1], r = cor.raw[1], rmsd = rmsd.raw[[1]]) %>% 
    select(-bias.raw, -cor.raw, -rmsd.raw) %>% na.omit() -> summ
    
    summ$parameter <- "ds"
    summ
}

require(cowplot)


realdat %>% filter(method %in% c("SLAC2", "FUBAR2"), type == "bias_gtr") %>% summarize_dn() -> dn
realdat %>% filter(method %in% c("SLAC2", "FUBAR2"), type == "bias_gtr") %>% summarize_ds() -> ds
rbind(dn,ds) %>% filter(method == "FUBAR2") %>% ggplot(aes(x = dataset, fill = parameter, y = rmsd)) + geom_boxplot()

realdat %>% filter(method == "FUBAR2", type == "gtr") %>% ggplot(aes(x = dataset, y = ds)) + 
                                                          geom_boxplot(fill = "cadetblue") + 
                                                          geom_hline(yintercept=1)



realdat.sum <- realdat %>% summarize_dnds() 
realdat.sum$method <- factor(realdat.sum$method, levels = c("SLAC1", "FUBAR1", "SLAC2", "FUBAR2"))
realdat.sum$dataset <- factor(realdat.sum$dataset, levels = c("amine", "camelid", "vertrho", "h3", "hivrt"))
fill_colors <- c("red", "orange", "cadetblue", "lightblue")
r    <- realdat.sum %>% filter(truetype == "true2") %>% ggplot(aes(x = dataset, fill = method, y = r)) + geom_boxplot() + facet_wrap(~type) + scale_fill_manual(values=fill_colors)
rmsd <- realdat.sum %>% filter(truetype == "true2") %>% ggplot(aes(x = dataset, fill = method, y = rmsd)) + geom_boxplot() + facet_wrap(~type) + scale_fill_manual(values=fill_colors)
plot_grid(r,rmsd, nrow=2)



