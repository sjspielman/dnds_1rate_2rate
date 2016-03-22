# SJS
# Functions used to summarize results.

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


# Function to create summary data frame with information about mean correlation and estimator bias
summarize_dnds_real <- function(dat)
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


summarize_dn_real <- function(dat)
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

summarize_ds_real <- function(dat)
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