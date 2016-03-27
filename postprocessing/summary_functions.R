# SJS
# Functions used to summarize results, specifically correlations, RMSD, and estimator bias between inferred and true dN/dS, dN, and dS.

summarize_dnds <- function(dat, type)
{
   
    dat %>% 
      na.omit() %>% filter(!is.infinite(dnds)) %>% 
      group_by(ntaxa, bl, method, rep) %>% 
      do( bias.raw = glm(dnds ~ offset(true), dat = .),  
          cor.raw  = cor(.$true, .$dnds),
          rmsd.raw = sqrt(mean((.$true - .$dnds)^2))) %>% 
      mutate(estbias = summary(bias.raw)$coeff[1], r = cor.raw[1], rmsd = rmsd.raw[[1]]) %>% 
      select(-bias.raw, -cor.raw, -rmsd.raw) %>% 
      na.omit() %>% mutate(type = type)-> summ
  
    summ
    
}

summarize_dn <- function(dat, type)
{

    dat %>% select(-true, -trueds, -ds) %>% 
      na.omit() %>% filter(!is.infinite(dn)) %>% 
      group_by(ntaxa, bl, method, rep) %>% 
      do( bias.raw = glm(dn ~ offset(truedn), dat = .),  
          cor.raw  = cor(.$dn, .$truedn),
          rmsd.raw = sqrt(mean((.$dn - .$truedn)^2))) %>% 
      mutate(estbias = summary(bias.raw)$coeff[1], r = cor.raw[1], rmsd = rmsd.raw[[1]]) %>% 
      select(-bias.raw, -cor.raw, -rmsd.raw) %>% 
      na.omit() %>% mutate(type = type, parameter = "dn") -> summ

    summ
    
}

summarize_ds <- function(dat, type)
{

  dat %>% select(-true, -truedn, -dn) %>% 
    na.omit() %>% filter(!is.infinite(ds)) %>% 
    group_by(ntaxa, bl, method, rep) %>% 
    do( bias.raw = glm(ds ~ offset(trueds), dat = .),  
        cor.raw  = cor(.$ds, .$trueds),
        rmsd.raw = sqrt(mean((.$ds - .$trueds)^2))) %>% 
    mutate(estbias = summary(bias.raw)$coeff[1], r = cor.raw[1], rmsd = rmsd.raw[[1]]) %>% 
    select(-bias.raw, -cor.raw, -rmsd.raw) %>% 
    na.omit() %>% mutate(type = type, parameter = "ds") -> summ
  
  summ

  }


# Function to create summary data frame with information about mean correlation and estimator bias
summarize_dnds_real <- function(dat)
{

    # Summarize, for true2
    dat  %>% 
      na.omit() %>% filter(!is.infinite(dnds)) %>% 
      group_by(dataset, method, type, rep) %>% 
      do( bias.raw = glm(dnds ~ offset(true), dat = .),  
          cor.raw  = cor(.$true, .$dnds),
          rmsd.raw = sqrt(mean((.$true - .$dnds)^2))) %>% 
      mutate(estbias = summary(bias.raw)$coeff[1], r = cor.raw[1], rmsd = rmsd.raw[[1]]) %>% 
      select(-bias.raw, -cor.raw, -rmsd.raw) %>% na.omit() -> summ
    
    summ
}


summarize_dn_real <- function(dat)
{

    dat %>% select(-true, -trueds -ds) %>% 
      na.omit() %>% filter(!is.infinite(dn)) %>% group_by(dataset, type, method, rep) %>% 
      do( bias.raw = glm(dn ~ offset(truedn), dat = .),  
          cor.raw  = cor(.$dn, .$truedn),
          rmsd.raw = sqrt(mean((.$dn - .$truedn)^2))) %>% 
      mutate(estbias = summary(bias.raw)$coeff[1], r = cor.raw[1], rmsd = rmsd.raw[[1]]) %>% 
      select(-bias.raw, -cor.raw, -rmsd.raw) %>% 
      na.omit() %>% mutate(parameter = "dn") -> summ

    summ
    
}

summarize_ds_real <- function(dat)
{

    dat %>% select(-true, -truedn, -dn) %>% 
      na.omit() %>% filter(!is.infinite(ds)) %>% group_by(dataset, type, method, rep) %>% 
      do( bias.raw = glm(ds ~ offset(trueds), dat = .),  
          cor.raw  = cor(.$ds, .$trueds),
          rmsd.raw = sqrt(mean((.$ds - .$trueds)^2))) %>% 
      mutate(estbias = summary(bias.raw)$coeff[1], r = cor.raw[1], rmsd = rmsd.raw[[1]]) %>% 
      select(-bias.raw, -cor.raw, -rmsd.raw) %>% 
      na.omit() %>% mutate(parameter = "ds") -> summ
    
    summ
}