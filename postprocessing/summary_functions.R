# SJS
# Functions used to summarize results, specifically correlations, RMSD, and estimator bias between inferred and true dN/dS, dN, and dS.


summarize_means <- function(dat)
{
  dat %>%
  na.omit() %>% filter(!is.infinite(rmsd), !is.infinite(resvar)) %>%
  summarize(r = mean(r), estbias = mean(estbias), rmsd = mean(rmsd)) -> dat2
  dat2
}

summarize_dnds_empirical <- function(dat)
{
  dat %>%
  na.omit() %>% filter(!is.infinite(dnds), !is.infinite(empdnds)) %>%
  do( bias.raw = glm(dnds ~ offset(empdnds), dat = .),
      cor.raw = cor(.$dnds, .$empdnds),
      lm.raw   = lm(dnds ~ empdnds, dat = .),
      rmsd.raw = sqrt(mean((.$empdnds - .$dnds)^2))) %>%
  mutate(estbias = summary(bias.raw)$coeff[1],
         r = cor.raw[1],
         rmsd = rmsd.raw[[1]],
         resvar = summary(lm.raw)$sigma^2) %>%
  select(-bias.raw, -rmsd.raw, -lm.raw, -cor.raw) %>%
  na.omit() -> sumdat
  sumdat$rmsd[sumdat$rmsd >= 1000] <- Inf
  sumdat$resvar[sumdat$resvar >= 1000] <- Inf
  sumdat
}

summarize_dnds <- function(dat)
{

    dat %>%
      na.omit() %>% filter(!is.infinite(dnds)) %>%
      group_by(ntaxa, bl, method, rep, biastype, pitype) %>%
      do( bias.raw = glm(dnds ~ offset(truednds), dat = .),
          cor.raw = cor(.$dnds, .$truednds),
          lm.raw   = lm(dnds ~ truednds, dat = .),
          rmsd.raw = sqrt(mean((.$truednds - .$dnds)^2))) %>%
      mutate(estbias = summary(bias.raw)$coeff[1], r = cor.raw[1], rmsd = rmsd.raw[[1]], resvar = summary(lm.raw)$sigma^2) %>%
      select(-bias.raw, -rmsd.raw, -lm.raw, -cor.raw) %>%
      na.omit() -> summ
    summ$rmsd[summ$rmsd >= 1000] <- Inf
    summ$resvar[summ$resvar >= 1000] <- Inf

    summ
}

summarize_dn  <- function(dat)
{

    dat %>% select(-truednds, -trueds, -ds) %>%
      na.omit() %>% filter(!is.infinite(dn)) %>%
      group_by(ntaxa, bl, method, rep, biastype, pitype) %>%
      do( bias.raw = glm(dn ~ offset(truedn), dat = .),
          cor.raw = cor(.$dn, .$truedn),
          lm.raw   = lm(dn ~ truedn, dat = .),
          rmsd.raw = sqrt(mean((.$dn - .$truedn)^2))) %>%
      mutate(estbias = summary(bias.raw)$coeff[1], r = cor.raw[1], rmsd = rmsd.raw[[1]], resvar = summary(lm.raw)$sigma^2) %>%
      select(-bias.raw, -rmsd.raw, -lm.raw, -cor.raw) %>%
      na.omit() %>% mutate(parameter = "dn") -> summ
    summ$rmsd[summ$rmsd >= 1000] <- Inf
    summ$resvar[summ$resvar >= 1000] <- Inf

    summ

}

summarize_ds <- function(dat)
{

  dat %>% select(-truednds, -truedn, -dn) %>%
    na.omit() %>% filter(!is.infinite(ds)) %>%
    group_by(ntaxa, bl, method, rep, biastype, pitype) %>%
    do( bias.raw = glm(ds ~ offset(trueds), dat = .),
        cor.raw = cor(.$ds, .$trueds),
        lm.raw   = lm(ds ~ trueds, dat = .),
        rmsd.raw = sqrt(mean((.$ds - .$trueds)^2))) %>%
      mutate(estbias = summary(bias.raw)$coeff[1], r = cor.raw[1], rmsd = rmsd.raw[[1]], resvar = summary(lm.raw)$sigma^2) %>%
      select(-bias.raw, -rmsd.raw, -lm.raw, -cor.raw) %>%
    na.omit() %>% mutate(parameter = "ds") -> summ
    summ$rmsd[summ$rmsd >= 1000] <- Inf
    summ$resvar[summ$resvar >= 1000] <- Inf

  summ

}


# Function to create summary data frame with information about mean correlation and estimator bias
summarize_dnds_real <- function(dat)
{

    dat  %>%
      na.omit() %>% filter(!is.infinite(dnds)) %>%
      group_by(dataset, method, biastype, rep) %>%
      do( bias.raw = glm(dnds ~ offset(truednds), dat = .),
          cor.raw = cor(.$dnds, .$truednds),
          lm.raw   = lm(dnds ~ truednds, dat = .),
          rmsd.raw = sqrt(mean((.$truednds - .$dnds)^2))) %>%
      mutate(estbias = summary(bias.raw)$coeff[1], r = cor.raw[1], rmsd = rmsd.raw[[1]], resvar = summary(lm.raw)$sigma^2) %>%
      select(-bias.raw, -rmsd.raw, -lm.raw, -cor.raw) %>%
      na.omit() -> summ
    summ$rmsd[summ$rmsd >= 1000] <- Inf
    summ$resvar[summ$resvar >= 1000] <- Inf

    summ
}


summarize_dn_real <- function(dat)
{

    dat %>% select(-truednds, -trueds -ds) %>%
      na.omit() %>% filter(!is.infinite(dn)) %>% group_by(dataset, biastype, method, rep) %>%
      do( bias.raw = glm(dn ~ offset(truedn), dat = .),
          cor.raw = cor(.$dn, .$truedn),
          lm.raw   = lm(dn ~ truedn, dat = .),
          rmsd.raw = sqrt(mean((.$dn - .$truedn)^2))) %>%
      mutate(estbias = summary(bias.raw)$coeff[1], r = cor.raw[1], rmsd = rmsd.raw[[1]], resvar = summary(lm.raw)$sigma^2) %>%
      select(-bias.raw, -rmsd.raw, -lm.raw, -cor.raw) %>%
      na.omit() %>% mutate(parameter = "dn") -> summ
    summ$rmsd[summ$rmsd >= 1000] <- Inf
    summ$resvar[summ$resvar >= 1000] <- Inf

    summ

}

summarize_ds_real <- function(dat)
{

    dat %>% select(-truednds, -truedn, -dn) %>%
      na.omit() %>% filter(!is.infinite(ds)) %>% group_by(dataset, biastype, method, rep) %>%
      do( bias.raw = glm(ds ~ offset(trueds), dat = .),
          cor.raw = cor(.$ds, .$trueds),
          lm.raw   = lm(ds ~ trueds, dat = .),
          rmsd.raw = sqrt(mean((.$ds - .$trueds)^2))) %>%
      mutate(estbias = summary(bias.raw)$coeff[1], r = cor.raw[1], rmsd = rmsd.raw[[1]], resvar = summary(lm.raw)$sigma^2) %>%
      select(-bias.raw, -rmsd.raw, -lm.raw, -cor.raw) %>%
      na.omit() %>% mutate(parameter = "ds") -> summ
    summ$rmsd[summ$rmsd >= 1000] <- Inf
    summ$resvar[summ$resvar >= 1000] <- Inf

    summ
}
