dnds.sum.mean %>% na.omit() %>%
  filter(pi == "unequal") %>%
  select(-pi, -rmsd, -estbias, -rmsd, -resvar) %>%
  group_by(ntaxa, bl) %>%
  spread(method,r) %>%
  mutate(FEL   = FEL1 - FEL2,
         FUBAR = FUBAR1 - FUBAR2,
         SLAC  = SLAC1 - SLAC2) %>%
  select(-FUBAR1, -FUBAR2, -FEL1, -FEL2, -SLAC1, -SLAC2) %>%
  gather(method, rdiff, FEL, FUBAR, SLAC) %>%
  ggplot(aes(x = as.factor(ntaxa), y = as.factor(bl), fill = rdiff)) +
    geom_tile() + geom_text(aes(label = round(rdiff,3)), fontface="bold") +
    scale_fill_gradient2(name = "Correlation Difference", limits=c(-0.2, 0.2)) + facet_grid(method~type) +
    theme(legend.position = "bottom", legend.box = "horizontal", legend.key.size = unit(0.6, "cm"), legend.key.width = unit(.8, "cm"))

    # guides(fill = guide_legend(title = "Correlation difference", title.position = "top"))



for (pi in c("unequal", "equal"))
{
  for (type in c("nobias", "bias"))
  {
    for (n in c(128, 256, 512, 1024, 2048))
    {
      for (bl in c(0.0025, 0.01, 0.04, 0.16, 0.64))
      {
          dnds.sum %>% filter(pi == pi, type == type, ntaxa == n, bl == bl) %>% na.omit() -> testme
          fit <- lmer(r ~ method + (1|rep), data = testme)
          fit.mc <- glht(fit, linfct=mcp(method='Tukey'))
          fits <- rbind(fits, extract_multcomp(fit.mc, "r", type, pi))


          dnds.sum %>% group_by(pi, type, ntaxa, bl) %>%
            do(fit = lmer(r ~ method + (1|rep), data=.)) %>%
            do(mc = glht(.$fit, linfct=mcp(method='Tukey'))) %>%
            mutate(p = summary(confint(mc))$test$pvalues)

      }
    }
  }
}

        dat.sum <- read_csv(paste0("dnds_summary_", pi, "pi_", type, ".csv"))

        dat.sum$method <- factor(dat.sum$method)


        fit <- lmer(r ~ method + (1|ntaxa:bl) + (1|rep), data = dat.sum)
        fit.mc <- glht(fit, linfct=mcp(method='Tukey'))
        fits <- rbind(fits, extract_multcomp(fit.mc, "r", type, pi))



# x is glht object
mini.extract.multcomp <- function(x, model, type)
{



    confsum <- confint(summary(x))$confint
    pvalue <- summary(confint(x))$test$pvalues[i]
        if (pvalue < 0.01){
            sig <- TRUE
        }
        else{
            sig <- FALSE
        }
        if (swap){
        temp <- data.frame("comp" = comp, "coeff" = coeff, "lowerCI" = lowerCI, "upperCI" = upperCI, "pvalue" = pvalue, "sig" = sig)
        dat <- rbind(dat, temp)
    }
    dat$model <- model
    dat$type  <- type
    dat$pi    <- pi
    dat
}


dnds.sum %>% filter(pi == "unequal") %>%
  select(ntaxa, bl, method, rep, r, type) %>%
  group_by(ntaxa, bl, method, type) %>%
