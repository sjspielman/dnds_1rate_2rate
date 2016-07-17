# SJS
# This script builds mixed effects linear models to determine relative accuracy of models/methods.
# Note that model fits *exclude* the BL=0.0025 condition, as this condition is mostly noise.

require(dplyr)
require(lme4)
require(multcomp)
require(readr)
OUTDIR="dataframes/"

# x is glht object
extract_multcomp <- function(x, model, btype)
{
    confsum <- confint(summary(x))$confint
    testnames <- names(summary(x)$test$tstat)
    dat <- data.frame("comp"    = character(),
                      "coeff"   = numeric(),
                      "lowerCI" = numeric(),
                      "upperCI" = numeric(),
                      "pvalue"  = numeric(),
                      "sig"     = character())
    for (i in 1:length(testnames)){
        #print(i)
        comp <- testnames[i]
        swap <- TRUE

        # One-rate must come first
        backwards <- grep("[A-Z]+2 - [A-Z]+1", comp) # 1rate before 2rate
        if (length(backwards) == 1)
        {
            comp <- gsub("([A-Z]+2) - ([A-Z]+1)", "\\2 - \\1", comp)
            swap <- TRUE
        }

        # SLAC must come second
        bad_order1 <- grep("SLAC[12] - F[A-Z]+[12]", comp)
        if (length(bad_order1) == 1)
        {
            comp <- gsub("(SLAC[12]) - ([A-Z]+[12])", "\\2 - \\1", comp)
            swap <- TRUE
        }

        # FEL must come first
        bad_order2<- grep("[A-Z]+1 - FEL1", comp)
        if (length(bad_order2) == 1)
        {
            comp <- gsub("([A-Z]+[12]) - (FEL[12])", "\\2 - \\1", comp)
            swap <- TRUE
        }

        coeff <- confsum[i]
        pvalue <- summary(confint(x))$test$pvalues[i]
        if (pvalue < 0.01){
            sig <- TRUE
        }
        else{
            sig <- FALSE
        }
        if (swap){
            coeff <- coeff * -1
            lowerCI <- confsum[i,3] * -1
            upperCI <- confsum[i,2] * -1
        }
        else{
            upperCI <- confsum[i,3]
            lowerCI <- confsum[i,2]
        }

        temp <- data.frame("comp" = comp, "coeff" = coeff, "lowerCI" = lowerCI, "upperCI" = upperCI, "pvalue" = pvalue, "sig" = sig)
        dat <- rbind(dat, temp)
    }
    dat$model    <- model
    dat$biastype <- btype
    dat
}



fits <- data.frame("comp"    = character(),
                   "coeff"   = numeric(),
                   "lowerCI" = numeric(),
                   "upperCI" = numeric(),
                   "pvalue"  = numeric(),
                   "sig"     = character(),
                   "model"   = character(),
                   "pitype"  = character())

dat <- read.csv("summary_dnds_empirical.csv")
for (btype in c("nobias", "bias"))
{
    dat %>% filter(bl >= 0.01, biastype == btype) %>%
      na.omit() %>%
      filter(!is.infinite(rmsd), !is.infinite(resvar)) %>%
      filter(rmsd <= 1e5, resvar <= 1e5)-> dat.sum
    dat.sum$method <- factor(dat.sum$method)

    fit <- lmer(r ~ method + (1|ntaxa:bl) + (1|rep), data = dat.sum)
    fit.mc <- glht(fit, linfct=mcp(method='Tukey'))
    fits <- rbind(fits, extract_multcomp(fit.mc, "r", btype))

    fit <- lmer(rmsd ~ method + (1|ntaxa:bl) + (1|rep), data = dat.sum)
    fit.mc <- glht(fit, linfct=mcp(method='Tukey'))
    fits <- rbind(fits, extract_multcomp(fit.mc, "rmsd", btype))

    fit <- lmer(resvar ~ method + (1|ntaxa:bl) + (1|rep), data = dat.sum)
    fit.mc <- glht(fit, linfct=mcp(method='Tukey'))
    fits <- rbind(fits, extract_multcomp(fit.mc, "resvar", btype))
}


comp.order <- c("FUBAR1 - SLAC1", "FEL1 - SLAC1", "FEL1 - FUBAR1", "SLAC1 - SLAC2", "FUBAR1 - FUBAR2", "FEL1 - FEL2")
fits <- filter(fits, comp %in% comp.order)
fits$model <- factor(fits$model, levels = c("r", "rmsd", "resvar"))
fits$comp <- factor(fits$comp, levels = comp.order)
fits$biastype <- factor(fits$biastype, levels=c("nobias", "bias"), labels=c("No codon bias", "Codon bias"))


sig_colors <- c("grey60", "black")
ptsize = 1.5
linesize = 0.8

fits %>% filter(model == "r") %>% arrange(comp) %>%
  ggplot(aes(x = coeff, y = comp, color = sig)) + geom_point(size=ptsize) +
  geom_segment(aes(x=lowerCI,xend=upperCI,y=comp,yend=comp),size=linesize) +
  geom_vline(xintercept=0, size=0.5) +
  facet_grid(~biastype) + background_grid() +
  xlab("Average Correlation Difference") + ylab("Comparison") +
  scale_color_manual(values = sig_colors)  + theme(legend.position = "none", axis.text.y = element_text(size=10)) -> linmodel.r

fits %>% filter(model == "rmsd") %>% arrange(comp) %>%
  ggplot(aes(x = coeff, y = comp, color = sig)) + geom_point(size=ptsize) +
  geom_segment(aes(x=lowerCI,xend=upperCI,y=comp,yend=comp),size=linesize) +
  geom_vline(xintercept=0, size=0.5) +
  facet_grid(~biastype) + background_grid() +
  xlab("Average RMSD Difference") + ylab("Comparison") +
  scale_color_manual(values = sig_colors) + theme(legend.position = "none", axis.text.y = element_text(size=10)) -> linmodel.rmsd


fits %>% filter(model == "resvar") %>% arrange(comp) %>%
  ggplot(aes(x = coeff, y = comp, color = sig)) + geom_point(size=ptsize) +
  geom_segment(aes(x=lowerCI,xend=upperCI,y=comp,yend=comp),size=linesize) +
  geom_vline(xintercept=0, size=0.5) +
  facet_grid(~biastype) + background_grid() +
  xlab("Average Residual Variance Difference") + ylab("Comparison") +
  scale_color_manual(values = sig_colors) + theme(legend.position = "none", axis.text.y = element_text(size=10)) -> linmodel.resvar

x <- plot_grid(linmodel.r, linmodel.rmsd, nrow=2, labels=c("A", "B"))
print(x)
