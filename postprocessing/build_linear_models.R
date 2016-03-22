# SJS
# This script builds mixed effects linear models to determine relative accuracy of models/methods
# NOTE: Condition BL=0.0025 is excluded from linear models since most results are meaningless at that low of a divergence level.


require(dplyr)
require(lme4)
require(multcomp)
require(readr)


# x is glht object
extract_multcomp <- function(x, type)
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
    dat$type <- type
    dat
}
      

dat.sum <- read_csv("dnds_summary.csv")

nobias <- dat.sum %>% na.omit() %>% filter(type == "nobias") # %>% filter(bl >= 0.01, ntaxa <= 1024)
nobias$method <- factor(nobias$method)
bias <- dat.sum %>% na.omit() %>% filter(type == "bias")    # %>% filter(bl >= 0.01, ntaxa <= 1024)
bias$method <- factor(bias$method)



##################### Correlation linear models #######################
# NOTE: BL = 0.0025 excluded since results are largely meaningless at this low of a divergence level

nobias %>% filter(truetype == "true1", bl >= 0.01) -> nobias.true1
bias %>% filter(truetype == "true1", bl >= 0.01) -> bias.true1
bias %>% filter(truetype == "true2", bl >= 0.01) -> bias.true2

r.fits <- data.frame("comp"    = character(), 
                  "coeff"   = numeric(),
                  "lowerCI" = numeric(),
                  "upperCI" = numeric(),
                  "pvalue"  = numeric(),
                  "sig"     = character())


fit.nobias <- lmer(r ~ method + (1|ntaxa:bl) + (1|rep), data = nobias.true1)
fit.nobias.mc <- glht(fit.nobias, linfct=mcp(method='Tukey'))
r.fits <- rbind(r.fits, extract_multcomp(fit.nobias.mc, "nobias"))

fit.bias1 <- lmer(r ~ method + (1|ntaxa:bl) + (1|rep), data = bias.true1)
fit.bias1.mc <- glht(fit.bias1, linfct=mcp(method='Tukey'))
r.fits <- rbind(r.fits, extract_multcomp(fit.bias1.mc, "bias1"))

fit.bias2 <- lmer(r ~ method + (1|ntaxa:bl) + (1|rep), data = bias.true2)
fit.bias2.mc <- glht(fit.bias2, linfct=mcp(method='Tukey'))
r.fits <- rbind(r.fits, extract_multcomp(fit.bias2.mc, "bias2"))

write.csv(r.fits, "linear_model_results_correlation.csv", quote = F, row.names = F)



############### RMSD linear models #################
# NOTE: BL = <0.0025,0.01,0.04> and N = <128,256> excluded since results are largely meaningless at this low of a divergence level

nobias %>% filter(truetype == "true1", bl >= 0.16, ntaxa >= 512) -> nobias.true1
bias %>% filter(truetype == "true1", bl >= 0.16, ntaxa >= 512) -> bias.true1
bias %>% filter(truetype == "true2", bl >= 0.16, ntaxa >= 512) -> bias.true2



rmsd.fits <- data.frame("comp"    = character(), 
                  "coeff"   = numeric(),
                  "lowerCI" = numeric(),
                  "upperCI" = numeric(),
                  "pvalue"  = numeric(),
                  "sig"     = character())


fit.nobias <- lmer(rmsd ~ method + (1|ntaxa:bl) + (1|rep), data = nobias.true1)
fit.nobias.mc <- glht(fit.nobias, linfct=mcp(method='Tukey'))
rmsd.fits <- rbind(rmsd.fits, extract_multcomp(fit.nobias.mc, "nobias"))

fit.bias1 <- lmer(rmsd ~ method + (1|ntaxa:bl) + (1|rep), data = bias.true1)
fit.bias1.mc <- glht(fit.bias1, linfct=mcp(method='Tukey'))
rmsd.fits <- rbind(rmsd.fits, extract_multcomp(fit.bias1.mc, "bias1"))

fit.bias2 <- lmer(rmsd ~ method + (1|ntaxa:bl) + (1|rep), data = bias.true2)
fit.bias2.mc <- glht(fit.bias2, linfct=mcp(method='Tukey'))
rmsd.fits <- rbind(rmsd.fits, extract_multcomp(fit.bias2.mc, "bias2"))

write.csv(rmsd.fits, "linear_model_results_rmsd.csv", quote = F, row.names = F)

