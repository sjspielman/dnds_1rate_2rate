# SJS
# This script builds mixed effects linear models to determine relative accuracy of models/methods


require(dplyr)
require(lme4)
require(multcomp)
require(readr)


# x is glht object
extract_multcomp <- function(x, type, model)
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
    dat$model <- model
    dat
}
      

dat.sum <- read_csv("dnds_summary.csv")

nobias <- dat.sum %>% na.omit() %>% filter(type == "nobias") # %>% filter(bl >= 0.01, ntaxa <= 1024)
nobias$method <- factor(nobias$method)
bias <- dat.sum %>% na.omit() %>% filter(type == "bias")    # %>% filter(bl >= 0.01, ntaxa <= 1024)
bias$method <- factor(bias$method)

fits <- data.frame("comp"    = character(), 
                  "coeff"   = numeric(),
                  "lowerCI" = numeric(),
                  "upperCI" = numeric(),
                  "pvalue"  = numeric(),
                  "sig"     = character(),
                  "model"   = character())

##################### Correlation linear models #######################
# NOTE: BL = 0.0025 excluded since results are largely meaningless at this low of a divergence level

nobias %>% filter(bl >= 0.01) -> nobias.sub
bias %>% filter(bl >= 0.01) -> bias.sub


fit.nobias <- lmer(r ~ method + (1|ntaxa:bl) + (1|rep), data = nobias.sub)
fit.nobias.mc <- glht(fit.nobias, linfct=mcp(method='Tukey'))
fits <- rbind(fits, extract_multcomp(fit.nobias.mc, "nobias", "r"))

fit.bias1 <- lmer(r ~ method + (1|ntaxa:bl) + (1|rep), data = bias.sub)
fit.bias1.mc <- glht(fit.bias1, linfct=mcp(method='Tukey'))
fits <- rbind(fits, extract_multcomp(fit.bias1.mc, "bias", "r"))




############### RMSD linear models #################
# NOTE: BL = <0.0025,0.01> and N = 128 excluded since results are largely meaningless at this low of a divergence level

nobias %>% filter(bl >= 0.04, ntaxa >= 256, rmsd <= 100 ) -> nobias.sub
bias %>% filter(bl >= 0.04, ntaxa >= 256, rmsd <= 100) -> bias.sub


fit.nobias <- lmer(rmsd ~ method + (1|ntaxa:bl) + (1|rep), data = nobias.sub)
fit.nobias.mc <- glht(fit.nobias, linfct=mcp(method='Tukey'))
fits <- rbind(fits, extract_multcomp(fit.nobias.mc, "nobias", "rmsd"))

fit.bias1 <- lmer(rmsd ~ method + (1|ntaxa:bl) + (1|rep), data = bias.sub)
fit.bias1.mc <- glht(fit.bias1, linfct=mcp(method='Tukey'))
fits <- rbind(fits, extract_multcomp(fit.bias1.mc, "bias", "rmsd"))


write.csv(fits, "linear_model_results.csv", quote = F, row.names = F)

