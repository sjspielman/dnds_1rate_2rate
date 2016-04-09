# SJS
# This script builds mixed effects linear models to determine relative accuracy of models/methods. 
# Note that model fits *exclude* the BL=0.0025 condition, as this condition is mostly noise.

require(dplyr)
require(lme4)
require(multcomp)
require(readr)


# x is glht object
extract_multcomp <- function(x, model, type, pi)
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
    dat$model <- model
    dat$type  <- type
    dat$pi    <- pi
    dat
}



fits <- data.frame("comp"    = character(),
                   "coeff"   = numeric(),
                   "lowerCI" = numeric(),
                   "upperCI" = numeric(),
                   "pvalue"  = numeric(),
                   "sig"     = character(),
                   "model"   = character(),
                   "pi"      = character())
      
for (pi in c("unequalpi", "equalpi"))
{
    
    for (type in c("nobias", "bias"))
    {

        dat.sum <- read_csv(paste0("dnds_summary_", pi, "_", type, ".csv"))
        dat.sum %>% filter(bl >= 0.01) %>% na.omit() %>% filter(!is.infinite(rmsd)) -> dat.sum
        dat.sum$method <- factor(dat.sum$method)

    
        fit <- lmer(r ~ method + (1|ntaxa:bl) + (1|rep), data = dat.sum)
        fit.mc <- glht(fit, linfct=mcp(method='Tukey'))
        fits <- rbind(fits, extract_multcomp(fit.mc, "r", type, pi))
        
        fit <- lmer(rmsd ~ method + (1|ntaxa:bl) + (1|rep), data = dat.sum)
        fit.mc <- glht(fit, linfct=mcp(method='Tukey'))
        fits <- rbind(fits, extract_multcomp(fit.mc, "rmsd", type, pi))

        fit <- lmer(resvar ~ method + (1|ntaxa:bl) + (1|rep), data = dat.sum)
        fit.mc <- glht(fit, linfct=mcp(method='Tukey'))
        fits <- rbind(fits, extract_multcomp(fit.mc, "resvar", type, pi))

  
    }   
} 
    write.csv(fits, "linear_model_results.csv", quote = F, row.names = F)
