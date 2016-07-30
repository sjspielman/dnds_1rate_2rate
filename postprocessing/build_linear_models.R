# SJS
# This script builds mixed effects linear models to determine relative accuracy of models/methods.
# Note that model fits *exclude* the BL=0.0025 condition.

require(dplyr)
require(lme4)
require(multcomp)
require(readr)
require(lmerTest) # summary() function with p-value, for heatmap linear models
OUTDIR="dataframes/"

# x is glht object
extract_multcomp <- function(x, model, biastype, pitype, alpha)
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
        if (pvalue <= alpha){
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
    dat$biastype <- biastype
    dat$pitype   <- pitype
    dat
}


#### Linear models for comparison with true (stationary) dN/dS ####

dat <- read_csv("dataframes/summary_balanced_dnds.csv") %>% na.omit() %>% filter(bl >= 0.01)
dat$method <- factor(dat$method)
alpha <- 0.01
fits <- data.frame("comp"     = character(),
                   "coeff"    = numeric(),
                   "lowerCI"  = numeric(),
                   "upperCI"  = numeric(),
                   "pvalue"   = numeric(),
                   "sig"      = character(),
                   "model"    = character(),
                   "pitype"   = character(),
                   "biastype" = character())

for (pitype in c("unequalpi", "equalpi"))
{

    for (biastype in c("nobias", "bias"))
    {

        subdat <- dat %>% filter(pitype == pitype, biastype == biastype)

        fit <- lmer(r ~ method + (1|ntaxa:bl) + (1|rep), data = subdat)
        fit.mc <- glht(fit, linfct=mcp(method='Tukey'))
        fits <- rbind(fits, extract_multcomp(fit.mc, "r", biastype, pitype, alpha))

        subdat <- subdat %>% filter(!is.infinite(rmsd))

        fit <- lmer(rmsd ~ method + (1|ntaxa:bl) + (1|rep), data = subdat)
        fit.mc <- glht(fit, linfct=mcp(method='Tukey'))
        fits <- rbind(fits, extract_multcomp(fit.mc, "rmsd", biastype, pitype, alpha))
    }
}
write_csv(fits, paste0(OUTDIR, "linear_models.csv"))




################################################################################
################# Linear models for heat map figures ###########################

alpha <- 0.05/300 # Correct for 100 models for 3 methods


build.heatmap.models <- function(data, methods, title, alpha)
{
  data %>%
    filter(method %in% methods) %>%
    dplyr::select(-rmsd, -resvar, -estbias) %>%
    group_by(ntaxa, bl, biastype, pitype) %>%
    spread(method, r) %>%
    setNames(c("ntaxa", "bl", "replicate", "biastype", "pitype", "onerate", "tworate")) %>%
    na.omit() %>%
    do(tt = t.test(.$onerate, .$tworate, paired=T)) %>%
    tidy(tt) %>%
    mutate(r.diff = estimate,
           sig = as.numeric(p.value <= alpha)) %>%
    dplyr::select(ntaxa, bl, biastype, pitype, r.diff, sig) %>%
    mutate(method = title) -> fitted.dat

    fitted.dat
}

dat <- read_csv("dataframes/summary_balanced_dnds.csv") %>% na.omit() %>% filter(!is.infinite(rmsd))

# Filter uninformative conditions from FEL to avoid lm errors
fel.dat <- dat %>% filter(bl >= 0.01, !(bl == 0.01 & ntaxa <= 512))
# NA rows for FEL
fel.na.rows <-    data.frame(ntaxa = c( rep(128,4),rep(256,4), rep(512,4), rep(1024, 4), rep(2048,4)),
                             bl = 0.0025,
                             biastype = rep(c("nobias", "bias"), 10),
                             pitype = rep(c("equalpi", "unequalpi", "unequalpi", "equalpi"), 5)) %>%
            rbind(data.frame(ntaxa = c( rep(128,4),rep(256,4), rep(512,4)),
                             bl = 0.01,
                             biastype = rep(c("nobias", "bias"), 6),
                             pitype = rep(c("equalpi", "unequalpi", "unequalpi", "equalpi"), 3))) %>%
            mutate(r.diff = NA, sig = NA, method = "FEL")

fel.heatmap.data.true <- build.heatmap.models(fel.dat, c("FEL1", "FEL2"), "FEL", alpha) %>% rbind(fel.na.rows) %>% arrange(ntaxa, bl, biastype, pitype)
heatmap.data <-  rbind(fel.heatmap.data.true, build.heatmap.models(dat, c("SLAC1", "SLAC2"), "SLAC", alpha)) %>%
                      rbind(build.heatmap.models(dat, c("FUBAR1", "FUBAR2"), "FUBAR", alpha))%>%
                      ungroup()
write_csv(heatmap.data, paste0(OUTDIR, "linmodels_heatmap_version.csv"))
#
# # finding out which FEL comparisons are nonfunctional
# for (bl in c(0.0025, 0.01, 0.04)){
#   for (nt in c(128, 256, 512)){
#     for (bias in c("bias", "nobias")){
#       for (pi in c("unequalpi", "uequalpi")){
#         print(c(bl,nt,bias,pi))
#         if (bl == 0.0025 & ntaxa == 128 & pitype == "unequalpi")next
#         if (bl == 0.01 & ntaxa == 128 & pitype == "unequalpi")next
#         if (bl == 0.04 & ntaxa == 128 & pitype == "unequalpi")next
#         d2 <- fel.dat %>% filter(bl == bl, ntaxa == nt, method %in% c("FEL1", "FEL2")) %>% spread(method,r)
#         tt <- t.test(d2$FEL1, d2$FEL2, paired=T)
#         print(tt)
#       }
#     }
#   }
#
# }
#


########## Code dump! Difference between equal, unequal, bias, nobias simulations, overall? ###########
## TL;DR -  Equal and unequal are basically the same (all P>0.01). Bias has lower correlations than does nobias, by all measures.
#
# # Whether the simulation sets show differences?
# a <- read_csv("dnds_summary_equalpi_nobias.csv")
# a$pi <- "equal"
# b <- read_csv("dnds_summary_equalpi_bias.csv")
# b$pi <- "equal"
# d <- read_csv("dnds_summary_unequalpi_nobias.csv")
# d$pi <- "unequal"
# e <- read_csv("dnds_summary_unequalpi_bias.csv")
# e$pi <- "unequal"
# full <- rbind(a,b)
# full <- rbind(full,d)
# full <- rbind(full,e)
# # Note, the following was done for FEL1, FUBAR1, SLAC1. For SLAC1 only, the pi effects were significant, but all magnitudes were on order of 10^-3, so differences are exceedingly minimal.
# full.fel1 <- full %>% filter(method == "FEL1", bl >= 0.01) %>% na.omit() %>% filter(!is.infinite(rmsd), !is.infinite(resvar))
# fit <- lmer(r ~ type + pi + (1|ntaxa:bl) + (1|rep), data = full.fel1)
# summary(fit)
# #               Estimate Std. Error         df t value Pr(>|t|)
# # (Intercept)  6.645e-01  3.943e-02  1.900e+01  16.852 6.68e-13 ***
# # typenobias   5.395e-02  1.766e-03  3.929e+03  30.540  < 2e-16 ***
# # piunequal   -7.176e-04  1.766e-03  3.929e+03  -0.406    0.685
#
# fit <- lmer(rmsd ~ type + pi + (1|ntaxa:bl) + (1|rep), data = full.fel1)
# summary(fit)
# #               Estimate Std. Error         df t value Pr(>|t|)
# # (Intercept)  3.558e-01  4.079e-02  1.900e+01   8.723 4.43e-08 ***
# # typenobias  -6.771e-02  2.196e-03  3.978e+03 -30.830  < 2e-16 ***
# # piunequal   -4.455e-03  2.196e-03  3.978e+03  -2.029   0.0426 *
#
# fit <- lmer(resvar ~ type + pi + (1|ntaxa:bl) + (1|rep), data = full.fel1)
# summary(fit)
# #               Estimate Std. Error         df t value Pr(>|t|)
# # (Intercept)  1.648e-01  3.363e-02  1.900e+01   4.900 9.75e-05 ***
# # typenobias  -5.524e-02  2.913e-03  3.978e+03 -18.966  < 2e-16 ***
# # piunequal   -3.628e-03  2.913e-03  3.978e+03  -1.246    0.213
