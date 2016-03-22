# SJS
# Script to merge dN/dS inferences for each dataset into a single data frame.

require(dplyr)
require(readr)
require(cowplot)

max_threshold = 1e2 # if dN/dS greater than this, clearly we're in an uninformative situation. Also, 5e-9 is HyPhy's bailed 0.

# Clean dN/dS values for those methods which return separate dN and dS (hence divide)
# Methods include SLAC, FUBAR1, FUBAR2
clean_dnds_divide <- function(w, numcol)
{
    for (i in 1:numcol)
    {
        if (is.na(w[i]))
        {
            w[i] <- NA
        }
        else if (w[i] >= max_threshold | w[i] == Inf)
        {
            w[i] <- NA
        }
    }
    w
}

# Clean dN/dS values for fel1 and fel2 (detect NAs)
clean_dnds_fel <- function(df.fel, numcol)
{   
    w <- df.fel$dN.dS
    for (i in 1:numcol)
    {
        if (df.fel$LRT[i] == 0 & df.fel$p.value[i] == 1 & (w[i] == 0 | w[i] == 1)) #uninformative
        {
            w[i] <- NA
        }
        else if (w[i] >= max_threshold | w[i] == Inf)
        {
            w[i] <- NA
        }
    }
    w
} 
numcol <- 100
ntaxa <- 7:11
branch_lengths <- c(0.0025, 0.01, 0.04, 0.16, 0.64)
nreps <- 50


# Initialize results data frame
dat <- data.frame(       "type"     = factor(),
                         "ntaxa"     = numeric(),
                         "bl"    = numeric(),
                         "site"  = numeric(),
                         "rep"  = numeric(),
                         "ncount"= numeric(),
                         "scount"   = numeric())
levels(dat$type) <- as.factor(c("bias", "nobias"))

for (type in c("gtr", "bias_gtr")){
for (n in ntaxa){
    print(n)
    for (bl in branch_lengths){
        print(bl)
        for (repl in 1:nreps){

            name <- paste0("../results/raw_results/rep", repl, "_n", n, "_bl", bl, "_", type, "_SLAC_GTR.txt")
            counted <- read.table(name, header=T)        
            ns = counted$ObservedNSChanges
            ss = counted$ObservedSChanges
	
            # create data frame per method  
            temp  <- data.frame("type" = type, "ntaxa" = rep(2^n, numcol), "bl" = rep(bl, numcol), "site" = 1:numcol, "rep" = repl, "ncount" = ns, "scount" = ss )    
 
            # Merge
            dat <- rbind(dat, temp)
        }
    }
}}
write_csv(dat, "slac_counts.csv", delim = ",")
#dat %>% group_by(ntaxa, bl, rep) %>% mutate(ratio_ns = ncount/scount) %>% na.omit() %>% summarize(mean_ratio = mean(ratio_ns)) %>% filter(!is.infinite(mean_ratio)) -> ns_count_ratio

dat %>% group_by(ntaxa, bl, site, type) %>% mutate(ratio_ns = ncount/scount) %>% na.omit() %>% filter(!is.infinite(ratio_ns)) %>% summarize(mean_ratio = mean(ratio_ns))  -> ns_count_ratio


ggplot(ns_count_ratio, aes(x = as.factor(bl), y = mean_ratio, fill = as.factor(ntaxa))) + geom_violin(scale="width") + geom_hline(yintercept=1) + facet_grid(~type) + theme_bw()

# 
# summary(lm(mean_ratio ~ ntaxa+bl+type, data=ns_count_ratio))
# 
# Call:
# lm(formula = mean_ratio ~ ntaxa + bl + type, data = ns_count_ratio)
# 
# Residuals:
#      Min       1Q   Median       3Q      Max 
# -1.59881 -0.41463  0.04511  0.44805  2.41683 
# 
# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  1.518e+00  2.508e-02  60.526   <2e-16 ***
# ntaxa        3.186e-04  1.845e-05  17.268   <2e-16 ***
# bl           4.844e-01  5.335e-02   9.080   <2e-16 ***
# typebias_gtr 3.879e-02  2.577e-02   1.506    0.132    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.6441 on 2496 degrees of freedom
# Multiple R-squared:  0.133,	Adjusted R-squared:  0.132 
# F-statistic: 127.6 on 3 and 2496 DF,  p-value: < 2.2e-16
