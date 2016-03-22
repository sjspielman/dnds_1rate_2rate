require(dplyr)
require(tidyr)
require(grid)
require(cowplot)
require(readr)



method_order  <- c("SLAC1", "FEL1", "FUBAR1", "SLAC2", "FEL2", "FUBAR2")
method_colors <- c("deeppink3", "red", "darkred", "skyblue", "dodgerblue", "blue")

dnds.sum <- read.csv("dnds_summary.csv")
dn.ds.sum <- read.csv("dn_ds_summary.csv")
dnds.sum$method <- factor(dnds.sum$method, levels = method_order)
dn.ds.sum$method <- factor(dn.ds.sum$method, levels = method_order)

dnds.sum.mean <- dnds.sum %>% group_by(method, ntaxa, bl, truetype, type) %>% summarize(r = mean(r), estbias = mean(estbias), rmsd = mean(rmsd))
dn.ds.sum.mean <- dn.ds.sum %>% group_by(method, ntaxa, bl, parameter, type) %>% summarize(r = mean(r), estbias = mean(estbias), rmsd = mean(rmsd))
dnds.sum.mean$method <- factor(dnds.sum.mean$method, levels = method_order)
dn.ds.sum.mean$method <- factor(dn.ds.sum.mean$method, levels = method_order)




r.nobias <- dnds.sum.mean %>% filter(type == "nobias", truetype == "true1", bl >0.0025, ntaxa<2048) %>% ggplot(aes(x = factor(ntaxa), y = r, group = method, color=method)) +  geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 1), breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) + xlab("Number of Taxa") + ylab("Pearson Correlation") + scale_color_manual(values = method_colors, name = "Inference Method")
r.bias <- dnds.sum.mean %>% filter(type == "bias", truetype == "true2", bl >0.0025, ntaxa<2048) %>% ggplot(aes(x = factor(ntaxa), y = r, group = method, color=method)) +  geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 1), breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) + xlab("Number of Taxa") + ylab("Pearson Correlation") + scale_color_manual(values = method_colors, name = "Inference Method")
plot_grid(r.nobias,r.bias,nrow=2)


rmsd.nobias <- dnds.sum.mean %>% filter(type == "nobias", truetype == "true1", bl >=0.04, ntaxa >= 512) %>% ggplot(aes(x = factor(ntaxa), y = rmsd, group = method, color=method)) +  geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 1.25)) + xlab("Number of Taxa") + ylab("RMSD") + scale_color_manual(values = method_colors, name = "Inference Method")
rmsd.bias <- dnds.sum.mean %>% filter(type == "bias", truetype == "true2", bl >=0.04, ntaxa >= 512) %>% ggplot(aes(x = factor(ntaxa), y = rmsd, group = method, color=method)) +  geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 1.25)) + xlab("Number of Taxa") + ylab("RMSD") + scale_color_manual(values = method_colors, name = "Inference Method")
plot_grid(rmsd.nobias,rmsd.bias,nrow=2)



############################################################################
######### dN is virtually the same between 1rate and 2rate models ##########
############################################################################
full.nobias <- read_csv("../results/processed_results/full_results_gtr_dataset.csv")
full.bias <- read_csv("../results/processed_results/full_results_gtr_dataset.csv")

full.nobias %>% select(-true1, -true2, -truedn, -trueds, -dnds, -ds) %>% 
               na.omit() %>%
               group_by(ntaxa, bl, site, method) %>%
               summarize(mean.dn = mean(dn)) %>%
               spread(method, mean.dn) -> nobias.compare.dn
nobias.fubar.dn <- nobias.compare.dn %>% ggplot(aes(x = FUBAR1, y = FUBAR2)) + geom_point() + facet_grid(ntaxa~bl) + geom_abline(slope=1, intercept=0, color="mediumvioletred") + scale_x_continuous(limits=c(0,1)) + scale_y_continuous(limits=c(0,1)) + xlab("dN, FUBAR1") + ylab("dN, FUBAR2")
nobias.fel.dn <-nobias.compare.dn %>% ggplot(aes(x = FEL1, y = FEL2)) + geom_point() + facet_grid(ntaxa~bl) + geom_abline(slope=1, intercept=0, color="mediumvioletred") + scale_x_continuous(limits=c(0,1)) + scale_y_continuous(limits=c(0,1)) + xlab("dN, FEL1") + ylab("dN, FEL2")

full.bias  %>% select(-true1, -true2, -truedn, -trueds, -dnds, -ds) %>% 
               na.omit() %>%
               group_by(ntaxa, bl, site, method) %>%
               summarize(mean.dn = mean(dn)) %>%
               spread(method, mean.dn) -> bias.compare.dn
bias.fubar.dn <- nobias.compare.dn %>% ggplot(aes(x = FUBAR1, y = FUBAR2)) + geom_point() + facet_grid(ntaxa~bl) + geom_abline(slope=1, intercept=0, color="mediumvioletred") + scale_x_continuous(limits=c(0,1)) + scale_y_continuous(limits=c(0,1)) + xlab("dN, FEL1") + ylab("dN, FEL2")
bias.fel.dn <-nobias.compare.dn %>% ggplot(aes(x = FEL1, y = FEL2)) + geom_point() + facet_grid(ntaxa~bl) + geom_abline(slope=1, intercept=0, color="mediumvioletred") + scale_x_continuous(limits=c(0,1)) + scale_y_continuous(limits=c(0,1)) + xlab("dN, FEL1") + ylab("dN, FEL2")



dn.r.fel2    <- dn.ds.sum %>% filter(type == "bias", method == "FEL2", bl >= 0.01, ntaxa<=1024) %>% ggplot(aes(x = factor(ntaxa), y = r, fill = parameter)) + geom_violin(scale = "width") + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(-0.1, 1), breaks = c(-0.1, 0.1, 0.3, 0.5, 0.7, 0.9)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("Pearson Correlation") + scale_fill_manual(values = c("grey80", "grey30"), name = "Parameter", labels = c("dN", "dS")) + theme(legend.position = "none")
dn.rmsd.fel2 <- dn.ds.sum %>% filter(type == "bias", method == "FEL2", bl >= 0.01, ntaxa<=1024, rmsd <= 2) %>% ggplot(aes(x = factor(ntaxa), y = rmsd, fill = parameter)) + geom_violin(scale = "width") + facet_grid(~bl, scales = "free_x") + coord_cartesian(ylim=c(0,2)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("RMSD") + scale_fill_manual(values = c("grey80", "grey30"), name = "Parameter", labels = c("dN", "dS"))+ theme(legend.position = "none")
plot_grid(dn.r.fel2 , dn.rmsd.fel2, nrow=2)


##########################################################################################
sig_colors <- c("grey60", "black")
comp.order <- c("FUBAR1 - SLAC1", "FEL1 - SLAC1", "FEL1 - FUBAR1", "SLAC1 - SLAC2", "FUBAR1 - FUBAR2", "FEL1 - FEL2")
theme_set(theme_cowplot() + theme(axis.title = element_text(size=15), 
                                  panel.border = element_rect(size = 0.5), 
                                  panel.margin = unit(0.75, "lines"), 
                                  strip.background = element_rect(fill="white"), 
                                  strip.text = element_text(size=12)))


linmodels <- read.csv("linear_model_results_correlation.csv") # Note that this file was manually created generated using information found in the file "build_linear_models.R"
linmodels2 <- linmodels %>% filter(comp %in% comp.order, type %in% c("nobias", "bias2"))
linmodels2$comp <- factor(linmodels2$comp, levels = comp.order)
linmodels2$type <- factor(linmodels2$type, levels=c("nobias", "bias2"), labels=c("No codon bias", "Codon bias"))


linmodels2 %>% arrange(comp) %>% 
           ggplot(aes(x = coeff, y = comp, color = sig)) + 
                  geom_point(size=2.5) + 
                  geom_segment(aes(x=lowerCI,xend=upperCI,y=comp,yend=comp),size=1) + 
                  geom_vline(xintercept=0, size=0.5) + 
                  facet_grid(~type) + 
                  xlab("Average Correlation Difference") + 
                  ylab("Method Comparison") + 
                  background_grid() + 
                  scale_x_continuous(limits=c(-0.05, 0.25)) +  
                  scale_color_manual(values = sig_colors, name = "", labels = c("NS", "S")) -> linmodel.r


linmodels <- read.csv("linear_model_results_rmsd.csv") # Note that this file was manually created generated using information found in the file "build_linear_models.R"
linmodels2 <- linmodels %>% filter(comp %in% comp.order, type %in% c("nobias", "bias2"))
linmodels2$comp <- factor(linmodels2$comp, levels = comp.order)
linmodels2$type <- factor(linmodels2$type, levels=c("nobias", "bias2"), labels=c("No codon bias", "Codon bias"))
linmodels2 %>% arrange(comp) %>% 
           ggplot(aes(x = coeff, y = comp, color = sig)) + 
                  geom_point(size=1.75) + 
                  geom_segment(aes(x=lowerCI,xend=upperCI,y=comp,yend=comp),size=0.8) + 
                  geom_vline(xintercept=0, size=0.5) + 
                  facet_grid(~type) + 
                  xlab("Average RMSD Difference") + 
                  ylab("Method Comparison") + 
                  background_grid() + 
                  scale_x_continuous(limits=c(-0.15, 0.1)) +  
                  scale_color_manual(values = sig_colors, name = "", labels = c("NS", "S")) -> linmodel.rmsd



##########################################################################################
gtr.bias.count <- read.csv("gtr_bias_count.csv")
gtr.bias.count %>% group_by(ntaxa, bl, rep) %>% mutate(ratio_ns = ncount/scount) %>% na.omit() %>% filter(!is.infinite(ratio_ns)) %>% summarize(mean_ratio = mean(ratio_ns))  -> ns_count_ratio
count.ratio.box <- ggplot(ns_count_ratio, aes(x = as.factor(bl), y = mean_ratio, fill = as.factor(ntaxa))) + geom_boxplot(outlier.size = 0.75) + geom_hline(yintercept=1) + scale_fill_brewer(palette = "YlOrRd", name = "Number of Taxa") + xlab("Branch Length") + ylab("Mean N:S count ratio")


##########################################################################################
source("summarize_results.R")
realdat <- read.csv("full_results_realtrees.csv")

realdat %>% filter(method %in% c("SLAC2", "FUBAR2"), type == "bias_gtr") %>% summarize_dn_real() -> dn
realdat %>% filter(method %in% c("SLAC2", "FUBAR2"), type == "bias_gtr") %>% summarize_ds_real() -> ds
rmsd.dn.ds.real.box<- rbind(dn,ds) %>% filter(method == "FUBAR2") %>% ggplot(aes(x = dataset, fill = parameter, y = rmsd)) + geom_boxplot()
r.dn.ds.real.box <- rbind(dn,ds) %>% filter(method == "FUBAR2") %>% ggplot(aes(x = dataset, fill = parameter, y = r)) + geom_boxplot()


realdat.sum <- realdat %>% summarize_dnds_real() 
realdat.sum$method <- factor(realdat.sum$method, levels = c("SLAC1", "FUBAR1", "SLAC2", "FUBAR2"))
realdat.sum$dataset <- factor(realdat.sum$dataset, levels = c("amine", "camelid", "vertrho", "h3", "hivrt"))
fill_colors <- c("red", "orange")
r    <- realdat.sum %>% filter(truetype == "true2", method %in% c("SLAC1", "SLAC2")) %>% ggplot(aes(x = dataset, fill = method, y = r)) + geom_boxplot() + facet_wrap(~type) + scale_fill_manual(values=fill_colors)
rmsd <- realdat.sum %>% filter(truetype == "true2", method %in% c("SLAC1", "SLAC2")) %>% ggplot(aes(x = dataset, fill = method, y = rmsd)) + geom_boxplot() + facet_wrap(~type) + scale_fill_manual(values=fill_colors)
real.r.rmsd <- plot_grid(r,rmsd, nrow=2)




