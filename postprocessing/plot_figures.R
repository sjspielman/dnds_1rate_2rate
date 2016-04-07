require(dplyr)
require(tidyr)
require(grid)
require(cowplot)
require(readr)
source("summary_functions.R")

PLOTDIR <- "figures/"
method_order  <- c("SLAC1", "FEL1", "FUBAR1", "SLAC2", "FEL2", "FUBAR2")
method_colors <- c("deeppink3", "red", "darkred", "skyblue", "dodgerblue", "blue")
felfubar_colors <- c("red", "darkred", "dodgerblue", "blue")
dn_ds_colors <- c("red", "dodgerblue")

### Read in and process data for simulations performed along balanced trees
dnds.sum1 <- read.csv("dnds_summary_equalpi.csv")
dnds.sum2 <- read.csv("dnds_summary_unequalpi.csv")
dnds.sum1$pi <- "equalpi"
dnds.sum2$pi <- "unequalpi"
dnds.sum <- rbind(dnds.sum1, dnds.sum2)

dn.ds.sum1 <- read.csv("dn_ds_summary_equalpi.csv")
dn.ds.sum2 <- read.csv("dn_ds_summary_unequalpi.csv")
dn.ds.sum1$pi <- "equalpi"
dn.ds.sum2$pi <- "unequalpi"
dn.ds.sum <- rbind(dn.ds.sum1, dn.ds.sum2)


dnds.sum$method <- factor(dnds.sum$method, levels = method_order)
dn.ds.sum$method <- factor(dn.ds.sum$method, levels = method_order)
dnds.sum$pi <- factor(dnds.sum$mu, levels = c("equalpi", "unequalpi"))
dn.ds.sum$method <- factor(dn.ds.sum$method, levels = method_order)
dn.ds.sum$pi <- factor(dn.ds.sum$mu, levels = c("equalpi", "unequalpi"))

dnds.sum.mean <- dnds.sum %>% group_by(method, ntaxa, bl, type, pi) %>% summarize(r = mean(r), estbias = mean(estbias), rmsd = mean(rmsd))
dn.ds.sum.mean <- dn.ds.sum %>% group_by(method, ntaxa, bl, parameter, type, pi) %>% summarize(r = mean(r), estbias = mean(estbias), rmsd = mean(rmsd))
dnds.sum.mean$method <- factor(dnds.sum.mean$method, levels = method_order)
dn.ds.sum.mean$method <- factor(dn.ds.sum.mean$method, levels = method_order)
treelen.dat <- dnds.sum %>% mutate(treelen = 2*(ntaxa-1)*bl) %>% filter(method == "FEL1", treelen > 162, treelen < 164) %>% na.omit()

## Read in data of actual counted substitutions
equal.count.nobias <- read.csv("substitution_counts_equalpi_nobias.csv")
equal.count.nobias <- read.csv("substitution_counts_equalpi_bias.csv")
unequal.count.bias <- read.csv("substitution_counts_unequalpi_nobias.csv")
unequal.count.nobias <- read.csv("substitution_counts_unequalpi_bias.csv")

### Read in and factorize linear model results
comp.order <- c("FUBAR1 - SLAC1", "FEL1 - SLAC1", "FEL1 - FUBAR1", "SLAC1 - SLAC2", "FUBAR1 - FUBAR2", "FEL1 - FEL2")
linmodels <- read.csv("linear_model_results.csv") # Note that this file was manually created generated using information found in the file "build_linear_models.R"
linmodels <- filter(linmodels, comp %in% comp.order)
linmodels$model <- factor(linmodels$model, levels = c("r", "rmsd"))
linmodels$comp <- factor(linmodels$comp, levels = comp.order)
linmodels$type <- factor(linmodels$type, levels=c("nobias", "bias"), labels=c("No codon bias", "Codon bias"))
linmodels$pi   <- factor(linmodels$mutype, levels=c("unequalpi", "equalpi"), labels=c("Unequal F", "Equal F"))

### Frequently used legend
dn.ds.legend.grob <- dn.ds.sum %>% filter(pi == "equalpi", type == "bias", method == "FUBAR2") %>% ggplot(aes(x = factor(ntaxa), y = r, fill = parameter)) + geom_violin(scale = "width") + scale_fill_manual(values = dn_ds_colors, name = "Parameter", labels = c("dN    ", "dS")) + theme(legend.position = "bottom", legend.key.size = unit(.3, "cm"))
grobs <- ggplotGrob(dn.ds.legend.grob)$grobs
dn.ds.legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]



### Read in and process data for simulations performed along real trees
realdat <- read_csv("full_results_realtrees.csv")
realdat.sum <- realdat %>% filter(method %in% c("FUBAR1", "FUBAR2")) %>% summarize_dnds_real()
realdat.sum$method <- factor(realdat.sum$method, levels = c("FUBAR1", "FUBAR2"))
realdat.sum$dataset <- factor(realdat.sum$dataset, levels = c("amine", "vertrho", "camelid", "h3", "hivrt"))
realdat.sum$type <- factor(realdat.sum$type, levels=c("nobias", "bias"), labels=c("No codon bias", "Codon bias"))
realdat %>% filter(method %in% c("FUBAR2"), type == "bias") %>% summarize_dn_real() -> dn
realdat %>% filter(method %in% c("FUBAR2"), type == "bias") %>% summarize_ds_real() -> ds
dn.ds.real <- rbind(dn,ds)  
dn.ds.real$parameter <- factor(dn.ds.real$parameter, levels=c("dn", "ds"), labels=c("dN", "dS"))
dn.ds.real$dataset <- factor(dn.ds.real$dataset, levels = c("amine", "vertrho", "camelid", "h3", "hivrt"))




####################################################################################################
############ Lineplots of dN/dS correlations and RMSD for balanced tree simulations ################
####################################################################################################

theme_set(theme_cowplot() + theme(axis.text.y = element_text(size = 12), 
                                  axis.text.x = element_text(size = 9, angle=30),  
                                  axis.title = element_text(size = 12), 
                                  plot.title = element_text(size=13),
                                  panel.border = element_rect(size = 0.5), 
                                  panel.margin = unit(0.75, "lines"), 
                                  strip.background = element_rect(fill="white"), 
                                  strip.text = element_text(size=11), 
                                  legend.text = element_text(size=10), 
                                  legend.title = element_text(size=11)))

r.nobias.gtr <- dnds.sum.mean %>% filter(pi == "equalpi", type == "nobias") %>% ggplot(aes(x = factor(ntaxa), y = r, group = method, color=method)) +  geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 1), breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) + xlab("Number of Taxa") + ylab("Correlation") + scale_color_manual(values = method_colors, name = "Inference Method") + background_grid("xy") + ggtitle("No codon bias")
r.bias.gtr <- dnds.sum.mean %>% filter(pi == "equalpi", type == "bias") %>% ggplot(aes(x = factor(ntaxa), y = r, group = method, color=method)) +  geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 1), breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) + xlab("Number of Taxa") + ylab("Correlation") + scale_color_manual(values = method_colors, name = "Inference Method") + background_grid("xy") + ggtitle("Codon bias")
corr.lineplots <- plot_grid(r.nobias.equalpi, r.bias.equalpi, nrow=2, labels=c("A", "B"))
save_plot(paste0(PLOTDIR, "correlation_lineplots_equalpi.pdf"), corr.lineplots, base_width = 10, base_height=5)

r.nobias.hky <- dnds.sum.mean %>% filter(pi == "unequalpi", type == "nobias") %>% ggplot(aes(x = factor(ntaxa), y = r, group = method, color=method)) +  geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 1), breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) + xlab("Number of Taxa") + ylab("Correlation") + scale_color_manual(values = method_colors, name = "Inference Method") + background_grid("xy") + ggtitle("No codon bias")
r.bias.hky <- dnds.sum.mean %>% filter(pi == "unequalpi", type == "bias") %>% ggplot(aes(x = factor(ntaxa), y = r, group = method, color=method)) +  geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 1), breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) + xlab("Number of Taxa") + ylab("Correlation") + scale_color_manual(values = method_colors, name = "Inference Method") + background_grid("xy") + ggtitle("Codon bias")
corr.lineplots <- plot_grid(r.nobias.unequalpi, r.bias.unequalpi, nrow=2, labels=c("A", "B"))
save_plot(paste0(PLOTDIR, "correlation_lineplots_equalpi.pdf"), corr.lineplots, base_width = 10, base_height=5)



theme_set(theme_cowplot() + theme(axis.text = element_text(size = 10), 
                                  axis.title = element_text(size = 11), 
                                  plot.title = element_text(size=11, face="bold"),
                                  panel.border = element_rect(size = 0.5), 
                                  panel.margin = unit(0.75, "lines"), 
                                  strip.text = element_text(size=11), 
                                  strip.background = element_rect(fill="white"), 
                                  legend.text = element_text(size=10), 
                                  legend.title = element_text(size=11)))

rmsd.nobias.gtr <- dnds.sum %>% filter(pi == "equalpi", type == "nobias", bl >=0.04, ntaxa >= 256, method %in% c("FEL1", "FEL2", "FUBAR1", "FUBAR2"), rmsd<=100) %>% ggplot(aes(x = factor(ntaxa), y = rmsd, fill=method)) + geom_violin(scale="width", lwd=0.2) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 1.3 )) + xlab("Number of Taxa") + ylab("RMSD") + scale_fill_manual(values = felfubar_colors, name = "Inference Method") + background_grid("xy") + ggtitle("No codon bias")
rmsd.bias.gtr <- dnds.sum %>% filter(pi == "equalpi", type == "bias", bl >=0.04, ntaxa >= 256, method %in% c("FEL1", "FEL2", "FUBAR1", "FUBAR2"), rmsd<=100) %>% ggplot(aes(x = factor(ntaxa), y = rmsd, fill=method)) +  geom_violin(scale="width", lwd=0.2) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 2)) + xlab("Number of Taxa") + ylab("RMSD") + scale_fill_manual(values = felfubar_colors, name = "Inference Method") + background_grid("xy") + ggtitle("Codon bias")
rmsd.violins.gtr <- plot_grid(rmsd.nobias.equalpi, rmsd.bias.equalpi, nrow=2, labels=c("A", "B"))
save_plot(paste0(PLOTDIR, "rmsd_violins_equalpi.pdf"), rmsd.violins.gtr, base_width = 9, base_height=5)

rmsd.nobias.hky <- dnds.sum %>% filter(pi == "unequalpi", type == "nobias", bl >=0.04, ntaxa >= 256, method %in% c("FEL1", "FEL2", "FUBAR1", "FUBAR2"), rmsd<=100) %>% ggplot(aes(x = factor(ntaxa), y = rmsd, fill=method)) + geom_violin(scale="width", lwd=0.2) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 1.3 )) + xlab("Number of Taxa") + ylab("RMSD") + scale_fill_manual(values = felfubar_colors, name = "Inference Method") + background_grid("xy") + ggtitle("No codon bias")
rmsd.bias.hky <- dnds.sum %>% filter(pi == "unequalpi", type == "bias", bl >=0.04, ntaxa >= 256, method %in% c("FEL1", "FEL2", "FUBAR1", "FUBAR2"), rmsd<=100) %>% ggplot(aes(x = factor(ntaxa), y = rmsd, fill=method)) +  geom_violin(scale="width", lwd=0.2) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 2)) + xlab("Number of Taxa") + ylab("RMSD") + scale_fill_manual(values = felfubar_colors, name = "Inference Method") + background_grid("xy") + ggtitle("Codon bias")
rmsd.violins.hky <- plot_grid(rmsd.nobias.unequalpi, rmsd.bias.unequalpi, nrow=2, labels=c("A", "B"))
save_plot(paste0(PLOTDIR, "rmsd_violins_unequalpi.pdf"), rmsd.violins.unequalpi, base_width = 9, base_height=5)





####################################################################################################
######### Violin plots  of dN & dS correlations and RMSD for balanced tree simulations #############
####################################################################################################

theme_set(theme_cowplot() + theme(axis.text.y = element_text(size = 10),
                                  axis.text.x = element_text(size = 9, angle=30),
                                  axis.title = element_text(size = 11), 
                                  legend.title = element_text(size = 10),
                                  legend.text = element_text(size=9),
                                  strip.background = element_rect(fill="white"), 
                                  strip.text = element_text(size=10), 
                                  panel.border = element_rect(size = 0.5), 
                                  panel.margin = unit(1.0, "lines")))

param.r.fel2.sub    <- dn.ds.sum %>% filter(mu == "unequalpi", type == "bias", method == "FEL2", bl >= 0.01, ntaxa<=1024) %>% ggplot(aes(x = factor(ntaxa), y = r, fill = parameter)) + geom_violin(scale = "width") + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(-0.1, 1), breaks = c(-0.1, 0.1, 0.3, 0.5, 0.7, 0.9)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("Correlation")  + scale_fill_manual(values = dn_ds_colors) + theme(legend.position = "none")
param.rmsd.fel2y.sub <- dn.ds.sum %>% filter(mu == "unequalpi", type == "bias", method == "FEL2", bl >= 0.01, ntaxa<=1024, rmsd <= 3) %>% ggplot(aes(x = factor(ntaxa), y = rmsd, fill = parameter)) + geom_violin(scale = "width") + facet_grid(~bl, scales = "free_x") + coord_cartesian(ylim=c(0,3)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("RMSD") + scale_fill_manual(values = dn_ds_colors) + theme(legend.position = "none")

dn.ds.fel2.sub.legend <- plot_grid(param.r.fel2.sub , param.rmsd.fel2.sub, dn.ds.legend, nrow=3, labels = c("A", "B"), rel_heights=c(1,1,0.1))
save_plot(paste0(PLOTDIR, "dn_ds_corr_violin_fel2_subset.pdf"), dn.ds.fel2.sub.legend, base_width = 6.5, base_height=4.5)


param.r.fel2.gtr.sub    <- dn.ds.sum %>% filter(mu == "GTR", type == "bias", method == "FEL2", bl >= 0.01, ntaxa<=1024) %>% ggplot(aes(x = factor(ntaxa), y = r, fill = parameter)) + geom_violin(scale = "width") + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(-0.1, 1), breaks = c(-0.1, 0.1, 0.3, 0.5, 0.7, 0.9)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("Correlation")  + scale_fill_manual(values = dn_ds_colors) + theme(legend.position = "none")
param.rmsd.fel2.gtr.sub <- dn.ds.sum %>% filter(mu == "GTR", type == "bias", method == "FEL2", bl >= 0.01, ntaxa<=1024, rmsd <= 3) %>% ggplot(aes(x = factor(ntaxa), y = rmsd, fill = parameter)) + geom_violin(scale = "width") + facet_grid(~bl, scales = "free_x") + coord_cartesian(ylim=c(0,3)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("RMSD") + scale_fill_manual(values = dn_ds_colors) + theme(legend.position = "none")

dn.ds.fel2.gtr.sub.legend <- plot_grid(param.r.fel2.gtr.sub , param.rmsd.fel2.gtr.sub, dn.ds.legend, nrow=3, labels = c("A", "B"), rel_heights=c(1,1,0.1))
save_plot(paste0(PLOTDIR, "dn_ds_corr_violin_fel2_gtr_subset.pdf"), dn.ds.fel2.gtr.sub.legend, base_width = 6.5, base_height=4.5)


theme_set(theme_cowplot() + theme(axis.text.y = element_text(size = 10),
                                  axis.text.x = element_text(size = 7, angle=30),
                                  axis.title = element_text(size = 12), 
                                  legend.title = element_text(size = 10),
                                  legend.text = element_text(size=9),
                                  strip.background = element_rect(fill="white"), 
                                  strip.text = element_text(size=10), 
                                  panel.border = element_rect(size = 0.5), 
                                  panel.margin = unit(1.0, "lines")))




param.r.fubar2.gtr <- dn.ds.sum %>% filter(mu == "GTR", type == "bias", method == "FUBAR2") %>% ggplot(aes(x = factor(ntaxa), y = r, fill = parameter)) + geom_violin(scale = "width") + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(-0.1, 1), breaks = c(-0.1, 0.1, 0.3, 0.5, 0.7, 0.9)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("Correlation") + scale_fill_manual(values = dn_ds_colors) + ggtitle("GTR-bias") + theme(legend.position = "none")
param.r.fubar2.hky <- dn.ds.sum %>% filter(mu == "HKY", type == "bias", method == "FUBAR2") %>% ggplot(aes(x = factor(ntaxa), y = r, fill = parameter)) + geom_violin(scale = "width") + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(-0.1, 1), breaks = c(-0.1, 0.1, 0.3, 0.5, 0.7, 0.9)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("Correlation") + scale_fill_manual(values = dn_ds_colors) + ggtitle("HKY-bias") + theme(legend.position = "none")
param.rmsd.fubar2.hky <- dn.ds.sum %>% filter(mu == "HKY", type == "bias", method == "FUBAR2") %>% ggplot(aes(x = factor(ntaxa), y = rmsd, fill = parameter)) + geom_violin(scale = "width") + facet_grid(~bl, scales = "free_x") + coord_cartesian(ylim=c(0,5)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("RMSD") + scale_fill_manual(values = dn_ds_colors) + ggtitle("HKY-bias") + theme(legend.position = "none")
param.rmsd.fubar2.gtr <- dn.ds.sum %>% filter(mu == "GTR", type == "bias", method == "FUBAR2") %>% ggplot(aes(x = factor(ntaxa), y = rmsd, fill = parameter)) + geom_violin(scale = "width") + facet_grid(~bl, scales = "free_x") + coord_cartesian(ylim=c(0,5)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("RMSD") + scale_fill_manual(values = dn_ds_colors) + ggtitle("GTR-bias") + theme(legend.position = "none")
param.r.fubar2 <- plot_grid(param.r.fubar2.hky, param.r.fubar2.gtr, nrow=1, labels=c("A", "B"))
param.rmsd.fubar2 <- plot_grid(param.rmsd.fubar2.hky, param.rmsd.fubar2.gtr, nrow=1, labels=c("C", "D"))
param.r.rmsd.fubar2 <- plot_grid(param.r.fubar2, param.rmsd.fubar2, dn.ds.legend, nrow=3, rel_heights=c(1,1,0.1))
save_plot(paste0(PLOTDIR, "dn_ds_r_rmsd_fubar2.pdf"), param.r.rmsd.fubar2, base_width = 12, base_height=5)



param.r.fel2.gtr <- dn.ds.sum %>% filter(mu == "GTR", type == "bias", method == "FEL2") %>% ggplot(aes(x = factor(ntaxa), y = r, fill = parameter)) + geom_violin(scale = "width") + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(-0.1, 1), breaks = c(-0.1, 0.1, 0.3, 0.5, 0.7, 0.9)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("Correlation") + scale_fill_manual(values = dn_ds_colors) + ggtitle("GTR-bias") + theme(legend.position = "none")
param.r.fel2.hky <- dn.ds.sum %>% filter(mu == "HKY", type == "bias", method == "FEL2") %>% ggplot(aes(x = factor(ntaxa), y = r, fill = parameter)) + geom_violin(scale = "width") + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(-0.1, 1), breaks = c(-0.1, 0.1, 0.3, 0.5, 0.7, 0.9)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("Correlation") + scale_fill_manual(values = dn_ds_colors) + ggtitle("HKY-bias") + theme(legend.position = "none")
param.rmsd.fel2.hky <- dn.ds.sum %>% filter(mu == "HKY", type == "bias", method == "FEL2") %>% ggplot(aes(x = factor(ntaxa), y = rmsd, fill = parameter)) + geom_violin(scale = "width") + facet_grid(~bl, scales = "free_x") + coord_cartesian(ylim=c(0,5)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("RMSD") + scale_fill_manual(values = dn_ds_colors) + ggtitle("HKY-bias") + theme(legend.position = "none")
param.rmsd.fel2.gtr <- dn.ds.sum %>% filter(mu == "GTR", type == "bias", method == "FEL2") %>% ggplot(aes(x = factor(ntaxa), y = rmsd, fill = parameter)) + geom_violin(scale = "width") + facet_grid(~bl, scales = "free_x") + coord_cartesian(ylim=c(0,5)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("RMSD") + scale_fill_manual(values = dn_ds_colors) + ggtitle("GTR-bias") + theme(legend.position = "none")
param.r.fel2 <- plot_grid(param.r.fel2.hky, param.r.fel2.gtr, nrow=1, labels=c("A", "B"))
param.rmsd.fel2 <- plot_grid(param.rmsd.fel2.hky, param.rmsd.fel2.gtr, nrow=1, labels=c("C", "D"))
param.r.rmsd.fel2 <- plot_grid(param.r.fel2, param.rmsd.fel2, dn.ds.legend, nrow=3, rel_heights=c(1,1,0.1))
save_plot(paste0(PLOTDIR, "dn_ds_r_rmsd_fel2.pdf"), param.r.rmsd.fel2, base_width = 12, base_height=5)



####################################################################################################
###################################### Linear models plot ##########################################
####################################################################################################
sig_colors <- c("grey60", "black")
ptsize = 1.35
linesize = 0.8
theme_set(theme_cowplot() + theme(axis.text.x = element_text(size=10),
                                  axis.text.y = element_text(size=9),
                                  axis.title.x = element_text(size=11), 
                                  axis.title.y = element_text(size=12), 
                                  panel.border = element_rect(size = 0.5), 
                                  panel.margin = unit(1.0, "lines"), 
                                  strip.background = element_rect(fill="white"), 
                                  strip.text = element_text(size=11),
                                  legend.position = "none"))


linmodels %>% filter(model == "r") %>% arrange(comp) %>%
  ggplot(aes(x = coeff, y = comp, color = sig)) + geom_point(size=ptsize) + 
  geom_segment(aes(x=lowerCI,xend=upperCI,y=comp,yend=comp),size=linesize) + 
  geom_vline(xintercept=0, size=0.5) + 
  facet_grid(~type) + background_grid() +
  xlab("Average Correlation Difference") + ylab("Comparison") + 
  scale_x_continuous(limits=c(-0.05, 0.3)) +  
  scale_color_manual(values = sig_colors) -> linmodel.r

linmodels %>% filter(model == "rmsd") %>% arrange(comp) %>%
  ggplot(aes(x = coeff, y = comp, color = sig)) + geom_point(size=ptsize) + 
  geom_segment(aes(x=lowerCI,xend=upperCI,y=comp,yend=comp),size=linesize) + 
  geom_vline(xintercept=0, size=0.5) + 
  facet_grid(~type) + background_grid() +
  xlab("Average RMSD Difference") + ylab("Comparison") + 
  scale_x_continuous(limits=c(-0.3,0.1)) +  
  scale_color_manual(values = sig_colors) -> linmodel.rmsd

linmodel.plots <- plot_grid(linmodel.r, linmodel.rmsd, nrow=2, labels=c("A", "B"))
save_plot(paste0(PLOTDIR, "linmodels.pdf"), linmodel.plots, base_width = 6.5, base_height = 4.5)

##########################################################################################


theme_set(theme_cowplot() + theme(axis.text = element_text(size=13),
                                  axis.title = element_text(size=14)))



bias.gtr.count %>% group_by(ntaxa, bl, rep) %>% mutate(ratio_ns = ncount/scount) %>% na.omit() %>% filter(!is.infinite(ratio_ns)) %>% summarize(mean_ratio = mean(ratio_ns))  -> ns_count_ratio
count.ratio.box <- ggplot(ns_count_ratio, aes(x = as.factor(bl), y = mean_ratio, fill = as.factor(ntaxa))) + geom_boxplot(outlier.size = 0.75) + geom_hline(yintercept=1) + scale_fill_brewer(palette = "YlGn", name = "Number of Taxa") + xlab("Branch Length") + ylab("Mean N/S count ratio") +  scale_y_continuous(limits=c(0,4))
save_plot(paste0(PLOTDIR, "counted_ratio_bias.pdf"), count.ratio.box, base_width = 8, base_height = 5)

gtr.count %>% group_by(ntaxa, bl, rep) %>% mutate(ratio_ns = ncount/scount) %>% na.omit() %>% filter(!is.infinite(ratio_ns)) %>% summarize(mean_ratio = mean(ratio_ns))  -> ns_count_ratio
count.ratio.box <- ggplot(ns_count_ratio, aes(x = as.factor(bl), y = mean_ratio, fill = as.factor(ntaxa))) + geom_boxplot(outlier.size = 0.75) + geom_hline(yintercept=1) + scale_fill_brewer(palette = "YlGn", name = "Number of Taxa") + xlab("Branch Length") + ylab("Mean N/S count ratio") + scale_y_continuous(limits=c(0,4))
save_plot(paste0(PLOTDIR, "counted_ratio_nobias.pdf"), count.ratio.box, base_width = 8, base_height = 5)


####################################################################################################
######################### Correlations and RMSD for real tree simulations ##########################
####################################################################################################

theme_set(theme_cowplot() + theme(axis.text    = element_text(size = 10),
                                  axis.title   = element_text(size = 12),
                                  panel.border = element_rect(size = 0.5), 
                                  panel.margin = unit(0.8, "lines"), 
                                  strip.background = element_rect(fill="white"), 
                                  strip.text = element_text(size=11), 
                                  legend.text = element_text(size=10), 
                                  legend.title = element_text(size=11)))



r.dnds <- realdat.sum %>% ggplot(aes(x = dataset, fill = method, y = r)) + geom_boxplot(position=position_dodge(0.77), outlier.size = 0.75, size=0.4) + facet_wrap(~type) + scale_fill_manual(values=dn_ds_colors, name = "Inference Method") +  xlab("Dataset") + ylab("Correlation") + scale_y_continuous(limits=c(0,1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.))
rmsd.dnds <- realdat.sum %>% ggplot(aes(x = dataset, fill = method, y = rmsd)) + geom_boxplot(position=position_dodge(0.77), outlier.size = 0.75, size=0.4) + facet_wrap(~type) + scale_fill_manual(values=dn_ds_colors, name = "Inference Method") + xlab("Dataset") + ylab("RMSD")
#r.dn.ds <- dn.ds.real %>% ggplot(aes(x = dataset, fill = parameter, y = r)) + geom_boxplot(outlier.size = 0.75, size=0.4) + scale_fill_manual(values = fill_colors, name = "Parameter") + ylab("Correlation") + xlab("Dataset")
#rmsd.dn.ds <- dn.ds.real %>% ggplot(aes(x = dataset, fill = parameter, y = rmsd)) + geom_boxplot(outlier.size = 0.75, size=0.4) + scale_fill_manual(values = fill_colors, name = "Parameter") + ylab("RMSD") + xlab("Dataset")

real.plots <- plot_grid(r.dnds, rmsd.dnds, nrow=2, labels=c("A", "B"))
save_plot(paste0(PLOTDIR, "real_r_rmsd.pdf"), real.plots, base_width = 7.5, base_height = 4.5)



####################################################################################################
########################### Scatterplots of dN vs. dN for FEL and FUBAR ############################
####################################################################################################

full.nobias <- read_csv("full_results_gtr_nobias.csv")
full.bias <- read_csv("full_results_gtr_bias.csv")
full.nobias %>% select(-true, -truedn, -trueds, -dnds, -ds) %>% 
  na.omit() %>%
  group_by(ntaxa, bl, site, method) %>%
  summarize(mean.dn = mean(dn)) %>%
  spread(method, mean.dn) -> nobias.compare.dn
full.bias  %>% select(-true, -truedn, -trueds, -dnds, -ds) %>% 
  na.omit() %>%
  group_by(ntaxa, bl, site, method) %>%
  summarize(mean.dn = mean(dn)) %>%
  spread(method, mean.dn) -> bias.compare.dn

theme_set(theme_cowplot() + theme(axis.text = element_text(size = 10), 
                                  axis.title = element_text(size = 13), 
                                  strip.text = element_text(size=12),
                                  panel.border = element_rect(size = 0.75), 
                                  panel.margin = unit(1.0, "lines")))

nobias.fubar.dn <- nobias.compare.dn %>% ggplot(aes(x = FUBAR1, y = FUBAR2)) + geom_point() + facet_grid(ntaxa~bl) + geom_abline(slope=1, intercept=0, color="mediumvioletred", size=0.7) + scale_x_continuous(limits=c(0,1), breaks=c(0, 0.5, 1)) + scale_y_continuous(limits=c(0,1), breaks=c(0, 0.5, 1)) + xlab("dN, FUBAR1") + ylab("dN, FUBAR2")
nobias.fel.dn <-nobias.compare.dn %>% ggplot(aes(x = FEL1, y = FEL2)) + geom_point() + facet_grid(ntaxa~bl) + geom_abline(slope=1, intercept=0, color="mediumvioletred", size=0.7) + scale_x_continuous(limits=c(0,1), breaks=c(0, 0.5, 1)) + scale_y_continuous(limits=c(0,1), breaks=c(0, 0.5, 1)) + xlab("dN, FEL1") + ylab("dN, FEL2")

bias.fubar.dn <- bias.compare.dn %>% ggplot(aes(x = FUBAR1, y = FUBAR2)) + geom_point() + facet_grid(ntaxa~bl) + geom_abline(slope=1, intercept=0, color="mediumvioletred", size=0.7) + scale_x_continuous(limits=c(0,1), breaks=c(0, 0.5, 1)) + scale_y_continuous(limits=c(0,1), breaks=c(0, 0.5, 1)) + xlab("dN, FEL1") + ylab("dN, FEL2")
bias.fel.dn <-bias.compare.dn %>% ggplot(aes(x = FEL1, y = FEL2)) + geom_point() + facet_grid(ntaxa~bl) + geom_abline(slope=1, intercept=0, color="mediumvioletred", size=0.7) +scale_x_continuous(limits=c(0,1), breaks=c(0, 0.5, 1)) + scale_y_continuous(limits=c(0,1), breaks=c(0, 0.5, 1)) + xlab("dN, FEL1") + ylab("dN, FEL2")

dn.same.grid <- plot_grid(nobias.fel.dn, nobias.fubar.dn, bias.fel.dn, nobias.fel.dn, nrow=2, labels=c("A", "B", "C", "D"), label_size=20)
save_plot(paste0(PLOTDIR, "samedn.pdf"), dn.same.grid, base_width = 10, base_height=9)





####################################################################################################
############################## Divergence vs. number of taxa #######################################
####################################################################################################

theme_set(theme_cowplot() + theme(axis.text.x = element_text(size = 9), 
                                  axis.text.y = element_text(size = 11),  
                                  axis.title = element_text(size = 12), 
                                  legend.text = element_text(size=10)))

legend.grob <- treelen.dat %>% ggplot(aes(x = as.factor(bl), y = r, fill = as.factor(type))) + geom_boxplot(position = position_dodge(0.9)) + xlab("Parameterization") + scale_fill_manual(name = "", labels = c("No Codon bias     ", "Codon bias"), values = c("grey40", "grey80")) + theme(legend.position = "bottom")
grobs <- ggplotGrob(legend.grob)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

ntaxa.bl.r <- treelen.dat %>% ggplot(aes(x = as.factor(bl), y = r, fill = as.factor(type))) + geom_boxplot(position = position_dodge(0.75), lwd=0.3, outlier.size=0.8) + xlab("Parameterization") + scale_x_discrete(labels = c("N = 2048\nB = 0.04","N = 512\nB = 0.16","N = 128\nB = 0.64")) + ylab("Correlation") + scale_y_continuous(limits=c(0.45, 0.95)) + scale_fill_manual(values = c("grey30", "grey80")) + theme(legend.position = "none")
ntaxa.bl.rmsd <- treelen.dat %>% ggplot(aes(x = as.factor(bl), y = rmsd, fill = as.factor(type))) + geom_boxplot(position = position_dodge(0.75), lwd=0.3, outlier.size=0.8) + xlab("Parameterization") + scale_x_discrete(labels = c("N = 2048\nB = 0.04","N = 512\nB = 0.16","N = 128\nB = 0.64")) + ylab("RMSD") + scale_y_continuous(limits=c(0,0.75)) + scale_fill_manual(values = c("grey30", "grey80")) + theme(legend.position = "none")
ntaxa.bl.r.rmsd <- plot_grid(ntaxa.bl.r, ntaxa.bl.rmsd, nrow=1, labels=c("A", "B"))
ntaxa.bl.r.rmsd.legend <- plot_grid(ntaxa.bl.r.rmsd, legend, nrow=2, rel_heights=c(1,0.1))
save_plot(paste0(PLOTDIR, "ntaxa_bl_r_rmsd.pdf"), ntaxa.bl.r.rmsd.legend, base_width=7, base_height=3)