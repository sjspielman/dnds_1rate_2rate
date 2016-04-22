require(dplyr)
require(tidyr)
require(grid)
require(cowplot)
require(readr)
source("summary_functions.R")

PLOTDIR <- "figures/"
DATADIR <- "dataframes/"
method_order  <- c("SLAC1", "FEL1", "FUBAR1", "SLAC2", "FEL2", "FUBAR2")
method_colors <- c("deeppink3", "red", "darkred", "skyblue", "dodgerblue", "blue")
method_colors_nofel2 <- c("deeppink3", "red", "darkred", "skyblue", "blue")

felfubar_colors <- c("red", "darkred", "dodgerblue", "blue")
dn_ds_colors <- c("red", "dodgerblue")


theme_set(theme_cowplot() + theme(panel.border = element_rect(size = 0.5), 
                                  panel.margin = unit(0.75, "lines"), 
                                  strip.background = element_rect(fill="white"), 
                                  strip.text = element_text(size=12)))

### Read in and process data for simulations performed along balanced trees
dnds.sum.eq.nobias <- read.csv(paste0(DATADIR, "dnds_summary_equalpi_nobias.csv"))
dnds.sum.eq.bias <- read.csv(paste0(DATADIR, "dnds_summary_equalpi_bias.csv"))
dnds.sum.eq <- rbind(dnds.sum.eq.nobias, dnds.sum.eq.bias)
dnds.sum.eq$pi <- "equal"
dnds.sum.uneq.nobias <- read.csv(paste0(DATADIR, "dnds_summary_unequalpi_nobias.csv"))
dnds.sum.uneq.bias <- read.csv(paste0(DATADIR, "dnds_summary_unequalpi_bias.csv"))
dnds.sum.uneq <- rbind(dnds.sum.uneq.nobias, dnds.sum.uneq.bias)
dnds.sum.uneq$pi <- "unequal"
dnds.sum <- rbind(dnds.sum.eq, dnds.sum.uneq)
dnds.sum$method <- factor(dnds.sum$method, levels = method_order)
dnds.sum$pi <- factor(dnds.sum$pi, levels = c("equal", "unequal"))



dn.ds.sum.eq.nobias <- read.csv(paste0(DATADIR, "dn_ds_summary_equalpi_nobias.csv"))
dn.ds.sum.eq.bias <- read.csv(paste0(DATADIR, "dn_ds_summary_equalpi_bias.csv"))
dn.ds.sum.eq <- rbind(dn.ds.sum.eq.nobias, dn.ds.sum.eq.bias)
dn.ds.sum.eq$pi <- "equal"
dn.ds.sum.uneq.nobias <- read.csv(paste0(DATADIR, "dn_ds_summary_unequalpi_nobias.csv"))
dn.ds.sum.uneq.bias <- read.csv(paste0(DATADIR, "dn_ds_summary_unequalpi_bias.csv"))
dn.ds.sum.uneq <- rbind(dn.ds.sum.uneq.nobias, dn.ds.sum.uneq.bias)
dn.ds.sum.uneq$pi <- "unequal"
dn.ds.sum <- rbind(dn.ds.sum.eq, dn.ds.sum.uneq)
dn.ds.sum$method <- factor(dn.ds.sum$method, levels = method_order)
dn.ds.sum$pi <- factor(dn.ds.sum$pi, levels = c("equal", "unequal"))

dnds.sum.mean <- dnds.sum %>% group_by(method, ntaxa, bl, type, pi) %>% na.omit() %>% filter(!is.infinite(rmsd), !is.infinite(resvar)) %>% summarize(r = mean(r), estbias = mean(estbias), rmsd = mean(rmsd), resvar = mean(resvar))
dn.ds.sum.mean <- dn.ds.sum %>% group_by(method, ntaxa, bl, parameter, type, pi) %>% na.omit() %>% filter(!is.infinite(rmsd), !is.infinite(resvar)) %>% summarize(r = mean(r), estbias = mean(estbias), rmsd = mean(rmsd))
dnds.sum.mean$method <- factor(dnds.sum.mean$method, levels = method_order)
dn.ds.sum.mean$method <- factor(dn.ds.sum.mean$method, levels = method_order)
treelen.dat <- dnds.sum %>% mutate(treelen = 2*(ntaxa-1)*bl) %>% filter(method == "FEL1", treelen > 162, treelen < 164) %>% na.omit()

## Read in data of actual counted substitutions
unequal.count.bias <- read.csv(paste0(DATADIR, "substitution_counts_unequalpi_nobias.csv"))
unequal.count.nobias <- read.csv(paste0(DATADIR, "substitution_counts_unequalpi_bias.csv"))

### Read in and factorize linear model results
comp.order <- c("FUBAR1 - SLAC1", "FEL1 - SLAC1", "FEL1 - FUBAR1", "SLAC1 - SLAC2", "FUBAR1 - FUBAR2", "FEL1 - FEL2")
linmodels <- read.csv(paste0(DATADIR, "linear_model_results.csv")) # Note that this file was manually created generated using information found in the file "build_linear_models.R"
linmodels <- filter(linmodels, comp %in% comp.order)
linmodels$model <- factor(linmodels$model, levels = c("r", "rmsd", "resvar"))
linmodels$comp <- factor(linmodels$comp, levels = comp.order)
linmodels$type <- factor(linmodels$type, levels=c("nobias", "bias"), labels=c("No codon bias", "Codon bias"))
linmodels$pi   <- factor(linmodels$pi, levels=c("equal", "unequal"))

### Frequently used legends
reordered.dnds.sum.mean <- dnds.sum.mean 
reordered.dnds.sum.mean$method <- factor(reordered.dnds.sum.mean$method, levels=c("SLAC1", "SLAC2", "FEL1", "FEL2", "FUBAR1", "FUBAR2"), labels=c("SLAC1  ", "SLAC2  ", "FEL1  ", "FEL2  ", "FUBAR1", "FUBAR2"))
reordered.method.colors <- c("deeppink3", "skyblue", "red", "dodgerblue", "darkred", "blue")
methods.legend.grob <- reordered.dnds.sum.mean %>% filter(pi == "equal", type == "nobias") %>% ggplot(aes(x = factor(ntaxa), y = estbias, group = method, color = method)) + geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + scale_color_manual(values = reordered.method.colors, name = "Inference Method") + theme(legend.position = "bottom",legend.text = element_text(size=10), legend.title = element_text(size=11), legend.key.size = unit(.3, "cm"))
grobs <- ggplotGrob(methods.legend.grob)$grobs
methods.legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

reordered.dnds.sum.mean2 <- dnds.sum.mean %>% filter(method %in% c("FEL1", "FEL2", "FUBAR1", "FUBAR2"))
reordered.dnds.sum.mean2$method <- factor(reordered.dnds.sum.mean2$method, levels=c("FEL1", "FEL2", "FUBAR1", "FUBAR2"), labels=c("FEL1  ", "FEL2  ", "FUBAR1", "FUBAR2"))
reordered.method.colors2 <- c("red", "dodgerblue", "darkred", "blue")
methods.legend.grob2 <- reordered.dnds.sum.mean2 %>% filter(pi == "equal", type == "nobias") %>% ggplot(aes(x = factor(ntaxa), y = estbias, group = method, color = method)) + geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + scale_color_manual(values = reordered.method.colors2, name = "Inference Method") + theme(legend.position = "bottom",legend.text = element_text(size=10), legend.title = element_text(size=11), legend.key.size = unit(.3, "cm"))
grobs2 <- ggplotGrob(methods.legend.grob2)$grobs
methods.legend2 <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]


dn.ds.legend.grob <- dn.ds.sum %>% filter(pi == "equal", type == "bias", method == "FUBAR2") %>% ggplot(aes(x = factor(ntaxa), y = r, fill = parameter)) + geom_violin(scale = "width") + scale_fill_manual(values = dn_ds_colors, name = "Parameter", labels = c("dN    ", "dS")) + theme(legend.position = "bottom", legend.key.size = unit(.3, "cm"))
grobs <- ggplotGrob(dn.ds.legend.grob)$grobs
dn.ds.legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]



### Read in and process data for simulations performed along real trees
realdat <- read_csv(paste0(OUTDIR,"full_results_realtrees.csv"))
realdat.sum <- realdat %>% summarize_dnds_real()
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

### CORRELATION
r.nobias.equalpi <- dnds.sum.mean %>% filter(pi == "equal", type == "nobias") %>% ggplot(aes(x = factor(ntaxa), y = r, group = method, color=method)) +  geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 1), breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) + xlab("Number of Taxa") + ylab("Correlation") + scale_color_manual(values = method_colors, name = "Inference Method") + background_grid("xy") + ggtitle("No codon bias") + theme(axis.text.x = element_text(size = 9, angle=30)) 
r.bias.equalpi <- dnds.sum.mean %>% filter(pi == "equal", type == "bias") %>% ggplot(aes(x = factor(ntaxa), y = r, group = method, color=method)) +  geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 1), breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) + xlab("Number of Taxa") + ylab("Correlation") + scale_color_manual(values = method_colors, name = "Inference Method") + background_grid("xy") + ggtitle("Codon bias") + theme(axis.text.x = element_text(size = 9, angle=30)) 
corr.lineplots <- plot_grid(r.nobias.equalpi, r.bias.equalpi, nrow=2, labels=c("A", "B"))
save_plot(paste0(PLOTDIR, "correlation_lineplots_equalpi.pdf"), corr.lineplots, base_width = 10, base_height=5)

r.nobias.unequalpi <- dnds.sum.mean %>% filter(pi == "unequal", type == "nobias") %>% ggplot(aes(x = factor(ntaxa), y = r, group = method, color=method)) +  geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(-0.1, 1), breaks = c(-0.1,0.1, 0.3, 0.5, 0.7, 0.9)) + xlab("Number of Taxa") + ylab("Correlation") + scale_color_manual(values = method_colors, name = "Inference Method") + background_grid("xy") + ggtitle("No codon bias") + theme(axis.text.x = element_text(size = 9, angle=30)) 
r.bias.unequalpi <- dnds.sum.mean %>% filter(pi == "unequal", type == "bias") %>% ggplot(aes(x = factor(ntaxa), y = r, group = method, color=method)) +  geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(-0.1, 1), breaks = c(-0.1, 0.1, 0.3, 0.5, 0.7, 0.9)) + xlab("Number of Taxa") + ylab("Correlation") + scale_color_manual(values = method_colors, name = "Inference Method") + background_grid("xy") + ggtitle("Codon bias") + theme(axis.text.x = element_text(size = 9, angle=30)) 
corr.lineplots <- plot_grid(r.nobias.unequalpi, r.bias.unequalpi, nrow=2, labels=c("A", "B"))
save_plot(paste0(PLOTDIR, "correlation_lineplots_unequalpi.pdf"), corr.lineplots, base_width = 10, base_height=5)




## Estimator Bias
estbias.nobias.equalpi <- dnds.sum.mean %>% filter(pi == "equal", type == "nobias") %>% ggplot(aes(x = factor(ntaxa), y = estbias, group = method, color = method)) + geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(-0.5, 0.75)) + xlab("Number of Taxa") + ylab("Est. Bias") + scale_color_manual(values = method_colors, name = "Inference Method") + background_grid("xy") + ggtitle("No codon bias")  + theme(axis.text.x = element_text(size = 8, angle=30), legend.position = "none") 
estbias.bias.equalpi <- dnds.sum.mean %>% filter(pi == "equal", type == "bias") %>% ggplot(aes(x = factor(ntaxa), y = estbias, group = method, color = method)) + geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(-0.5, 0.75)) + xlab("Number of Taxa") + ylab("Est. Bias") + scale_color_manual(values = method_colors, name = "Inference Method") + background_grid("xy") + ggtitle("Codon bias") + theme(axis.text.x = element_text(size = 8, angle=30), legend.position = "none") 
estbias.nobias.unequalpi <- dnds.sum.mean %>% filter(pi == "unequal", type == "nobias") %>% ggplot(aes(x = factor(ntaxa), y = estbias, group = method, color = method)) + geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(-0.5, 0.75)) + xlab("Number of Taxa") + ylab("Est. Bias") + scale_color_manual(values = method_colors, name = "Inference Method") + background_grid("xy") + ggtitle("No codon bias")  + theme(axis.text.x = element_text(size = 8, angle=30), legend.position = "none") 
estbias.bias.unequalpi <- dnds.sum.mean %>% filter(pi == "unequal", type == "bias") %>% ggplot(aes(x = factor(ntaxa), y = estbias, group = method, color = method)) + geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(-0.5, 0.75)) + xlab("Number of Taxa") + ylab("Est. Bias") + scale_color_manual(values = method_colors, name = "Inference Method") + background_grid("xy") + ggtitle("Codon bias")  + theme(axis.text.x = element_text(size = 8, angle=30), legend.position = "none") 
estbias.grid <- plot_grid(estbias.nobias.equalpi, estbias.bias.equalpi, estbias.nobias.unequalpi, estbias.bias.unequalpi, methods.legend,  nrow=5, labels=c("A", "B", "C", "D"), rel_heights=c(1,1,1,1,0.17)) + theme(axis.text.x = element_text(size = 9, angle=30)) 
save_plot(paste0(PLOTDIR, "estbias_lineplots.pdf"), estbias.grid, base_width = 7, base_height=9)




## RMSD
rmsd.nobias.equalpi <- dnds.sum.mean %>% filter(pi == "equal", type == "nobias", bl >=0.01) %>% ggplot(aes(x = factor(ntaxa), y = rmsd, group = method, color = method)) + geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 1.5)) + xlab("Number of Taxa") + ylab("RMSD") + scale_color_manual(values = method_colors, name = "Inference Method") + background_grid("xy") + ggtitle("No codon bias") + theme(axis.text.x = element_text(size = 9, angle=30), legend.position = "none")
rmsd.bias.equalpi <- dnds.sum.mean %>% filter(pi == "equal", type == "bias", bl >=0.01) %>% ggplot(aes(x = factor(ntaxa), y = rmsd, group = method, color = method)) + geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 1.5)) + xlab("Number of Taxa") + ylab("RMSD") + scale_color_manual(values = method_colors, name = "Inference Method") + background_grid("xy") + ggtitle("Codon bias")+ theme(axis.text.x = element_text(size = 9, angle=30), legend.position = "none")
rmsd.nobias.unequalpi <- dnds.sum.mean %>% filter(pi == "unequal", type == "nobias", bl >=0.01) %>% ggplot(aes(x = factor(ntaxa), y = rmsd, group = method, color = method)) + geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 1.5)) + xlab("Number of Taxa") + ylab("RMSD") + scale_color_manual(values = method_colors, name = "Inference Method") + background_grid("xy") + ggtitle("No codon bias")+ theme(axis.text.x = element_text(size = 9, angle=30), legend.position = "none")
rmsd.bias.unequalpi <- dnds.sum.mean %>% filter(pi == "unequal", type == "bias", bl >=0.01) %>% ggplot(aes(x = factor(ntaxa), y = rmsd, group = method, color = method)) + geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 1.5)) + xlab("Number of Taxa") + ylab("RMSD") + scale_color_manual(values = method_colors, name = "Inference Method") + background_grid("xy") + ggtitle("Codon bias")+ theme(axis.text.x = element_text(size = 9, angle=30), legend.position = "none")
rmsd.lineplots <- plot_grid(rmsd.nobias.equalpi, rmsd.bias.equalpi, rmsd.nobias.unequalpi, rmsd.bias.unequalpi,methods.legend2, nrow=5, labels=c("A", "B", "C", "D"), rel_heights=c(1,1,1,1,0.17))
save_plot(paste0(PLOTDIR, "rmsd_lineplots.pdf"), rmsd.lineplots, base_width = 7, base_height=8.5)


## Residual variance
resvar.nobias.equalpi <- dnds.sum.mean %>% filter(pi == "equal", type == "nobias", bl >=0.01) %>% ggplot(aes(x = factor(ntaxa), y = resvar, group = method, color = method)) + geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 2.5)) + xlab("Number of Taxa") + ylab("Res. Variance") + scale_color_manual(values = method_colors, name = "Inference Method") + background_grid("xy") + ggtitle("No codon bias") + theme(axis.text.x = element_text(size = 9, angle=30), legend.position = "none")
resvar.bias.equalpi <- dnds.sum.mean %>% filter(pi == "equal", type == "bias", bl >=0.01) %>% ggplot(aes(x = factor(ntaxa), y = resvar, group = method, color = method)) + geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 2.5)) + xlab("Number of Taxa") + ylab("Res. Variance") + scale_color_manual(values = method_colors, name = "Inference Method") + background_grid("xy") + ggtitle("Codon bias") + theme(axis.text.x = element_text(size = 9, angle=30), legend.position = "none")
resvar.nobias.unequalpi <- dnds.sum.mean %>% filter(pi == "unequal", type == "nobias", bl >=0.01) %>% ggplot(aes(x = factor(ntaxa), y = resvar, group = method, color = method)) + geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 2)) + xlab("Number of Taxa") + ylab("Res. Variance") + scale_color_manual(values = method_colors, name = "Inference Method") + background_grid("xy") + ggtitle("No codon bias") + theme(axis.text.x = element_text(size = 9, angle=30), legend.position = "none")
resvar.bias.unequalpi <- dnds.sum.mean %>% filter(pi == "unequal", type == "bias", bl >=0.01) %>% ggplot(aes(x = factor(ntaxa), y = resvar, group = method, color = method)) + geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 2)) + xlab("Number of Taxa") + ylab("Res. Variance") + scale_color_manual(values = method_colors, name = "Inference Method") + background_grid("xy") + ggtitle("Codon bias") + theme(axis.text.x = element_text(size = 9, angle=30), legend.position = "none")

resvar.lineplots <- plot_grid(resvar.nobias.equalpi, resvar.bias.equalpi, resvar.nobias.unequalpi, resvar.bias.unequalpi, nrow=5, methods.legend, labels=c("A", "B", "C", "D"), rel_heights=c(1,1,1,1,0.15))
save_plot(paste0(PLOTDIR, "resvar_lineplots.pdf"), resvar.lineplots, base_width = 7, base_height=8.5)



####################################################################################################
######### Violin plots  of dN & dS correlations and RMSD for balanced tree simulations #############
####################################################################################################

param.r.fel2.sub    <- dn.ds.sum %>% filter(pi == "unequal", type == "bias", method == "FEL2", bl >= 0.01, ntaxa<=1024) %>% ggplot(aes(x = factor(ntaxa), y = r, fill = parameter)) + geom_violin(scale = "width") + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(-0.1, 1), breaks = c(-0.1, 0.1, 0.3, 0.5, 0.7, 0.9)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("Correlation")  + scale_fill_manual(values = dn_ds_colors) + theme(legend.position = "none")
#param.estbias.fel2.sub <- dn.ds.sum %>% filter(pi == "unequal", type == "bias", method == "FEL2", bl >= 0.01, ntaxa<=1024) %>% ggplot(aes(x = factor(ntaxa), y = estbias, fill = parameter)) + geom_violin(scale = "width") + facet_grid(~bl, scales = "free_x") + coord_cartesian(ylim=c(0,2)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("Estimator Bias") + scale_fill_manual(values = dn_ds_colors) + theme(legend.position = "none")
param.rmsd.fel2.sub <- dn.ds.sum %>% filter(pi == "unequal", type == "bias", method == "FEL2", bl >= 0.01, ntaxa<=1024, rmsd <= 5) %>% ggplot(aes(x = factor(ntaxa), y = rmsd, fill = parameter)) + geom_violin(scale = "width") + facet_grid(~bl, scales = "free_x") + coord_cartesian(ylim=c(0,3)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("RMSD") + scale_fill_manual(values = dn_ds_colors) + theme(legend.position = "none")

g.corr <- ggplotGrob(param.r.fel2.sub)
g.rmsd <- ggplotGrob(param.rmsd.fel2.sub)
g.corr$widths <- grid:::unit.list(g.corr$widths)   
g.rmsd$widths <-  grid:::unit.list(g.rmsd$widths)
corr.widths <- g.corr$widths[1:3]
rmsd.widths <- g.rmsd$widths[1:3] 
max.widths <- unit.pmax(corr.widths, rmsd.widths) # calculate maximum widths
g.corr$widths[1:3] <- max.widths
g.rmsd$widths[1:3] <- max.widths


dn.ds.fel2.sub.legend <- plot_grid(g.corr, g.rmsd, dn.ds.legend, nrow=3, labels = c("A", "B"), rel_heights=c(1,1,0.1))
save_plot(paste0(PLOTDIR, "dn_ds_unequal_violin_fel2_subset.pdf"), dn.ds.fel2.sub.legend, base_width = 6.75, base_height=5.5)



### No more fubar2 here, since only FEL2 analysis is performed. ###
# param.r.fubar2.equal <- dn.ds.sum %>% filter(pi == "equal", type == "bias", method == "FUBAR2") %>% ggplot(aes(x = factor(ntaxa), y = r, fill = parameter)) + geom_violin(scale = "width") + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(-0.1, 1), breaks = c(-0.1, 0.1, 0.3, 0.5, 0.7, 0.9)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("Correlation") + scale_fill_manual(values = dn_ds_colors) + theme(legend.position = "none")
# param.r.fubar2.unequal <- dn.ds.sum %>% filter(pi == "unequal", type == "bias", method == "FUBAR2") %>% ggplot(aes(x = factor(ntaxa), y = r, fill = parameter)) + geom_violin(scale = "width") + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(-0.1, 1), breaks = c(-0.1, 0.1, 0.3, 0.5, 0.7, 0.9)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("Correlation") + scale_fill_manual(values = dn_ds_colors) + theme(legend.position = "none")
# param.estbias.fubar2.equal <- dn.ds.sum %>% filter(pi == "equal", type == "bias", method == "FUBAR2") %>% ggplot(aes(x = factor(ntaxa), y = estbias, fill = parameter)) + geom_violin(scale = "width") + facet_grid(~bl, scales = "free_x") + coord_cartesian(ylim=c(0,2)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("Estimator Bias") + scale_fill_manual(values = dn_ds_colors)  + theme(legend.position = "none")
# #param.estbias.fubar2.unequal <- dn.ds.sum %>% filter(pi == "unequal", type == "bias", method == "FUBAR2") %>% ggplot(aes(x = factor(ntaxa), y = estbias, fill = parameter)) + geom_violin(scale = "width") + facet_grid(~bl, scales = "free_x") + coord_cartesian(ylim=c(0,2)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("Estimator Bias") + scale_fill_manual(values = dn_ds_colors) + theme(legend.position = "none")
# param.rmsd.fubar2.equal <- dn.ds.sum %>% filter(pi == "equal", type == "bias", method == "FUBAR2") %>% ggplot(aes(x = factor(ntaxa), y = rmsd, fill = parameter)) + geom_violin(scale = "width") + facet_grid(~bl, scales = "free_x") + coord_cartesian(ylim=c(0,5)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("RMSD") + scale_fill_manual(values = dn_ds_colors)  + theme(legend.position = "none")
# param.rmsd.fubar2.unequal <- dn.ds.sum %>% filter(pi == "unequal", type == "bias", method == "FUBAR2") %>% ggplot(aes(x = factor(ntaxa), y = rmsd, fill = parameter)) + geom_violin(scale = "width") + facet_grid(~bl, scales = "free_x") + coord_cartesian(ylim=c(0,5)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("RMSD") + scale_fill_manual(values = dn_ds_colors) + theme(legend.position = "none")
# 
# 
# 
# param.r.fubar2 <- plot_grid(param.r.fubar2.equal, param.r.fubar2.unequal, nrow=1, labels=c("A", "B"), scale=0.98)
# param.estbias.fubar2 <- plot_grid(param.estbias.fubar2.equal, param.estbias.fubar2.unequal, nrow=1, labels=c("C", "D"))
# param.rmsd.fubar2 <- plot_grid(param.rmsd.fubar2.equal, param.rmsd.fubar2.unequal, nrow=1, labels=c("C", "D"),scale=0.98)
# param.r.rmsd.fubar2 <- plot_grid(param.r.fubar2, param.rmsd.fubar2, dn.ds.legend, nrow=3, rel_heights=c(1,1,0.1))
# save_plot(paste0(PLOTDIR, "dn_ds_r_rmsd_fubar2.pdf"), param.r.rmsd.fubar2, base_width = 11.5, base_height=5.5)
# 



param.r.fel2.equal <- dn.ds.sum %>% filter(pi == "equal", type == "bias", method == "FEL2", bl>=0.01) %>% ggplot(aes(x = factor(ntaxa), y = r, fill = parameter)) + geom_violin(scale = "width") + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(-0.1, 1), breaks = c(-0.1, 0.1, 0.3, 0.5, 0.7, 0.9)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("Correlation") + scale_fill_manual(values = dn_ds_colors)  + theme(axis.text.x = element_text(size = 8, angle=30), legend.position = "none") 
param.r.fel2.unequal <- dn.ds.sum %>% filter(pi == "unequal", type == "bias", method == "FEL2", bl>=0.01) %>% ggplot(aes(x = factor(ntaxa), y = r, fill = parameter)) + geom_violin(scale = "width") + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(-0.1, 1), breaks = c(-0.1, 0.1, 0.3, 0.5, 0.7, 0.9)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("Correlation") + scale_fill_manual(values = dn_ds_colors)  + theme(axis.text.x = element_text(size = 8, angle=30), legend.position = "none") 
#param.estbias.fel2.equal <- dn.ds.sum %>% filter(pi == "equal", type == "bias", method == "FEL2", bl>=0.01) %>% ggplot(aes(x = factor(ntaxa), y = estbias, fill = parameter)) + geom_violin(scale = "width") + facet_grid(~bl, scales = "free_x") + coord_cartesian(ylim=c(0,2)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("Estimator Bias") + scale_fill_manual(values = dn_ds_colors)  + theme(legend.position = "none")
#param.estbias.fel2.unequal <- dn.ds.sum %>% filter(pi == "unequal", type == "bias", method == "FEL2", bl>=0.01) %>% ggplot(aes(x = factor(ntaxa), y = estbias, fill = parameter)) + geom_violin(scale = "width") + facet_grid(~bl, scales = "free_x") + coord_cartesian(ylim=c(0,2)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("Estimator Bias") + scale_fill_manual(values = dn_ds_colors) + theme(legend.position = "none")
param.rmsd.fel2.equal <- dn.ds.sum %>% filter(pi == "equal", type == "bias", method == "FEL2", bl>=0.01) %>% ggplot(aes(x = factor(ntaxa), y = rmsd, fill = parameter)) + geom_violin(scale = "width") + facet_grid(~bl, scales = "free_x") + coord_cartesian(ylim=c(0,5)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("RMSD") + scale_fill_manual(values = dn_ds_colors)  + theme(axis.text.x = element_text(size = 8, angle=30), legend.position = "none") 
param.rmsd.fel2.unequal <- dn.ds.sum %>% filter(pi == "unequal", type == "bias", method == "FEL2", bl>=0.01) %>% ggplot(aes(x = factor(ntaxa), y = rmsd, fill = parameter)) + geom_violin(scale = "width") + facet_grid(~bl, scales = "free_x") + coord_cartesian(ylim=c(0,5)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("RMSD") + scale_fill_manual(values = dn_ds_colors)  + theme(axis.text.x = element_text(size = 8, angle=30), legend.position = "none") 

param.r.fel2 <- plot_grid(param.r.fel2.equal, param.r.fel2.unequal, nrow=1, labels=c("A", "B"), scale=0.98)
#param.estbias.fel2 <- plot_grid(param.estbias.fel2.equal, param.estbias.fel2.unequal, nrow=1, labels=c("C", "D"))
param.rmsd.fel2 <- plot_grid(param.rmsd.fel2.equal, param.rmsd.fel2.unequal, nrow=1, labels=c("C", "D"),scale=0.98)
param.r.rmsd.fel2 <- plot_grid(param.r.fel2, param.rmsd.fel2, dn.ds.legend, nrow=3, rel_heights=c(1,1,0.1))
save_plot(paste0(PLOTDIR, "dn_ds_r_rmsd_fel2.pdf"), param.r.rmsd.fel2, base_width = 10, base_height=5)

####################################################################################################
###################################### Linear models plot ##########################################
####################################################################################################
sig_colors <- c("grey60", "black")
ptsize = 1.5
linesize = 0.8

linmodels %>% filter(model == "r", pi == "unequal") %>% arrange(comp) %>%
  ggplot(aes(x = coeff, y = comp, color = sig)) + geom_point(size=ptsize) + 
  geom_segment(aes(x=lowerCI,xend=upperCI,y=comp,yend=comp),size=linesize) + 
  geom_vline(xintercept=0, size=0.5) + 
  facet_grid(~type) + background_grid() +
  xlab("Average Correlation Difference") + ylab("Comparison") + 
  scale_x_continuous(limits=c(-0.025, 0.15)) +  
  scale_color_manual(values = sig_colors)  + theme(legend.position = "none", axis.text.y = element_text(size=10)) -> linmodel.unequal.r

linmodels %>% filter(model == "rmsd", pi == "unequal") %>% arrange(comp) %>%
  ggplot(aes(x = coeff, y = comp, color = sig)) + geom_point(size=ptsize) + 
  geom_segment(aes(x=lowerCI,xend=upperCI,y=comp,yend=comp),size=linesize) + 
  geom_vline(xintercept=0, size=0.5) + 
  facet_grid(~type) + background_grid() +
  xlab("Average RMSD Difference") + ylab("Comparison") + 
  scale_x_continuous(limits=c(-0.25,0.1)) +  
  scale_color_manual(values = sig_colors) + theme(legend.position = "none", axis.text.y = element_text(size=10)) -> linmodel.unequal.rmsd

linmodel.unequal <- plot_grid(linmodel.unequal.r, linmodel.unequal.rmsd, nrow=2, labels=c("A", "B"))
save_plot(paste0(PLOTDIR, "linmodels_unequal.pdf"), linmodel.unequal, base_width = 6.5, base_height = 5)


linmodels %>% filter(model == "r", pi == "equal") %>% arrange(comp) %>%
    ggplot(aes(x = coeff, y = comp, color = sig)) + geom_point(size=ptsize) + 
    geom_segment(aes(x=lowerCI,xend=upperCI,y=comp,yend=comp),size=linesize) + 
    geom_vline(xintercept=0, size=0.5) + 
    facet_grid(~type) + background_grid() +
    xlab("Average Correlation Difference") + ylab("Comparison") + 
    scale_x_continuous(limits=c(-0.025, 0.15)) +  
    scale_color_manual(values = sig_colors)  + theme(legend.position = "none", axis.text.y = element_text(size=10))  -> linmodel.equal.r

linmodels %>% filter(model == "rmsd", pi == "equal") %>% arrange(comp) %>%
    ggplot(aes(x = coeff, y = comp, color = sig)) + geom_point(size=ptsize) + 
    geom_segment(aes(x=lowerCI,xend=upperCI,y=comp,yend=comp),size=linesize) + 
    geom_vline(xintercept=0, size=0.5) + 
    facet_grid(~type) + background_grid() +
    xlab("Average RMSD Difference") + ylab("Comparison") + 
    scale_x_continuous(limits=c(-0.25,0.1)) +  
    scale_color_manual(values = sig_colors)  + theme(legend.position = "none", axis.text.y = element_text(size=10))   -> linmodel.equal.rmsd

linmodel.equal <- plot_grid(linmodel.equal.r, linmodel.equal.rmsd, nrow=1, labels=c("A", "B"), scale=0.98)
save_plot(paste0(PLOTDIR, "linmodels_equal.pdf"), linmodel.equal, base_width = 12, base_height = 3)




####################################################################################################
###################################### Counted n/s plots ###########################################
####################################################################################################


unequal.count.nobias %>% group_by(ntaxa, bl, rep) %>% mutate(ratio_ns = ncount/scount) %>% na.omit() %>% filter(!is.infinite(ratio_ns)) %>% summarize(mean_ratio = mean(ratio_ns))  -> ns_count_ratio
count.grob1 <- ggplot(ns_count_ratio, aes(x = as.factor(bl), y = mean_ratio, fill = as.factor(ntaxa))) + geom_boxplot(outlier.size = 0.75) +  scale_fill_brewer(palette = "YlGn", name = "Number of Taxa   ") + theme(legend.position = "bottom") 
count.grob <- ggplotGrob(count.grob1)$grobs
count.legend <- count.grob[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

unequal.nobias <- ggplot(ns_count_ratio, aes(x = as.factor(bl), y = mean_ratio, fill = as.factor(ntaxa))) + geom_boxplot(outlier.size = 0.75) + geom_hline(yintercept=1) + scale_fill_brewer(palette = "YlGn", name = "Number of Taxa") + xlab("Branch Length") + ylab("Mean N/S count ratio") + scale_y_continuous(limits=c(0,4)) + background_grid("xy") + ggtitle("No codon bias") + theme(legend.position = "none") 
unequal.count.bias %>% group_by(ntaxa, bl, rep) %>% mutate(ratio_ns = ncount/scount) %>% na.omit() %>% filter(!is.infinite(ratio_ns)) %>% summarize(mean_ratio = mean(ratio_ns))  -> ns_count_ratio
unequal.bias <- ggplot(ns_count_ratio, aes(x = as.factor(bl), y = mean_ratio, fill = as.factor(ntaxa))) + geom_boxplot(outlier.size = 0.75) + geom_hline(yintercept=1) + scale_fill_brewer(palette = "YlGn", name = "Number of Taxa") + xlab("Branch Length") + ylab("Mean N/S count ratio") +  scale_y_continuous(limits=c(0,4)) + background_grid("xy") + ggtitle("Codon bias") + theme(legend.position = "none") 
unequal.count <- plot_grid(unequal.nobias, unequal.bias, nrow=1, labels=c("A", "B"), scale=0.95)
unequal.count.withlegend <- plot_grid(unequal.count, count.legend, nrow=2, rel_heights=c(1,0.1))

save_plot(paste0(PLOTDIR, "counted_ratio_unequal.pdf"), unequal.count.withlegend, base_width = 7.5, base_height = 3.25)


####################################################################################################
######################### Correlations and RMSD for real tree simulations ##########################
####################################################################################################

fubar.legend.grob <- realdat.sum %>% ggplot(aes(x = dataset, fill = method, y = r)) + geom_boxplot(position=position_dodge(0.77), outlier.size = 0.75, size=0.4) + facet_wrap(~type) + scale_fill_manual(values=dn_ds_colors, name = "Inference Method   ", labels=c("FUBAR1 ", "FUBAR2")) + theme(legend.position = "bottom")
grobs <- ggplotGrob(fubar.legend.grob)$grobs
fubar.legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

r.dnds <- realdat.sum %>% ggplot(aes(x = dataset, fill = method, y = r)) + geom_boxplot(position=position_dodge(0.77), outlier.size = 0.75, size=0.4) + facet_wrap(~type) + scale_fill_manual(values=dn_ds_colors, name = "Inference Method") +  xlab("Dataset") + ylab("Correlation") + scale_y_continuous(limits=c(0,1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.)) + theme(legend.position = "none")
rmsd.dnds <- realdat.sum %>% ggplot(aes(x = dataset, fill = method, y = rmsd)) + geom_boxplot(position=position_dodge(0.77), outlier.size = 0.75, size=0.4) + facet_wrap(~type) + scale_fill_manual(values=dn_ds_colors, name = "Inference Method") + xlab("Dataset") + ylab("RMSD") + scale_y_continuous(limits=c(0,2), breaks=c(0, 0.5, 1, 1.5, 2)) + theme(legend.position = "none")
#resvar.dnds <- realdat.sum %>% ggplot(aes(x = dataset, fill = method, y = resvar)) + geom_boxplot(position=position_dodge(0.77), outlier.size = 0.75, size=0.4) + facet_wrap(~type) + scale_fill_manual(values=dn_ds_colors, name = "Inference Method") + xlab("Dataset") + ylab("Residual Variance") + scale_y_continuous(limits=c(0,2.5), breaks=c(0, 0.5, 1, 1.5, 2, 2.5)) + theme(legend.position = "none")

real.plots <- plot_grid(r.dnds, rmsd.dnds, nrow=2, labels=c("A", "B"))
real.plots.withlegend <- plot_grid(real.plots, fubar.legend, nrow=2, rel_heights=c(1,0.08) )
save_plot(paste0(PLOTDIR, "real_r_rmsd.pdf"), real.plots.withlegend, base_width = 7.25, base_height = 4.75)

r.dn.ds <- dn.ds.real %>% ggplot(aes(x = dataset, fill = parameter, y = r)) + geom_boxplot(outlier.size = 0.75, size=0.4) + scale_fill_manual(values = dn_ds_colors, name = "Parameter") + ylab("Correlation") + xlab("Dataset") + theme(legend.position = "none")
rmsd.dn.ds <- dn.ds.real %>% ggplot(aes(x = dataset, fill = parameter, y = rmsd)) + geom_boxplot(outlier.size = 0.75, size=0.4) + scale_fill_manual(values = dn_ds_colors, name = "Parameter") + ylab("RMSD") + xlab("Dataset") + theme(legend.position = "none")
real.plots.dn.ds <- plot_grid(r.dn.ds, rmsd.dn.ds, nrow=1, labels=c("A", "B"))
real.plots.dn.ds.legend <- plot_grid(real.plots.dn.ds, dn.ds.legend, nrow=2, rel_heights=c(1,0.1) )

save_plot(paste0(PLOTDIR, "real_dn_ds_r_rmsd.pdf"), real.plots.dn.ds.legend, base_width = 9, base_height = 2.5)


####################################################################################################
########################### Scatterplots of dN vs. dN for FEL and FUBAR ############################
####################################################################################################

full.nobias <- read_csv(paste0(OUTDIR,"full_results_unequalpi_nobias.csv"))
full.bias <- read_csv(paste0(OUTDIR,"full_results_unequalpi_bias.csv"))
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

nobias.fel.dn <-nobias.compare.dn %>% ggplot(aes(x = FEL1, y = FEL2)) + geom_point() + facet_grid(ntaxa~bl) + geom_abline(slope=1, intercept=0, color="red", size=0.7) + scale_x_continuous(limits=c(0,1), breaks=c(0, 0.5, 1)) + scale_y_continuous(limits=c(0,1), breaks=c(0, 0.5, 1)) + xlab("dN, FEL1") + ylab("dN, FEL2") + ggtitle("No codon bias")
bias.fel.dn <- bias.compare.dn %>% ggplot(aes(x = FEL1, y = FEL2)) + geom_point() + facet_grid(ntaxa~bl) + geom_abline(slope=1, intercept=0, color="red", size=0.7) + scale_x_continuous(limits=c(0,1), breaks=c(0, 0.5, 1)) + scale_y_continuous(limits=c(0,1), breaks=c(0, 0.5, 1)) + xlab("dN, FEL1") + ylab("dN, FEL2") + ggtitle("Codon bias")

dn.same.grid <- plot_grid(nobias.fel.dn, bias.fel.dn, nrow=1, labels=c("A", "B"), label_size=20)
save_plot(paste0(PLOTDIR, "samedn.pdf"), dn.same.grid, base_width = 10, base_height=5)




####################################################################################################
############################## Divergence vs. number of taxa #######################################
####################################################################################################


legend.grob <- treelen.dat %>% ggplot(aes(x = as.factor(bl), y = r, fill = as.factor(type))) + geom_boxplot(position = position_dodge(0.9)) + xlab("Parameterization") + scale_fill_manual(name = "", labels = c("No Codon bias     ", "Codon bias"), values = c("grey40", "grey80")) + theme(legend.position = "bottom")
grobs <- ggplotGrob(legend.grob)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

ntaxa.bl.r <- treelen.dat %>% filter(pi == "unequal") %>% ggplot(aes(x = as.factor(bl), y = r, fill = as.factor(type))) + geom_boxplot(position = position_dodge(0.75), lwd=0.3, outlier.size=0.8) + xlab("Parameterization") + scale_x_discrete(labels = c("N = 2048\nB = 0.04","N = 512\nB = 0.16","N = 128\nB = 0.64")) + ylab("Correlation") + scale_y_continuous(limits=c(0.45, 0.95)) + scale_fill_manual(values = c("grey40", "grey80")) + theme(legend.position = "none")
ntaxa.bl.rmsd <- treelen.dat %>% filter(pi == "unequal") %>% ggplot(aes(x = as.factor(bl), y = rmsd, fill = as.factor(type))) + geom_boxplot(position = position_dodge(0.75), lwd=0.3, outlier.size=0.8) + xlab("Parameterization") + scale_x_discrete(labels = c("N = 2048\nB = 0.04","N = 512\nB = 0.16","N = 128\nB = 0.64")) + ylab("RMSD") + scale_y_continuous(limits=c(0,0.75)) + scale_fill_manual(values = c("grey40", "grey80")) + theme(legend.position = "none")
ntaxa.bl.r.rmsd <- plot_grid(ntaxa.bl.r, ntaxa.bl.rmsd, nrow=1, labels=c("A", "B"))
ntaxa.bl.r.rmsd.legend <- plot_grid(ntaxa.bl.r.rmsd, legend, nrow=2, rel_heights=c(1,0.1))
save_plot(paste0(PLOTDIR, "ntaxa_bl_r_rmsd_unequalpi.pdf"), ntaxa.bl.r.rmsd.legend, base_width=8, base_height=3.25)

ntaxa.bl.r <- treelen.dat %>% filter(pi == "equal") %>% ggplot(aes(x = as.factor(bl), y = r, fill = as.factor(type))) + geom_boxplot(position = position_dodge(0.75), lwd=0.3, outlier.size=0.8) + xlab("Parameterization") + scale_x_discrete(labels = c("N = 2048\nB = 0.04","N = 512\nB = 0.16","N = 128\nB = 0.64")) + ylab("Correlation") + scale_y_continuous(limits=c(0.45, 0.95)) + scale_fill_manual(values = c("grey40", "grey80")) + theme(legend.position = "none")
ntaxa.bl.rmsd <- treelen.dat %>% filter(pi == "equal") %>% ggplot(aes(x = as.factor(bl), y = rmsd, fill = as.factor(type))) + geom_boxplot(position = position_dodge(0.75), lwd=0.3, outlier.size=0.8) + xlab("Parameterization") + scale_x_discrete(labels = c("N = 2048\nB = 0.04","N = 512\nB = 0.16","N = 128\nB = 0.64")) + ylab("RMSD") + scale_y_continuous(limits=c(0,0.75)) + scale_fill_manual(values = c("grey40", "grey80")) + theme(legend.position = "none")
ntaxa.bl.r.rmsd <- plot_grid(ntaxa.bl.r, ntaxa.bl.rmsd, nrow=1, labels=c("A", "B"))
ntaxa.bl.r.rmsd.legend <- plot_grid(ntaxa.bl.r.rmsd, legend, nrow=2, rel_heights=c(1,0.1))
save_plot(paste0(PLOTDIR, "ntaxa_bl_r_rmsd_equalpi.pdf"), ntaxa.bl.r.rmsd.legend, base_width=8, base_height=3.25)
