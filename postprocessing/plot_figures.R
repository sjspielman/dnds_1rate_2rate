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
dnds.sum <- read.csv("dnds_summary.csv")
dn.ds.sum <- read.csv("dn_ds_summary.csv")
dnds.sum$method <- factor(dnds.sum$method, levels = method_order)
dn.ds.sum$method <- factor(dn.ds.sum$method, levels = method_order)

dnds.sum.mean <- dnds.sum %>% group_by(method, ntaxa, bl, type) %>% summarize(r = mean(r), estbias = mean(estbias), rmsd = mean(rmsd))
dn.ds.sum.mean <- dn.ds.sum %>% group_by(method, ntaxa, bl, parameter, type) %>% summarize(r = mean(r), estbias = mean(estbias), rmsd = mean(rmsd))
dnds.sum.mean$method <- factor(dnds.sum.mean$method, levels = method_order)
dn.ds.sum.mean$method <- factor(dn.ds.sum.mean$method, levels = method_order)
treelen.dat <- dnds.sum %>% mutate(treelen = 2*(ntaxa-1)*bl) %>% filter(method == "FEL1", treelen > 162, treelen < 164) %>% na.omit() %>% select(-estbias) %>% group_by(treelen, type, ntaxa, bl, rep)

## Read in data of actual counted substitutions
gtr.count <- read.csv("substitution_counts_gtr_nobias.csv")
bias.gtr.count <- read.csv("substitution_counts_gtr_bias.csv")

### Read in and process linear model results
comp.order <- c("FUBAR1 - SLAC1", "FEL1 - SLAC1", "FEL1 - FUBAR1", "SLAC1 - SLAC2", "FUBAR1 - FUBAR2", "FEL1 - FEL2")
linmodels.correlations <- read.csv("linear_model_results_correlation.csv") # Note that this file was manually created generated using information found in the file "build_linear_models.R"
linmodels.correlations$modeltype <- "r"
linmodels.rmsd <- read.csv("linear_model_results_rmsd.csv") # Note that this file was manually created generated using information found in the file "build_linear_models.R"
linmodels.rmsd$modeltype <- "rmsd"
linmodels.all <-rbind(linmodels.correlations, linmodels.rmsd) %>% filter(comp %in% comp.order)
linmodels.all$modeltype <- factor(linmodels.all$modeltype, levels = c("r", "rmsd"))
linmodels.all$comp <- factor(linmodels.all$comp, levels = comp.order)
linmodels.all$type <- factor(linmodels.all$type, levels=c("nobias", "bias"), labels=c("No codon bias", "Codon bias"))





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

r.nobias <- dnds.sum.mean %>% filter(type == "nobias") %>% ggplot(aes(x = factor(ntaxa), y = r, group = method, color=method)) +  geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 1), breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) + xlab("Number of Taxa") + ylab("Correlation") + scale_color_manual(values = method_colors, name = "Inference Method") + background_grid("xy") + ggtitle("No codon bias")
r.bias <- dnds.sum.mean %>% filter(type == "bias") %>% ggplot(aes(x = factor(ntaxa), y = r, group = method, color=method)) +  geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 1), breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) + xlab("Number of Taxa") + ylab("Correlation") + scale_color_manual(values = method_colors, name = "Inference Method") + background_grid("xy") + ggtitle("Codon bias")
corr.lineplots <- plot_grid(r.nobias, r.bias, nrow=2, labels=c("A", "B"))
save_plot(paste0(PLOTDIR, "correlation_lineplots.pdf"), corr.lineplots, base_width = 10, base_height=5)

theme_set(theme_cowplot() + theme(axis.text = element_text(size = 10), 
                                  axis.title = element_text(size = 11), 
                                  plot.title = element_text(size=11, face="bold"),
                                  panel.border = element_rect(size = 0.5), 
                                  panel.margin = unit(0.75, "lines"), 
                                  strip.text = element_text(size=11), 
                                  legend.text = element_text(size=10), 
                                  legend.title = element_text(size=11)))

rmsd.nobias <- dnds.sum %>% filter(type == "nobias", bl >=0.04, ntaxa >= 256, method %in% c("FEL1", "FEL2", "FUBAR1", "FUBAR2")) %>% ggplot(aes(x = factor(ntaxa), y = rmsd, fill=method)) +  geom_violin(scale="width", lwd=0.2) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 1.25)) + xlab("Number of Taxa") + ylab("RMSD") + scale_fill_manual(values = felfubar_colors, name = "Inference Method") + background_grid("xy") + ggtitle("No codon bias")
rmsd.bias <- dnds.sum %>% filter(type == "bias", bl >=0.04, ntaxa >= 256, method %in% c("FEL1", "FEL2", "FUBAR1", "FUBAR2")) %>% ggplot(aes(x = factor(ntaxa), y = rmsd, fill=method)) +  geom_violin(scale="width", lwd=0.2) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 1.25)) + xlab("Number of Taxa") + ylab("RMSD") + scale_fill_manual(values = felfubar_colors, name = "Inference Method") + background_grid("xy") + ggtitle("Codon bias")

rmsd.boxplots <- plot_grid(rmsd.nobias, rmsd.bias, nrow=2, labels=c("A", "B"))
save_plot(paste0(PLOTDIR, "rmsd_violins.pdf"), rmsd.boxplots, base_width = 9, base_height=5)





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

dn.r.fel2    <- dn.ds.sum %>% filter(type == "bias", method == "FEL2", bl >= 0.01, ntaxa<=1024) %>% ggplot(aes(x = factor(ntaxa), y = r, fill = parameter)) + geom_violin(scale = "width") + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(-0.1, 1), breaks = c(-0.1, 0.1, 0.3, 0.5, 0.7, 0.9)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("Correlation") + scale_fill_manual(values = dn_ds_colors, name = "Parameter", labels = c("dN", "dS"))
dn.rmsd.fel2 <- dn.ds.sum %>% filter(type == "bias", method == "FEL2", bl >= 0.01, ntaxa<=1024, rmsd <= 3) %>% ggplot(aes(x = factor(ntaxa), y = rmsd, fill = parameter)) + geom_violin(scale = "width") + facet_grid(~bl, scales = "free_x") + coord_cartesian(ylim=c(0,3)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("RMSD") + scale_fill_manual(values = dn_ds_colors, name = "Parameter", labels = c("dN", "dS"))
dn.ds.corr.violin.fel2 <- plot_grid(dn.r.fel2 , dn.rmsd.fel2, nrow=2, labels = c("A", "B"))
save_plot(paste0(PLOTDIR, "dn_ds_corr_violin_fel2.pdf"), dn.ds.corr.violin.fel2, base_width = 6.5, base_height=4.5)


dn.r.fubar2    <- dn.ds.sum %>% filter(type == "bias", method == "FUBAR2", bl >= 0.01, ntaxa<=1024) %>% ggplot(aes(x = factor(ntaxa), y = r, fill = parameter)) + geom_violin(scale = "width") + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(-0.1, 1), breaks = c(-0.1, 0.1, 0.3, 0.5, 0.7, 0.9)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("Correlation") + scale_fill_manual(values = dn_ds_colors, name = "Parameter", labels = c("dN", "dS"))
dn.rmsd.fubar2 <- dn.ds.sum %>% filter(type == "bias", method == "FUBAR2", bl >= 0.01, ntaxa<=1024, rmsd <= 3) %>% ggplot(aes(x = factor(ntaxa), y = rmsd, fill = parameter)) + geom_violin(scale = "width") + facet_grid(~bl, scales = "free_x") + coord_cartesian(ylim=c(0,3)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("RMSD") + scale_fill_manual(values = dn_ds_colors, name = "Parameter", labels = c("dN", "dS"))
dn.ds.corr.violin.fubar2 <- plot_grid(dn.r.fubar2 , dn.rmsd.fubar2, nrow=2, labels = c("A", "B"))
save_plot(paste0(PLOTDIR, "dn_ds_corr_violin_fubar2.pdf"), dn.ds.corr.violin.fubar2, base_width = 7, base_height=6)





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


linmodels.all %>% filter(modeltype == "r") %>% arrange(comp) %>%
  ggplot(aes(x = coeff, y = comp, color = sig)) + geom_point(size=ptsize) + 
  geom_segment(aes(x=lowerCI,xend=upperCI,y=comp,yend=comp),size=linesize) + 
  geom_vline(xintercept=0, size=0.5) + 
  facet_grid(~type) + background_grid() +
  xlab("Average Correlation Difference") + ylab("Comparison") + 
  scale_x_continuous(limits=c(-0.05, 0.3)) +  
  scale_color_manual(values = sig_colors) -> linmodel.r

linmodels.all %>% filter(modeltype == "rmsd") %>% arrange(comp) %>%
  ggplot(aes(x = coeff, y = comp, color = sig)) + geom_point(size=ptsize) + 
  geom_segment(aes(x=lowerCI,xend=upperCI,y=comp,yend=comp),size=linesize) + 
  geom_vline(xintercept=0, size=0.5) + 
  facet_grid(~type) + background_grid() +
  xlab("Average RMSD Difference") + ylab("Comparison") + 
  scale_x_continuous(limits=c(-0.15,0.12)) +  
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



r.dnds <- realdat.sum %>% ggplot(aes(x = dataset, fill = method, y = r)) + geom_boxplot(position=position_dodge(0.77), outlier.size = 0.75, size=0.4) + facet_wrap(~type) + scale_fill_manual(values=fill_colors, name = "Inference Method") +  xlab("Dataset") + ylab("Correlation") + scale_y_continuous(limits=c(0,1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.))
rmsd.dnds <- realdat.sum %>% ggplot(aes(x = dataset, fill = method, y = rmsd)) + geom_boxplot(position=position_dodge(0.77), outlier.size = 0.75, size=0.4) + facet_wrap(~type) + scale_fill_manual(values=fill_colors, name = "Inference Method") + xlab("Dataset") + ylab("RMSD")
#r.dn.ds <- dn.ds.real %>% ggplot(aes(x = dataset, fill = parameter, y = r)) + geom_boxplot(outlier.size = 0.75, size=0.4) + scale_fill_manual(values = fill_colors, name = "Parameter") + ylab("Correlation") + xlab("Dataset")
#rmsd.dn.ds <- dn.ds.real %>% ggplot(aes(x = dataset, fill = parameter, y = rmsd)) + geom_boxplot(outlier.size = 0.75, size=0.4) + scale_fill_manual(values = fill_colors, name = "Parameter") + ylab("RMSD") + xlab("Dataset")

real.plots <- plot_grid(r.dnds, rmsd.dnds, nrow=2, labels=c("A", "B"))
save_plot(paste0(PLOTDIR, "real_r_rmsd.pdf"), real.plots, base_width = 7.5, base_height = 4.5)



####################################################################################################
############### Lineplots of correlations and RMSD for balanced tree simulations ###################
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