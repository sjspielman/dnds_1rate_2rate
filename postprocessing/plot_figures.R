require(dplyr)
require(tidyr)
require(grid)
require(cowplot)
require(readr)

PLOTDIR <- "figures/"

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

r.nobias <- dnds.sum.mean %>% filter(type == "nobias", truetype == "true1") %>% ggplot(aes(x = factor(ntaxa), y = r, group = method, color=method)) +  geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 1), breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) + xlab("Number of Taxa") + ylab("Pearson Correlation") + scale_color_manual(values = method_colors, name = "Inference Method") + background_grid("xy") + ggtitle("No codon bias")
r.bias <- dnds.sum.mean %>% filter(type == "bias", truetype == "true2") %>% ggplot(aes(x = factor(ntaxa), y = r, group = method, color=method)) +  geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 1), breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) + xlab("Number of Taxa") + ylab("Pearson Correlation") + scale_color_manual(values = method_colors, name = "Inference Method") + background_grid("xy") + ggtitle("Codon bias")
corr.lineplots <- plot_grid(r.nobias, r.bias, nrow=2, labels=c("A", "B"))
save_plot(paste0(PLOTDIR, "correlation_lineplots.pdf"), corr.lineplots, base_width = 10, base_height=5)

rmsd.nobias <- dnds.sum.mean %>% filter(type == "nobias", truetype == "true1", bl >=0.04, ntaxa >= 512) %>% ggplot(aes(x = factor(ntaxa), y = rmsd, group = method, color=method)) +  geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 1.25)) + xlab("Number of Taxa") + ylab("RMSD") + scale_color_manual(values = method_colors, name = "Inference Method") + background_grid("xy") + ggtitle("No codon bias")
rmsd.bias <- dnds.sum.mean %>% filter(type == "bias", truetype == "true2", bl >=0.04, ntaxa >= 512) %>% ggplot(aes(x = factor(ntaxa), y = rmsd, group = method, color=method)) +  geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 1.25)) + xlab("Number of Taxa") + ylab("RMSD") + scale_color_manual(values = method_colors, name = "Inference Method") + background_grid("xy") + ggtitle("Codon bias")
rmsd.lineplots <- plot_grid(rmsd.nobias, rmsd.bias, nrow=2, labels=c("A", "B"))
save_plot(paste0(PLOTDIR, "rmsd_lineplots.pdf"), rmsd.lineplots, base_width = 8, base_height=5)



############################################################################
######### dN is virtually the same between 1rate and 2rate models ##########
############################################################################
full.nobias <- read_csv("full_results_gtr_dataset.csv")
full.bias <- read_csv("full_results_bias_gtr_dataset.csv")


theme_set(theme_cowplot() + theme(axis.text.y = element_text(size = 9), 
                                  axis.text.x = element_text(size = 9, angle = 30),   
                                  axis.title = element_text(size = 14), 
                                  panel.border = element_rect(size = 0.5), 
                                  panel.margin = unit(1.0, "lines")))

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
bias.fubar.dn <- bias.compare.dn %>% ggplot(aes(x = FUBAR1, y = FUBAR2)) + geom_point() + facet_grid(ntaxa~bl) + geom_abline(slope=1, intercept=0, color="mediumvioletred") + scale_x_continuous(limits=c(0,1)) + scale_y_continuous(limits=c(0,1)) + xlab("dN, FEL1") + ylab("dN, FEL2")
bias.fel.dn <-bias.compare.dn %>% ggplot(aes(x = FEL1, y = FEL2)) + geom_point() + facet_grid(ntaxa~bl) + geom_abline(slope=1, intercept=0, color="mediumvioletred") + scale_x_continuous(limits=c(0,1)) + scale_y_continuous(limits=c(0,1)) + xlab("dN, FEL1") + ylab("dN, FEL2")

save_plot(paste0(PLOTDIR, "fubar12_nobias_samedn.pdf"), nobias.fubar.dn, base_width = 7, base_height=6)
save_plot(paste0(PLOTDIR, "fel12_nobias_samedn.pdf"), nobias.fel.dn, base_width = 7, base_height=6)
save_plot(paste0(PLOTDIR, "fubar12_bias_samedn.pdf"), bias.fubar.dn, base_width = 7, base_height=6)
save_plot(paste0(PLOTDIR, "fel12_bias_samedn.pdf"), bias.fel.dn, base_width = 7, base_height=6)


#########################################################################################################
############################### dN is more accurately estimated than dS #################################

theme_set(theme_cowplot() + theme(axis.text = element_text(size = 11), 
                                  axis.title = element_text(size = 12), 
                                  strip.background = element_rect(fill="white"), 
                                  strip.text = element_text(size=11), 
                                  panel.border = element_rect(size = 0.5), 
                                  panel.margin = unit(1.0, "lines")))


dn.r.fel2    <- dn.ds.sum %>% filter(type == "bias", method == "FEL2", bl >= 0.01, ntaxa<=1024) %>% ggplot(aes(x = factor(ntaxa), y = r, fill = parameter)) + geom_violin(scale = "width") + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(-0.1, 1), breaks = c(-0.1, 0.1, 0.3, 0.5, 0.7, 0.9)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("Pearson Correlation") + scale_fill_manual(values = c("grey80", "grey30"), name = "Parameter", labels = c("dN", "dS"))
dn.rmsd.fel2 <- dn.ds.sum %>% filter(type == "bias", method == "FEL2", bl >= 0.01, ntaxa<=1024, rmsd <= 2) %>% ggplot(aes(x = factor(ntaxa), y = rmsd, fill = parameter)) + geom_violin(scale = "width") + facet_grid(~bl, scales = "free_x") + coord_cartesian(ylim=c(0,2)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("RMSD") + scale_fill_manual(values = c("grey80", "grey30"), name = "Parameter", labels = c("dN", "dS"))
dn.ds.corr.violin.fel2 <- plot_grid(dn.r.fel2 , dn.rmsd.fel2, nrow=2, labels = c("A", "B"))
save_plot(paste0(PLOTDIR, "dn_ds_corr_violin_fel2.pdf"), dn.ds.corr.violin, base_width = 7, base_height=6)

dn.r.fubar2    <- dn.ds.sum %>% filter(type == "bias", method == "FUBAR2", bl >= 0.01, ntaxa<=1024) %>% ggplot(aes(x = factor(ntaxa), y = r, fill = parameter)) + geom_violin(scale = "width") + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(-0.1, 1), breaks = c(-0.1, 0.1, 0.3, 0.5, 0.7, 0.9)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("Pearson Correlation") + scale_fill_manual(values = c("grey80", "grey30"), name = "Parameter", labels = c("dN", "dS"))
dn.rmsd.fubar2 <- dn.ds.sum %>% filter(type == "bias", method == "FUBAR2", bl >= 0.01, ntaxa<=1024, rmsd <= 2) %>% ggplot(aes(x = factor(ntaxa), y = rmsd, fill = parameter)) + geom_violin(scale = "width") + facet_grid(~bl, scales = "free_x") + coord_cartesian(ylim=c(0,2)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("RMSD") + scale_fill_manual(values = c("grey80", "grey30"), name = "Parameter", labels = c("dN", "dS"))
dn.ds.corr.violin.fubar2 <- plot_grid(dn.r.fubar2 , dn.rmsd.fubar2, nrow=2, labels = c("A", "B"))
save_plot(paste0(PLOTDIR, "dn_ds_corr_violin_fubar2.pdf"), dn.ds.corr.violin.fubar2, base_width = 7, base_height=6)





##########################################################################################
sig_colors <- c("grey60", "black")
ptsize = 1.9
linesize = 0.9
comp.order <- c("FUBAR1 - SLAC1", "FEL1 - SLAC1", "FEL1 - FUBAR1", "SLAC1 - SLAC2", "FUBAR1 - FUBAR2", "FEL1 - FEL2")
theme_set(theme_cowplot() + theme(axis.text.x = element_text(size=10),
                                  axis.text.y = element_text(size=9),
                                  axis.title.x = element_text(size=11), 
                                  axis.title.y = element_text(size=12), 
                                  panel.border = element_rect(size = 0.5), 
                                  panel.margin = unit(1.0, "lines"), 
                                  strip.background = element_rect(fill="white"), 
                                  strip.text = element_text(size=12),
                                  legend.position = "none"))


linmodels <- read.csv("linear_model_results_correlation.csv") # Note that this file was manually created generated using information found in the file "build_linear_models.R"
linmodels2 <- linmodels %>% filter(comp %in% comp.order, type %in% c("nobias", "bias2"))
linmodels2$comp <- factor(linmodels2$comp, levels = comp.order)
linmodels2$type <- factor(linmodels2$type, levels=c("nobias", "bias2"), labels=c("No codon bias", "Codon bias"))


linmodels2 %>% arrange(comp) %>% 
           ggplot(aes(x = coeff, y = comp, color = sig)) + 
                  geom_point(size=ptsize) + 
                  geom_segment(aes(x=lowerCI,xend=upperCI,y=comp,yend=comp),size=linesize) + 
                  geom_vline(xintercept=0, size=0.5) + 
                  facet_grid(~type) + 
                  xlab("Average Correlation Difference") + 
                  ylab("Comparison") + 
                  background_grid() + 
                  scale_x_continuous(limits=c(-0.05, 0.25)) +  
                  scale_color_manual(values = sig_colors) -> linmodel.r


linmodels <- read.csv("linear_model_results_rmsd.csv") # Note that this file was manually created generated using information found in the file "build_linear_models.R"
linmodels2 <- linmodels %>% filter(comp %in% comp.order, type %in% c("nobias", "bias2"))
linmodels2$comp <- factor(linmodels2$comp, levels = comp.order)
linmodels2$type <- factor(linmodels2$type, levels=c("nobias", "bias2"), labels=c("No codon bias", "Codon bias"))
linmodels2 %>% arrange(comp) %>% 
           ggplot(aes(x = coeff, y = comp, color = sig)) + 
                  geom_point(size=ptsize) + 
                  geom_segment(aes(x=lowerCI,xend=upperCI,y=comp,yend=comp),size=linesize) + 
                  geom_vline(xintercept=0, size=0.5) + 
                  facet_grid(~type) + 
                  xlab("Average RMSD Difference") + 
                  ylab("Comparison") + 
                  background_grid() + 
                  scale_x_continuous(limits=c(-0.15, 0.05)) +  
                  scale_color_manual(values = sig_colors) -> linmodel.rmsd


linmodel.plots <- plot_grid(linmodel.r, linmodel.rmsd, nrow=2, labels=c("A", "B"))
save_plot(paste0(PLOTDIR, "linmodels.pdf"), linmodel.plots, base_width = 6.5, base_height = 4.5)

##########################################################################################


theme_set(theme_cowplot() + theme(axis.text = element_text(size=13),
                                  axis.title = element_text(size=14)))


bias.gtr.count <- read.csv("substitution_counts_bias_gtr.csv")
bias.gtr.count %>% group_by(ntaxa, bl, rep) %>% mutate(ratio_ns = ncount/scount) %>% na.omit() %>% filter(!is.infinite(ratio_ns)) %>% summarize(mean_ratio = mean(ratio_ns))  -> ns_count_ratio
count.ratio.box <- ggplot(ns_count_ratio, aes(x = as.factor(bl), y = mean_ratio, fill = as.factor(ntaxa))) + geom_boxplot(outlier.size = 0.75) + geom_hline(yintercept=1) + scale_fill_brewer(palette = "YlOrRd", name = "Number of Taxa") + xlab("Branch Length") + ylab("Mean N/S count ratio") +  scale_y_continuous(limits=c(0,4))
save_plot(paste0(PLOTDIR, "counted_ratio_bias.pdf"), count.ratio.box, base_width = 8, base_height = 5)



gtr.count <- read.csv("substitution_counts_gtr.csv")
gtr.count %>% group_by(ntaxa, bl, rep) %>% mutate(ratio_ns = ncount/scount) %>% na.omit() %>% filter(!is.infinite(ratio_ns)) %>% summarize(mean_ratio = mean(ratio_ns))  -> ns_count_ratio
count.ratio.box <- ggplot(ns_count_ratio, aes(x = as.factor(bl), y = mean_ratio, fill = as.factor(ntaxa))) + geom_boxplot(outlier.size = 0.75) + geom_hline(yintercept=1) + scale_fill_brewer(palette = "YlOrRd", name = "Number of Taxa") + xlab("Branch Length") + ylab("Mean N/S count ratio") + scale_y_continuous(limits=c(0,4))
save_plot(paste0(PLOTDIR, "counted_ratio_nobias.pdf"), count.ratio.box, base_width = 8, base_height = 5)

##########################################################################################
source("summary_functions.R")
realdat <- read_csv("full_results_realtrees.csv")


theme_set(theme_cowplot() + theme(axis.text    = element_text(size = 10),
                                  panel.border = element_rect(size = 0.5), 
                                  panel.margin = unit(1., "lines"), 
                                  strip.background = element_rect(fill="white"), 
                                  strip.text = element_text(size=11), 
                                  legend.text = element_text(size=10), 
                                  legend.title = element_text(size=11)))


realdat.sum <- realdat %>% summarize_dnds_real() 
realdat.sum$method <- factor(realdat.sum$method, levels = c("SLAC1", "SLAC2"))
realdat.sum$dataset <- factor(realdat.sum$dataset, levels = c("amine", "vertrho", "camelid", "h3", "hivrt"))
realdat.sum$type <- factor(realdat.sum$type, levels=c("gtr", "bias_gtr"), labels=c("No codon bias", "Codon bias"))
fill_colors <- c("darkred", "dodgerblue")


theme_set(theme_cowplot() + theme(axis.text    = element_text(size = 10),
                                  axis.title   = element_text(size = 12),
                                  panel.border = element_rect(size = 0.5), 
                                  panel.margin = unit(0.8, "lines"), 
                                  strip.background = element_rect(fill="white"), 
                                  strip.text = element_text(size=11), 
                                  legend.text = element_text(size=10), 
                                  legend.title = element_text(size=11)))



r <- realdat.sum %>% ggplot(aes(x = dataset, fill = method, y = r)) + geom_boxplot(position=position_dodge(0.77), outlier.size = 0.75, size=0.4) + facet_wrap(~type) + scale_fill_manual(values=fill_colors, name = "Inference Method") +  xlab("Dataset") + ylab("Pearson Correlation") + scale_y_continuous(limits=c(0,1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.))
rmsd <- realdat.sum %>% ggplot(aes(x = dataset, fill = method, y = rmsd)) + geom_boxplot(position=position_dodge(0.77), outlier.size = 0.75, size=0.4) + facet_wrap(~type) + scale_fill_manual(values=fill_colors, name = "Inference Method") + xlab("Dataset") + ylab("RMSD")
real.r.rmsd <- plot_grid(r,rmsd, nrow=2, labels=c("A", "B"))
save_plot(paste0(PLOTDIR, "real_r_rmsd.pdf"), real.r.rmsd, base_width = 8.75, base_height = 4.5)


realdat %>% filter(method %in% c("FUBAR2"), type == "bias_gtr") %>% summarize_dn_real() -> dn
realdat %>% filter(method %in% c("FUBAR2"), type == "bias_gtr") %>% summarize_ds_real() -> ds
dn.ds.real <- rbind(dn,ds)  
dn.ds.real$parameter <- factor(dn.ds.real$parameter, levels=c("dn", "ds"), labels=c("dN", "dS"))
dn.ds.real$dataset <- factor(dn.ds.real$dataset, levels = c("amine", "vertrho", "camelid", "h3", "hivrt"))

r.dn.ds <- dn.ds.real %>% ggplot(aes(x = dataset, fill = parameter, y = r)) + geom_boxplot(outlier.size = 0.75, size=0.4) + scale_fill_manual(values = fill_colors, name = "Parameter") + ylab("Pearson Correlation") + xlab("Dataset")
rmsd.dn.ds <- dn.ds.real %>% ggplot(aes(x = dataset, fill = parameter, y = rmsd)) + geom_boxplot(outlier.size = 0.75, size=0.4) + scale_fill_manual(values = fill_colors, name = "Parameter") + ylab("RMSD") + xlab("Dataset")
real.dn.ds.r.rmsd <- plot_grid(r.dn.ds, rmsd.dn.ds, nrow=2, labels=c("A", "B"))
save_plot(paste0(PLOTDIR, "real_dn_ds_r_rmsd.pdf"), real.dn.ds.r.rmsd, base_width = 6, base_height = 5)


