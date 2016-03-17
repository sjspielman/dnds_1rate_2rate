# SJS
# Script to plot results

require(dplyr)
require(tidyr)
require(grid)
require(cowplot)
require(readr)


dat.mean$method <- factor(dat.mean$method, levels = c("FEL1", "SLAC1", "FUBAR1", "FEL2", "SLAC2", "FUBAR2"))
dat.mean <-dat %>% group_by(method, ntaxa, bl, truetype, type) %>% summarize(r.mean = mean(r), estbias.mean = mean(estbias), rmsd.mean = mean(rmsd))
dat.mean$method <- factor(dat.mean$method, levels = c("FEL1", "SLAC1", "FUBAR1", "FEL2", "SLAC2", "FUBAR2"))
line_colors <- c("red1", "black", "red4", "steelblue1", "green", "navy")
a<- dat.mean %>% filter(type == "nobias", truetype == "true1") %>% ggplot(aes(x = factor(ntaxa), y = r.mean, group = method, color=method)) +  geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 1), breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("Pearson Correlation") + scale_color_manual(values = line_colors, name = "Inference Method") + ggtitle("No bias")
b<- dat.mean %>% filter(type == "bias", truetype == "true2") %>% ggplot(aes(x = factor(ntaxa), y = r.mean, group = method, color=method)) +  geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 1), breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("Pearson Correlation") + scale_color_manual(values = line_colors, name = "Inference Method") + ggtitle("Bias")

plot_grid(a,b,nrow=2)
# 

# 
# dn.ds.corr %>% filter(type == "bias", method== "FEL2") %>% ggplot(aes(x = factor(ntaxa), y = r, fill=parameter)) +  geom_boxplot() + facet_grid(~bl, scales = "free_x") + coord_cartesian(ylim=c(-0.3, 1)) + scale_y_continuous(breaks = c(-0.1, 0.1, 0.3, 0.5, 0.7, 0.9)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("Pearson Correlation") + scale_fill_manual(values = c("blue", "yellow"))  + theme(legend.position = "none") -> p1
# 
# dn.ds.rmsd %>% filter(type == "bias", method== "FEL2", bl==0.16) %>% ggplot(aes(x = factor(ntaxa), y = rmsd, fill=parameter)) +  geom_boxplot() + coord_cartesian(ylim=c(0, 1.5)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("RMSD") + scale_fill_manual(values = c("blue", "yellow"))  + theme(legend.position = "none") -> p2
# 
# dn.ds.corr %>% filter(type == "bias", method== "FEL2") %>% ggplot(aes(x = factor(ntaxa), y = r, fill=parameter)) +  geom_boxplot() + scale_fill_manual(values = c("blue", "yellow"), name = "Parameter", labels=c("dN", "dS")) + theme(legend.position="bottom", legend.margin = unit(0, "cm"), legend.text = element_text(size=9), legend.key.size = unit(0.4, "cm")) -> legend.plot
# 
# grobs <- ggplotGrob(legend.plot)$grobs
# legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
# 
#           
# plot_grid(p1, p2, legend, nrow=2, rel_heights=c(1,1, 0.05))
# 
# 
# 
# 
# 
# 
# 
# line_colors <- c("red1", "red3", "red4", "steelblue1", "steelblue4", "navy")
# subline_colors <- c("red1", "red4", "steelblue1", "navy")
# sublineplot_methods <- c("FEL1", "FUBAR1", "FEL2",  "FUBAR2")
# theme_set(theme_cowplot() + theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 9, angle=30),  axis.title = element_text(size = 12), panel.border = element_rect(size = 0.5), panel.margin = unit(0.75, "lines"), strip.background = element_rect(fill="white"), strip.text = element_text(size=12)))
# 
# 
# 
# 
# dn.sum.mean <- dn.sum %>% group_by(method, ntaxa, bl, type) %>% summarize(r = mean(r1), estbias. = mean(estbias1), rmsd = mean(rmsd1))
# colnames(dn.sum.mean) <- c("method", "ntaxa", "bl", "type", "r.dn", "estbias.dn", "rmsd.dn")
# ds.sum.mean <- ds.sum %>% group_by(method, ntaxa, bl, type) %>% summarize(r = mean(r1), estbias = mean(estbias1), rmsd = mean(rmsd1))
# colnames(ds.sum.mean) <- c("method", "ntaxa", "bl", "type", "r.ds", "estbias.ds", "rmsd.ds")
# dnds.sum.mean$method <- factor(dnds.sum.mean$method, levels = c("FEL1", "SLAC1", "FUBAR1", "FEL2", "SLAC2", "FUBAR2"))
# dn.sum.mean$method <- factor(dn.sum.mean$method, levels = c("FEL1", "SLAC1", "FUBAR1", "FEL2", "SLAC2", "FUBAR2"))
# ds.sum.mean$method <- factor(ds.sum.mean$method, levels = c("FEL1", "SLAC1", "FUBAR1", "FEL2", "SLAC2", "FUBAR2"))
# 
# ### dN/dS correlation plots ###
# 
# dnds.sum.mean %>% filter(type == "nobias") %>% ggplot(aes(x = factor(ntaxa), y = r1, group = method, color=method)) +  geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 1), breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("Pearson Correlation") + scale_color_manual(values = line_colors, name = "Inference Method") + ggtitle("No bias")
# 
# dnds.sum.mean %>% filter(type == "bias") %>% ggplot(aes(x = factor(ntaxa), y = r1, group = method, color=method)) +  geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 1), breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("Pearson Correlation") + scale_color_manual(values = line_colors, name = "Inference Method") + ggtitle("Bias")
# 
# dnds.sum.mean %>% filter(type == "bias") %>% ggplot(aes(x = factor(ntaxa), y = r2, group = method, color=method)) +  geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 1), breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("Pearson Correlation") + scale_color_manual(values = line_colors, name = "Inference Method") + ggtitle("Bias")
# 
# 
# ### dN/dS RMSD plots ###
# 
# dnds.sum.mean %>% filter(type == "nobias") %>% ggplot(aes(x = factor(ntaxa), y = rmsd1, group = method, color=method)) +  geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + coord_cartesian(ylim = c(0,1)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("RMSD") + scale_color_manual(values = line_colors, name = "Inference Method") + ggtitle("No bias")
# 
# dnds.sum.mean %>% filter(type == "bias") %>% ggplot(aes(x = factor(ntaxa), y = rmsd1, group = method, color=method)) +  geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + coord_cartesian(ylim = c(0,1)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("RMSD1") + scale_color_manual(values = line_colors, name = "Inference Method") + ggtitle("Bias")
# 
# dnds.sum.mean %>% filter(type == "bias") %>% ggplot(aes(x = factor(ntaxa), y = rmsd2, group = method, color=method)) +  geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + coord_cartesian(ylim = c(0,1)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("RMSD2") + scale_color_manual(values = line_colors, name = "Inference Method") + ggtitle("Bias")
# 
# 
# ### dN correlation plots ###
# 
# dn.sum.mean %>% filter(type == "nobias", method %in% lineplot_methods) %>% ggplot(aes(x = factor(ntaxa), y = r.dn, group = method, color=method)) +  geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 1), breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("Pearson Correlation, dN") + scale_color_manual(values = subline_colors, name = "Inference Method") + ggtitle("No bias")
# 
# dn.sum.mean %>% filter(type == "bias", method %in% lineplot_methods) %>% ggplot(aes(x = factor(ntaxa), y = r.dn, group = method, color=method)) +  geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 1), breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("Pearson Correlation, dN") + scale_color_manual(values = subline_colors, name = "Inference Method") + ggtitle("Bias")
# 
# 
# ### dN RMSD plots ###
# 
# dn.sum.mean %>% filter(type == "nobias") %>% ggplot(aes(x = factor(ntaxa), y = rmsd.dn, group = method, color=method)) +  geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + coord_cartesian(ylim = c(0,1)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("RMSD, dN") + scale_color_manual(values = subline_colors, name = "Inference Method") + ggtitle("No bias")
# 
# dn.sum.mean %>% filter(type == "bias") %>% ggplot(aes(x = factor(ntaxa), y = rmsd.dn, group = method, color=method)) +  geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + coord_cartesian(ylim = c(0,1)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("RMSD, dN") + scale_color_manual(values = subline_colors, name = "Inference Method") + ggtitle("Bias")
# 
# 
# ## dN and dS separately correlation plots ##
# 
# 
# 
# 
# ### dS correlation plots ###
# 
# ds.sum.mean %>% filter(type == "bias", method %in% c("FEL2", "FUBAR2")) %>% ggplot(aes(x = factor(ntaxa), y = r.ds, group = method, color=method)) +  geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 1), breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("Pearson Correlation, dS") + scale_color_manual(values = c("orange", "darkgreen"), name = "Inference Method") + ggtitle("Bias")
# 
# 
# ### dS RMSD plots ###
# 
# ds.sum.mean %>% filter(type == "bias") %>% ggplot(aes(x = factor(ntaxa), y = rmsd.ds, group = method, color=method)) +  geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 4)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("RMSD, dS") + scale_color_manual(values = subline_colors, name = "Inference Method") + ggtitle("Bias")
# 
# 
# 
# 
# 



# Read in summary data frame results. 
results_directory <- "../results/processed_results/"
plot_directory <- "figures/"
dat.sum  <- read_csv(paste0(results_directory, "summarized_results.csv"))  # File containing summary dataframe with correlations and estimator bias for each method, dataset, replicate
dat.sum$method <- factor(dat.sum$method, levels = c("FEL1", "SLAC1", "FUBAR1", "FEL2", "SLAC2", "FUBAR2", "FEL2_1", "FUBAR2_1"))

dat.sum.mean <- dat.sum %>% group_by(method, ntaxa, bl, type) %>% summarize(r1 = mean(r1), r2 = mean(r2), estbias1 = mean(estbias1), estbias2 = mean(estbias2), rmsd1 = mean(rmsd1), rmsd2 = mean(rmsd2))

# Correlation line plots
line_colors <- c("red1", "red3", "red4", "steelblue1", "steelblue4", "navy")
lineplot_methods <- c("FEL1", "SLAC1", "FUBAR1", "FEL2", "SLAC2", "FUBAR2")
theme_set(theme_cowplot() + theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 9, angle=30),  axis.title = element_text(size = 12), panel.border = element_rect(size = 0.5), panel.margin = unit(0.75, "lines"), strip.background = element_rect(fill="white"), strip.text = element_text(size=12)))

nobias.corr <- dat.sum.mean %>% filter(type == "nobias", method %in% lineplot_methods) %>% ggplot(aes(x = factor(ntaxa), y = r2, group = method, color=method)) +  geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 1), breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("Pearson Correlation") + scale_color_manual(values = line_colors, name = "Inference Method") + ggtitle("No Codon bias")

bias.corr2 <- dat.sum.mean %>% filter(type == "bias", method %in% lineplot_methods) %>% ggplot(aes(x = factor(ntaxa), y = r2, group = method, color=method)) +  geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 1), breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("Pearson Correlation") + scale_color_manual(values = line_colors, name = "Inference Method") + ggtitle("Codon bias, True2")

bias.corr1 <- dat.sum.mean %>% filter(type == "bias", method %in% lineplot_methods) %>% ggplot(aes(x = factor(ntaxa), y = r1, group = method, color=method)) +  geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 1), breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("Pearson Correlation") + scale_color_manual(values = line_colors, name = "Inference Method") + ggtitle("Codon bias, True1")





corr.grid   <- plot_grid(nobias.corr, bias.corr2, bias.corr1, nrow = 3, labels = c("A", "B", "C"))
save_plot(paste0(plot_directory, "correlations_12rate.pdf"), corr.grid, base_width = 10, base_height=8)

# Estimator bias line plots
nobias.estbias <- dat.sum.mean %>% filter(type == "nobias", method %in% lineplot_methods) %>% ggplot(aes(x = factor(ntaxa), y = estbias2, group = method, color=method)) +  geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(-0.55, 2.1)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("Estimator Bias") + geom_hline(yintercept = 0) + scale_color_manual(values = line_colors, name = "Inference Method") + ggtitle("No Codon bias")

bias.estbias2  <- dat.sum.mean %>% filter(type == "bias", method %in% lineplot_methods) %>% ggplot(aes(x = factor(ntaxa), y = estbias2, group = method, color=method)) +  geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(-0.55, 2.1)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("Estimator Bias") + geom_hline(yintercept = 0) + scale_color_manual(values = line_colors, name = "Inference Method") + ggtitle("Codon bias, True2")

bias.estbias1  <- dat.sum.mean %>% filter(type == "bias", method %in% lineplot_methods) %>% ggplot(aes(x = factor(ntaxa), y = estbias1, group = method, color=method)) +  geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(-0.55, 2.1)) + background_grid(major = "xy") + xlab("Number of Taxa") + ylab("Estimator Bias") + geom_hline(yintercept = 0) + scale_color_manual(values = line_colors, name = "Inference Method") + ggtitle("Codon bias, True1")

estbias.grid   <- plot_grid(nobias.estbias, bias.estbias2, bias.estbias1, nrow = 3, labels = c("A", "B", "C"))
save_plot(paste0(plot_directory, "estimatorbias_12rate.pdf"), estbias.grid, base_width = 10, base_height=8)





# Linear model
dat <- read_csv("linear_model_results.csv") # Note that this file was manually created generated using information found in the file "build_linear_models.R"
comp_order <- c("FEL1 - FUBAR1", "SLAC1 - FUBAR1", "SLAC1 - FEL1", "FUBAR1 - FUBAR2_1", "FEL1 - FEL2_1", "FUBAR1 - FUBAR2", "SLAC1 - SLAC2", "FEL1 - FEL2")
dat$comp <- factor(dat$comp, levels = comp_order)

sig_colors <- c("grey60", "black")
theme_set(theme_cowplot() + theme(axis.title = element_text(size=15), axis.text.y = element_text(size=13), legend.position = "none"))
dat %>% filter(type == "nobias") %>% arrange(comp) %>% ggplot(aes(x = ave, y = comp, color = sig)) + geom_point(size=3) + geom_segment(aes(x=lwr,xend=upr,y=comp,yend=comp),size=0.75) + geom_vline(xintercept=0) + xlab("Average Correlation Difference") + ylab("Method Comparison") + background_grid() + scale_x_continuous(limits=c(-0.025, 0.2)) + ggtitle("No Codon bias") + scale_color_manual(values = sig_colors) -> multcomp.nobias

dat %>% filter(type == "bias2") %>% arrange(comp) %>% ggplot(aes(x = ave, y = comp, color = sig)) + geom_point(size=3) + geom_segment(aes(x=lwr,xend=upr,y=comp,yend=comp),size=0.75) + geom_vline(xintercept=0) + xlab("Average Correlation Difference") + ylab("Method Comparison") + background_grid() + scale_x_continuous(limits=c(-0.025, 0.2)) + ggtitle("Codon bias, True2") + scale_color_manual(values = sig_colors) -> multcomp.bias2

dat %>% filter(type == "bias1") %>% arrange(comp) %>% ggplot(aes(x = ave, y = comp, color = sig)) + geom_point(size=3) + geom_segment(aes(x=lwr,xend=upr,y=comp,yend=comp),size=0.75) + geom_vline(xintercept=0) + xlab("Average Correlation Difference") + ylab("Method Comparison") + background_grid() + scale_x_continuous(limits=c(-0.025, 0.2)) + ggtitle("Codon bias, True1") + scale_color_manual(values = sig_colors)-> multcomp.bias1

multcomp.plot <- plot_grid(multcomp.nobias, multcomp.bias2, multcomp.bias1, labels=c("A", "B", "C"),nrow=1, label_size=18, scale=0.95)
save_plot(paste0(plot_directory, "multcomp_plot_12rate.pdf"), multcomp.plot, base_width = 16, base_height=4)








# Mean error vs true dN/dS

nobias.full <- read_csv(paste0(results_directory, "full_results_nobias_dataset.csv"))
bias.full <- read_csv(paste0(results_directory, "full_results_bias_dataset.csv"))

theme_set(theme_cowplot() + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14), panel.border = element_rect(size = 0.5), panel.margin = unit(0.75, "lines"), strip.background = element_rect(fill="white"), strip.text = element_text(size=12), strip.text.y = element_text(angle=0)))

nobias.full %>% na.omit() %>% filter(method == "SLAC1") %>% group_by(ntaxa, bl, method, true2) %>%  summarize(meanerror = mean(abs(dnds-true2)/true2) ) -> slac1.error.nobias

bias.full %>% na.omit() %>% filter(method == "SLAC1") %>% group_by(ntaxa, bl, method, true2) %>% summarize(meanerror = mean(abs(dnds-true2)/true2) ) -> slac1.error.bias2

bias.full %>% na.omit() %>% filter(method == "SLAC1") %>% group_by(ntaxa, bl, method, true1) %>% summarize(meanerror = mean(abs(dnds-true1)/true1) ) -> slac1.error.bias1

slac1.nobias.slopes <- slac1.error.nobias %>% filter(ntaxa %in% c(256, 1024), bl %in% c(0.0025, 0.04, 0.64)) %>% group_by(ntaxa, bl) %>% do( raw.fit = lm(meanerror ~ true2, data=.)) %>% mutate(slope = paste0("slope = ",round(raw.fit[[1]][[2]],3))) %>% select(-raw.fit)

subset.slac1.error.plot <- slac1.error.nobias %>% filter(ntaxa %in% c(256, 1024), bl %in% c(0.0025, 0.04, 0.64)) %>% ggplot(aes(x = true2, y = meanerror)) + geom_point() + facet_grid(bl~ntaxa) + geom_hline(yintercept = 0, size = 0.75) + geom_smooth(color = "red", method = "lm", size=0.75, se = FALSE) + xlab("True dN/dS") + ylab("Average relative error\n") + panel_border() + geom_text(data = slac1.nobias.slopes, aes(label = slope), x = 0.75, y = 3.5, size = 3)
save_plot(paste0(plot_directory, "meanerror_vs_true_slac1_nobias_subset.pdf"), subset.slac1.error.plot, base_width = 6, base_height=5.5)

theme_set(theme_cowplot() + theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 12), axis.title = element_text(size = 14), panel.border = element_rect(size = 0.5), panel.margin = unit(0.75, "lines"), strip.background = element_rect(fill="white"), strip.text = element_text(size=12), strip.text.y = element_text(angle=0)))

full.slac1.error.nobias <- slac1.error.nobias %>% ggplot(aes(x = true2, y = meanerror)) + geom_point() + facet_grid(bl~ntaxa) + geom_hline(yintercept = 0, size = 0.75)  + geom_smooth(color = "red", method = "lm", size=0.75, se = FALSE) + xlab("True dN/dS") + ylab("Average relative error\n") + panel_border() + scale_y_continuous(limits=c(0,8)) + scale_x_continuous(limits=c(0,1))
save_plot(paste0(plot_directory, "meanerror_vs_true_slac1_nobias_full.pdf"), full.slac1.error.nobias, base_width = 8, base_height=6)

full.slac1.error.bias2 <- slac1.error.bias2 %>% ggplot(aes(x = true2, y = meanerror)) + geom_point() + facet_grid(bl~ntaxa) + geom_hline(yintercept = 0, size = 0.75) + geom_smooth(color = "red", method = "lm", size=0.75, se = FALSE) + xlab("True dN/dS") + ylab("Average relative error\n") + panel_border() + scale_y_continuous(limits=c(0,8)) + scale_x_continuous(limits=c(0,1))
save_plot(paste0(plot_directory, "meanerror_vs_true_slac1_bias2_full.pdf"), full.slac1.error.bias2 , base_width = 8, base_height=6)

full.slac1.error.bias1 <- slac1.error.bias1 %>% ggplot(aes(x = true1, y = meanerror)) + geom_point() + facet_grid(bl~ntaxa) + geom_hline(yintercept = 0, size = 0.75) + geom_smooth(color = "red", method = "lm", size=0.75, se = FALSE) + xlab("True dN/dS") + ylab("Average relative error\n") + panel_border() + scale_y_continuous(limits=c(0,8)) + scale_x_continuous(limits=c(0,1))
save_plot(paste0(plot_directory, "meanerror_vs_true_slac1_bias1_full.pdf"), full.slac1.error.bias1 , base_width = 8, base_height=6)







# Ntaxa vs branch length boxplot
theme_set(theme_cowplot() + theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 12),  axis.title = element_text(size = 13)))

dat.sum %>% mutate(treelen = 2*(ntaxa-1)*bl) %>% filter(method == "SLAC1", treelen > 162, treelen < 164) %>% na.omit() %>% group_by(treelen, type, ntaxa, bl, rep) -> treelen.dat


treelen.dat %>% select(-estbias1, -estbias2) %>% gather(rtype, r, r1, r2) %>% filter(type == "bias") -> treelen.dat.bias
treelen.dat %>% select(-estbias1, -estbias2) %>% gather(rtype, r, r1, r2) %>% filter(type == "nobias") -> treelen.dat.nobias
treelen.dat2 <- rbind(treelen.dat.nobias, treelen.dat.bias)
treelen.dat2$grouping <- "a"
treelen.dat2$grouping[treelen.dat2$type == "bias" & treelen.dat2$rtype == "r1"] <- "b"
treelen.dat2$grouping[treelen.dat2$type == "bias" & treelen.dat2$rtype == "r2"] <- "c"


bl.vs.ntaxa <- ggplot(treelen.dat2, aes(x = as.factor(bl), y=r, fill = as.factor(grouping))) + geom_boxplot(position = position_dodge(0.9)) + xlab("Parameterization") + scale_x_discrete(labels = c("N = 2048\nB = 0.04\n","N = 512\nB = 0.16\n","N = 128\nB = 0.64\n")) + ylab("Pearson Correlation") + scale_y_continuous(limits=c(0.45, 0.95)) + scale_fill_manual(name = "", labels = c("No Codon bias", "Codon bias, True1", "Codon bias, True2"), values = c("grey50", "grey80", "grey100"))
save_plot(paste0(plot_directory, "bl_ntaxa_boxplot_12rate.pdf"), bl.vs.ntaxa, base_width=7, base_height=5)

# summary(lm(r ~ grouping+as.factor(treelen), data=treelen.dat2))   # NOTE: interaction effect NS
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)               0.896063   0.003182 281.596   <2e-16 ***
# groupingb                -0.038968   0.003897  -9.999   <2e-16 ***
# groupingc                -0.088391   0.003897 -22.680   <2e-16 ***
# as.factor(treelen)163.52 -0.083929   0.003897 -21.536   <2e-16 ***
# as.factor(treelen)163.76 -0.189160   0.003897 -48.537   <2e-16 ***
