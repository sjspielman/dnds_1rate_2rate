require(dplyr)
require(tidyr)
require(grid)
require(cowplot)
require(readr)
require(ggrepel)
require(broom)


PLOTDIR <- "hayomfigures/"
DATADIR <- "dataframes/"
method_order  <- c("SLAC1", "FEL1", "FUBAR1", "SLAC2", "FEL2", "FUBAR2")
method_colors <- c("deeppink3", "red", "darkred", "skyblue", "dodgerblue", "blue")
method_colors_nofel2 <- c("deeppink3", "red", "darkred", "skyblue", "blue")

felfubar_colors <- c("red", "darkred", "dodgerblue", "blue")
dn_ds_colors <- c("red", "dodgerblue")


### Read in and factorize data ###
dnds.sum.mean  <- read_csv(paste0(DATADIR, "mean_summary_balanced_dnds.csv"))
dn.ds.sum.mean <- read_csv(paste0(DATADIR, "mean_summary_balanced_dn_ds.csv"))
dnds.sum       <- read_csv(paste0(DATADIR, "summary_balanced_dnds.csv"))
dn.ds.sum      <- read_csv(paste0(DATADIR, "summary_balanced_dn_ds.csv"))
real.dnds.sum  <- read_csv(paste0(DATADIR, "summary_real_dnds.csv"))
real.dn.ds.sum <- read_csv(paste0(DATADIR, "summary_real_dn_ds.csv"))
counted        <- read_csv(paste0(DATADIR, "substitution_counts.csv"))

dnds.sum.mean$method  <- factor(dnds.sum.mean$method, levels = method_order)
dn.ds.sum.mean$method <- factor(dn.ds.sum.mean$method, levels = method_order)
dnds.sum$method       <- factor(dnds.sum$method, levels = method_order)
dn.ds.sum$method      <- factor(dn.ds.sum$method, levels = method_order)
real.dnds.sum$method <- factor(real.dnds.sum$method, levels=c("FUBAR1", "FUBAR2") )
real.dn.ds.sum$method <- factor(real.dn.ds.sum$method, levels=c("FUBAR1", "FUBAR2") )
real.dnds.sum$dataset <- factor(real.dnds.sum$dataset, levels = c("amine", "vertrho", "camelid", "h3", "hivrt"))
real.dn.ds.sum$dataset <- factor(real.dn.ds.sum$dataset, levels = c("amine", "vertrho", "camelid", "h3", "hivrt"))

### Read in and factorize linear model results
comp.order <- c("FUBAR1 - SLAC1", "FEL1 - SLAC1", "FEL1 - FUBAR1", "SLAC1 - SLAC2", "FUBAR1 - FUBAR2", "FEL1 - FEL2")
linmodels <- read.csv(paste0(DATADIR, "linear_models.csv"))
linmodels <- filter(linmodels, comp %in% comp.order)
linmodels$model <- factor(linmodels$model, levels = c("r", "rmsd"))
linmodels$comp <- factor(linmodels$comp, levels = comp.order)
linmodels$biastype <- factor(linmodels$biastype, levels=c("nobias", "bias"), labels=c("No codon bias", "Codon bias"))
linmodels$pitype   <- factor(linmodels$pitype, levels=c("equalpi", "unequalpi"))

linmodels.heatmap <- read_csv(paste0(DATADIR, "linmodels_heatmap_version.csv"))
linmodels.heatmap$biastype <- factor(linmodels.heatmap$biastype, levels=c("nobias", "bias"), labels=c("No codon bias", "Codon bias"))

treelen.dat <- dnds.sum %>% mutate(treelen = 2*(ntaxa-1)*bl) %>% filter(method == "FEL1", treelen > 162, treelen < 164) %>% na.omit()

counted <- counted %>% mutate(nscount.ratio = ncount/scount) %>% select(-ncount, -scount) %>% na.omit() %>% filter(!is.infinite(nscount.ratio))

counted.sitewise <- counted %>%
                      group_by(pitype, biastype, site, truednds, ntaxa, bl) %>%
                      summarize(mean_ratio = mean(nscount.ratio))
counted.repwise <- counted %>%
                      group_by(pitype, biastype, rep, ntaxa, bl) %>%
                      summarize(mean_ratio = mean(nscount.ratio))

raw.bias <- read_csv(paste0(DATADIR, "results_balancedtrees_bias_unequalpi.csv"), col_types = list(
  pitype = col_character(),
  biastype = col_character(),
  dn = col_double(),
  ds = col_double(),
  dnds = col_double(),
  site = col_integer(),
  rep = col_integer(),
  ntaxa = col_integer(),
  bl = col_double(),
  method = col_factor(c("SLAC1", "FEL1", "FUBAR1", "SLAC2", "FEL2", "FUBAR2"))))
raw.nobias <- read_csv(paste0(DATADIR, "results_balancedtrees_nobias_unequalpi.csv"), col_types = list(
  pitype = col_character(),
  biastype = col_character(),
  dn = col_double(),
  ds = col_double(),
  dnds = col_double(),
  site = col_integer(),
  rep = col_integer(),
  ntaxa = col_integer(),
  bl = col_double(),
  method = col_factor(c("SLAC1", "FEL1", "FUBAR1", "SLAC2", "FEL2", "FUBAR2"))))

source("build_figure_legends.R") # Script **fully dependent on this one (including data loaded above!)** to make figure legends.





####################################################################################################
############ Lineplots of dN/dS correlations and RMSD for balanced tree simulations ################
####################################################################################################

### CORRELATION
r.nobias.equalpi <- dnds.sum.mean %>%
                     filter(pitype == "equalpi", biastype == "nobias") %>%
                     ggplot(aes(x = factor(ntaxa), y = r, group = method, color=method)) +
                     geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) +
                     facet_grid(~bl, scales = "free_x") +
                     scale_y_continuous(limits=c(0, 1), breaks = c(0.1, 0.3, 0.5, 0.7, 0.9), name = "Correlation") +
                     xlab("Number of Taxa") + ggtitle("No codon bias") +
                     scale_color_manual(values = method_colors, name = "Inference Method") +
                     background_grid("xy") +
                     theme(axis.text.x = element_text(size = 9, angle=30))
r.bias.equalpi <- dnds.sum.mean %>%
                   filter(pitype == "equalpi", biastype == "bias") %>%
                   ggplot(aes(x = factor(ntaxa), y = r, group = method, color=method)) +
                   geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) +
                   facet_grid(~bl, scales = "free_x") +
                   scale_y_continuous(limits=c(0, 1), breaks = c(0.1, 0.3, 0.5, 0.7, 0.9), name = "Correlation") +
                   xlab("Number of Taxa") + ggtitle("Codon bias") +
                   scale_color_manual(values = method_colors, name = "Inference Method") +
                   background_grid("xy") +
                   theme(axis.text.x = element_text(size = 9, angle=30))
corr.lineplots <- plot_grid(r.nobias.equalpi, r.bias.equalpi, nrow=2, labels=c("A", "B"))
save_plot(paste0(PLOTDIR, "correlation_lineplots_equalpi.pdf"), corr.lineplots, base_width = 10, base_height=5)

r.nobias.unequalpi <- dnds.sum.mean %>%
                       filter(pitype == "unequalpi", biastype == "nobias") %>%
                       ggplot(aes(x = factor(ntaxa), y = r, group = method, color=method)) +
                       geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(-0.1, 1), breaks = c(-0.1,0.1, 0.3, 0.5, 0.7, 0.9), name = "Correlation") +
                       xlab("Number of Taxa") + ggtitle("No codon bias") +
                       scale_color_manual(values = method_colors, name = "Inference Method") +
                       background_grid("xy") +
                       theme(axis.text.x = element_text(size = 9, angle=30))

r.bias.unequalpi <- dnds.sum.mean %>%
                    filter(pitype == "unequalpi", biastype == "bias") %>%
                    ggplot(aes(x = factor(ntaxa), y = r, group = method, color=method)) +
                    geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) +
                    facet_grid(~bl, scales = "free_x") +
                    scale_y_continuous(limits=c(-0.1, 1), breaks = c(-0.1, 0.1, 0.3, 0.5, 0.7, 0.9), name = "Correlation") +
                    xlab("Number of Taxa") + ggtitle("Codon bias") +
                    scale_color_manual(values = method_colors, name = "Inference Method") +
                    background_grid("xy") + theme(axis.text.x = element_text(size = 9, angle=30))
corr.lineplots <- plot_grid(r.nobias.unequalpi, r.bias.unequalpi, nrow=2, labels=c("A", "B"))
save_plot(paste0(PLOTDIR, "correlation_lineplots_unequalpi.pdf"), corr.lineplots, base_width = 10, base_height=5)




## Estimator Bias
estbias.nobias.equalpi <- dnds.sum.mean %>% filter(pitype == "equalpi", biastype == "nobias") %>% ggplot(aes(x = factor(ntaxa), y = estbias, group = method, color = method)) + geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(-0.5, 0.75)) + xlab("Number of Taxa") + ylab("Est. Bias") + scale_color_manual(values = method_colors, name = "Inference Method") + background_grid("xy") + ggtitle("No codon bias")  + theme(axis.text.x = element_text(size = 8, angle=30), legend.position = "none")
estbias.bias.equalpi <- dnds.sum.mean %>% filter(pitype == "equalpi", biastype == "bias") %>% ggplot(aes(x = factor(ntaxa), y = estbias, group = method, color = method)) + geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(-0.5, 0.75)) + xlab("Number of Taxa") + ylab("Est. Bias") + scale_color_manual(values = method_colors, name = "Inference Method") + background_grid("xy") + ggtitle("Codon bias") + theme(axis.text.x = element_text(size = 8, angle=30), legend.position = "none")
estbias.nobias.unequalpi <- dnds.sum.mean %>% filter(pitype == "unequalpi", biastype == "nobias") %>% ggplot(aes(x = factor(ntaxa), y = estbias, group = method, color = method)) + geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(-0.5, 0.75)) + xlab("Number of Taxa") + ylab("Est. Bias") + scale_color_manual(values = method_colors, name = "Inference Method") + background_grid("xy") + ggtitle("No codon bias")  + theme(axis.text.x = element_text(size = 8, angle=30), legend.position = "none")
estbias.bias.unequalpi <- dnds.sum.mean %>% filter(pitype == "unequalpi", biastype == "bias") %>% ggplot(aes(x = factor(ntaxa), y = estbias, group = method, color = method)) + geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(-0.5, 0.75)) + xlab("Number of Taxa") + ylab("Est. Bias") + scale_color_manual(values = method_colors, name = "Inference Method") + background_grid("xy") + ggtitle("Codon bias")  + theme(axis.text.x = element_text(size = 8, angle=30), legend.position = "none")
estbias.grid <- plot_grid(estbias.nobias.equalpi, estbias.bias.equalpi, estbias.nobias.unequalpi, estbias.bias.unequalpi, methods.legend,  nrow=5, labels=c("A", "B", "C", "D"), rel_heights=c(1,1,1,1,0.17)) + theme(axis.text.x = element_text(size = 9, angle=30))
save_plot(paste0(PLOTDIR, "estbias_lineplots.pdf"), estbias.grid, base_width = 7, base_height=9)




## RMSD
rmsd.nobias.equalpi <- dnds.sum.mean %>% filter(pitype == "equalpi", biastype == "nobias", bl >=0.01) %>% ggplot(aes(x = factor(ntaxa), y = rmsd, group = method, color = method)) + geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 1.5)) + xlab("Number of Taxa") + ylab("RMSD") + scale_color_manual(values = method_colors, name = "Inference Method") + background_grid("xy") + ggtitle("No codon bias") + theme(axis.text.x = element_text(size = 9, angle=30), legend.position = "none")
rmsd.bias.equalpi <- dnds.sum.mean %>% filter(pitype == "equalpi", biastype == "bias", bl >=0.01) %>% ggplot(aes(x = factor(ntaxa), y = rmsd, group = method, color = method)) + geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 1.5)) + xlab("Number of Taxa") + ylab("RMSD") + scale_color_manual(values = method_colors, name = "Inference Method") + background_grid("xy") + ggtitle("Codon bias")+ theme(axis.text.x = element_text(size = 9, angle=30), legend.position = "none")
rmsd.nobias.unequalpi <- dnds.sum.mean %>% filter(pitype == "unequalpi", biastype == "nobias", bl >=0.01) %>% ggplot(aes(x = factor(ntaxa), y = rmsd, group = method, color = method)) + geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 1.5)) + xlab("Number of Taxa") + ylab("RMSD") + scale_color_manual(values = method_colors, name = "Inference Method") + background_grid("xy") + ggtitle("No codon bias")+ theme(axis.text.x = element_text(size = 9, angle=30), legend.position = "none")
rmsd.bias.unequalpi <- dnds.sum.mean %>% filter(pitype == "unequalpi", biastype == "bias", bl >=0.01) %>% ggplot(aes(x = factor(ntaxa), y = rmsd, group = method, color = method)) + geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 1.5)) + xlab("Number of Taxa") + ylab("RMSD") + scale_color_manual(values = method_colors, name = "Inference Method") + background_grid("xy") + ggtitle("Codon bias")+ theme(axis.text.x = element_text(size = 9, angle=30), legend.position = "none")
rmsd.lineplots <- plot_grid(rmsd.nobias.equalpi, rmsd.bias.equalpi, rmsd.nobias.unequalpi, rmsd.bias.unequalpi,methods.legend2, nrow=5, labels=c("A", "B", "C", "D"), rel_heights=c(1,1,1,1,0.17))
save_plot(paste0(PLOTDIR, "rmsd_lineplots.pdf"), rmsd.lineplots, base_width = 7, base_height=8.5)


## Residual variance
resvar.nobias.equalpi <- dnds.sum.mean %>% filter(pitype == "equalpi", biastype == "nobias", bl >=0.01) %>% ggplot(aes(x = factor(ntaxa), y = resvar, group = method, color = method)) + geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 2.5)) + xlab("Number of Taxa") + ylab("Res. Variance") + scale_color_manual(values = method_colors, name = "Inference Method") + background_grid("xy") + ggtitle("No codon bias") + theme(axis.text.x = element_text(size = 9, angle=30), legend.position = "none")
resvar.bias.equalpi <- dnds.sum.mean %>% filter(pitype == "equalpi", biastype == "bias", bl >=0.01) %>% ggplot(aes(x = factor(ntaxa), y = resvar, group = method, color = method)) + geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 2.5)) + xlab("Number of Taxa") + ylab("Res. Variance") + scale_color_manual(values = method_colors, name = "Inference Method") + background_grid("xy") + ggtitle("Codon bias") + theme(axis.text.x = element_text(size = 9, angle=30), legend.position = "none")
resvar.nobias.unequalpi <- dnds.sum.mean %>% filter(pitype == "unequalpi", biastype == "nobias", bl >=0.01) %>% ggplot(aes(x = factor(ntaxa), y = resvar, group = method, color = method)) + geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 2)) + xlab("Number of Taxa") + ylab("Res. Variance") + scale_color_manual(values = method_colors, name = "Inference Method") + background_grid("xy") + ggtitle("No codon bias") + theme(axis.text.x = element_text(size = 9, angle=30), legend.position = "none")
resvar.bias.unequalpi <- dnds.sum.mean %>% filter(pitype == "unequalpi", biastype == "bias", bl >=0.01) %>% ggplot(aes(x = factor(ntaxa), y = resvar, group = method, color = method)) + geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + facet_grid(~bl, scales = "free_x") + scale_y_continuous(limits=c(0, 2)) + xlab("Number of Taxa") + ylab("Res. Variance") + scale_color_manual(values = method_colors, name = "Inference Method") + background_grid("xy") + ggtitle("Codon bias") + theme(axis.text.x = element_text(size = 9, angle=30), legend.position = "none")

resvar.lineplots <- plot_grid(resvar.nobias.equalpi, resvar.bias.equalpi, resvar.nobias.unequalpi, resvar.bias.unequalpi, nrow=5, methods.legend, labels=c("A", "B", "C", "D"), rel_heights=c(1,1,1,1,0.15))
save_plot(paste0(PLOTDIR, "resvar_lineplots.pdf"), resvar.lineplots, base_width = 7, base_height=8.5)



####################################################################################################
######### Violin plots  of dN & dS correlations and RMSD for balanced tree simulations #############
####################################################################################################

param.r.fel2.sub <- dn.ds.sum %>%
                      filter(pitype == "unequalpi", biastype == "bias", method == "FEL2", bl >= 0.01, ntaxa<=1024) %>%
                      ggplot(aes(x = factor(ntaxa), y = r, fill = parameter)) +
                      geom_violin(scale = "width") +
                      facet_grid(~bl, scales = "free_x") +
                      scale_y_continuous(limits=c(-0.1, 1), breaks = c(-0.1, 0.1, 0.3, 0.5, 0.7, 0.9), name = "Correlation") +
                      background_grid(major = "xy") +
                      xlab("Number of Taxa") +
                      scale_fill_manual(values = dn_ds_colors) +
                      theme(legend.position = "none")
param.rmsd.fel2.sub <- dn.ds.sum %>%
                        filter(pitype == "unequalpi", biastype == "bias", method == "FEL2", bl >= 0.01, ntaxa<=1024, rmsd <= 5) %>%
                        ggplot(aes(x = factor(ntaxa), y = rmsd, fill = parameter)) +
                        geom_violin(scale = "width") +
                        facet_grid(~bl, scales = "free_x") +
                        coord_cartesian(ylim=c(0,3)) +
                        background_grid(major = "xy") +
                        xlab("Number of Taxa") +
                        ylab("RMSD") +
                        scale_fill_manual(values = dn_ds_colors) +
                        theme(legend.position = "none")

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




param.r.fel2.equal <- dn.ds.sum %>%
                       filter(pitype == "equalpi", biastype == "bias", method == "FEL2", bl>=0.01) %>%
                       ggplot(aes(x = factor(ntaxa), y = r, fill = parameter)) +
                       geom_violin(scale = "width") +
                       facet_grid(~bl, scales = "free_x") +
                       scale_y_continuous(limits=c(-0.1, 1), breaks = c(-0.1, 0.1, 0.3, 0.5, 0.7, 0.9)) +
                       background_grid(major = "xy") +
                       xlab("Number of Taxa") + ylab("Correlation") +
                       scale_fill_manual(values = dn_ds_colors) +
                       theme(axis.text.x = element_text(size = 8, angle=30), legend.position = "none")
param.r.fel2.unequal <- dn.ds.sum %>%
                         filter(pitype == "unequalpi", biastype == "bias", method == "FEL2", bl>=0.01) %>%
                         ggplot(aes(x = factor(ntaxa), y = r, fill = parameter)) +
                         geom_violin(scale = "width") +
                         facet_grid(~bl, scales = "free_x") +
                         scale_y_continuous(limits=c(-0.1, 1), breaks = c(-0.1, 0.1, 0.3, 0.5, 0.7, 0.9)) +
                         background_grid(major = "xy") +
                         xlab("Number of Taxa") + ylab("Correlation") +
                         scale_fill_manual(values = dn_ds_colors) +
                         theme(axis.text.x = element_text(size = 8, angle=30), legend.position = "none")
param.rmsd.fel2.equal <- dn.ds.sum %>%
                          filter(pitype == "equalpi", biastype == "bias", method == "FEL2", bl>=0.01) %>%
                          ggplot(aes(x = factor(ntaxa), y = rmsd, fill = parameter)) +
                          geom_violin(scale = "width") +
                          facet_grid(~bl, scales = "free_x") +
                          coord_cartesian(ylim=c(0,5)) +
                          background_grid(major = "xy") +
                          xlab("Number of Taxa") + ylab("RMSD") +
                          scale_fill_manual(values = dn_ds_colors) +
                          theme(axis.text.x = element_text(size = 8, angle=30), legend.position = "none")
param.rmsd.fel2.unequal <- dn.ds.sum %>%
                            filter(pitype == "unequalpi", biastype == "bias", method == "FEL2", bl>=0.01) %>%
                            ggplot(aes(x = factor(ntaxa), y = rmsd, fill = parameter)) +
                            geom_violin(scale = "width") +
                            facet_grid(~bl, scales = "free_x") +
                            coord_cartesian(ylim=c(0,5)) +
                            background_grid(major = "xy") +
                            xlab("Number of Taxa") + ylab("RMSD") +
                            scale_fill_manual(values = dn_ds_colors) +
                            theme(axis.text.x = element_text(size = 8, angle=30), legend.position = "none")

param.r.fel2 <- plot_grid(param.r.fel2.equal, param.r.fel2.unequal, nrow=1, labels=c("A", "B"), scale=0.98)
param.rmsd.fel2 <- plot_grid(param.rmsd.fel2.equal, param.rmsd.fel2.unequal, nrow=1, labels=c("C", "D"),scale=0.98)
param.r.rmsd.fel2 <- plot_grid(param.r.fel2, param.rmsd.fel2, dn.ds.legend, nrow=3, rel_heights=c(1,1,0.1))
save_plot(paste0(PLOTDIR, "dn_ds_r_rmsd_fel2.pdf"), param.r.rmsd.fel2, base_width = 10, base_height=5)

####################################################################################################
###################################### Linear models lineplot ######################################
####################################################################################################
sig_colors <- c("grey60", "black")
ptsize = 1.5
linesize = 0.8

linmodels %>% filter(model == "r", pitype == "unequalpi") %>% arrange(comp) %>%
  ggplot(aes(x = coeff, y = comp, color = sig)) + geom_point(size=ptsize) +
  geom_segment(aes(x=lowerCI,xend=upperCI,y=comp,yend=comp),size=linesize) +
  geom_vline(xintercept=0, size=0.5) +
  facet_grid(~biastype) + background_grid() +
  xlab("Average Correlation Difference") + ylab("Comparison") +
  scale_x_continuous(limits=c(-0.025, 0.2)) +
  scale_color_manual(values = sig_colors)  + theme(legend.position = "none", axis.text.y = element_text(size=10)) -> linmodel.unequal.r

linmodels %>% filter(model == "rmsd", pitype == "unequalpi") %>% arrange(comp) %>%
  ggplot(aes(x = coeff, y = comp, color = sig)) + geom_point(size=ptsize) +
  geom_segment(aes(x=lowerCI,xend=upperCI,y=comp,yend=comp),size=linesize) +
  geom_vline(xintercept=0, size=0.5) +
  facet_grid(~biastype) + background_grid() +
  xlab("Average RMSD Difference") + ylab("Comparison") +
  scale_x_continuous(limits=c(-0.2,0.1)) +
  scale_color_manual(values = sig_colors) + theme(legend.position = "none", axis.text.y = element_text(size=10)) -> linmodel.unequal.rmsd

linmodel.unequal <- plot_grid(linmodel.unequal.r, linmodel.unequal.rmsd, nrow=2, labels=c("A", "B"))
save_plot(paste0(PLOTDIR, "linmodels_unequal.pdf"), linmodel.unequal, base_width = 6.5, base_height = 5)


linmodels %>% filter(model == "r", pitype == "equalpi") %>% arrange(comp) %>%
    ggplot(aes(x = coeff, y = comp, color = sig)) + geom_point(size=ptsize) +
    geom_segment(aes(x=lowerCI,xend=upperCI,y=comp,yend=comp),size=linesize) +
    geom_vline(xintercept=0, size=0.5) +
    facet_grid(~biastype) + background_grid() +
    xlab("Average Correlation Difference") + ylab("Comparison") +
    scale_x_continuous(limits=c(-0.025, 0.2)) +
    scale_color_manual(values = sig_colors)  + theme(legend.position = "none", axis.text.y = element_text(size=10))  -> linmodel.equal.r

linmodels %>% filter(model == "rmsd", pitype == "equalpi") %>% arrange(comp) %>%
    ggplot(aes(x = coeff, y = comp, color = sig)) + geom_point(size=ptsize) +
    geom_segment(aes(x=lowerCI,xend=upperCI,y=comp,yend=comp),size=linesize) +
    geom_vline(xintercept=0, size=0.5) +
    facet_grid(~biastype) + background_grid() +
    xlab("Average RMSD Difference") + ylab("Comparison") +
    scale_x_continuous(limits=c(-0.25,0.1)) +
    scale_color_manual(values = sig_colors)  + theme(legend.position = "none", axis.text.y = element_text(size=10))   -> linmodel.equal.rmsd

linmodel.equal <- plot_grid(linmodel.equal.r, linmodel.equal.rmsd, nrow=1, labels=c("A", "B"), scale=0.98)
save_plot(paste0(PLOTDIR, "linmodels_equal.pdf"), linmodel.equal, base_width = 12, base_height = 3)



####################################################################################################
################################### Linear models heatmaps #########################################
####################################################################################################

na.value <- "grey90"
linmodels.heatmap %>%
  filter(pitype == "unequalpi") %>%
  ggplot(aes(x = as.factor(ntaxa), y = as.factor(bl), fill = r.diff)) +
    geom_tile() +
    geom_text(aes(label = round(r.diff,3), fontface = sig+1)) +
    scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "red", na.value = na.value) +
    facet_grid(method~biastype) +
    xlab("Number of Taxa") + ylab("Branch lengths") +
    theme(legend.position = "none")-> unequal.heatmaps.true

linmodels.heatmap %>%
  filter(pitype == "equalpi") %>%
  ggplot(aes(x = as.factor(ntaxa), y = as.factor(bl), fill = r.diff)) +
    geom_tile() +
    geom_text(aes(label = round(r.diff,3), fontface = sig+1)) +
    scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "red", na.value = na.value) +
    facet_grid(method~biastype) +
    xlab("Number of Taxa") + ylab("Branch lengths") +
    theme(legend.position = "none")-> equal.heatmaps.true

heatmap.true.grid <- plot_grid(equal.heatmaps.true, unequal.heatmaps.true, nrow=1, labels="AUTO", label_size=17)
heatmap.true.grid2 <- plot_grid(heatmap.true.grid, heatmap.legend, nrow=2, rel_heights=c(1,0.1), rel_widths=c(1,0.5))
save_plot(paste0(PLOTDIR, "heatmap.pdf"), heatmap.true.grid2, base_width = 16, base_height=7)



####################################################################################################
############################### Violins of optimized branch lengths  ###############################
####################################################################################################

bl <- read.csv(paste0(DATADIR, "optimized_bl.csv"))
bl$biastype <- factor(bl$biastype, levels=c("nobias", "bias"), labels=c("No codon bias", "Codon bias"))
opt.bl <- ggplot(bl, aes(y = meanbl, x = factor(ntaxa))) +
           facet_grid(biastype~bl) +
           geom_violin(alpha=0.8, scale = "width", fill="grey70") +
           xlab("Number of Taxa") + ylab("Mean optimized branch length") +
           theme(axis.text.x = element_text(angle=30, size=10))
save_plot(paste0(PLOTDIR, "optimized_bl_violins.pdf"), opt.bl, base_width = 8, base_height=4.5)

####################################################################################################
########################### Scatterplots of dN vs. dN and dN/dS vs. dN/dS  #########################
####################################################################################################

raw.nobias %>% dplyr::select(-dnds, -ds) %>%
  na.omit() %>%
  group_by(ntaxa, bl, site, method) %>%
  summarize(mean.dn = mean(dn)) %>%
  spread(method, mean.dn) -> nobias.compare.dn
raw.bias %>% dplyr::select(-dnds, -ds) %>%
  na.omit() %>%
  group_by(ntaxa, bl, site, method) %>%
  summarize(mean.dn = mean(dn)) %>%
  spread(method, mean.dn) -> bias.compare.dn

nobias.fel.dn <- nobias.compare.dn %>%
                  ggplot(aes(x = FEL1, y = FEL2)) +
                  geom_point() +
                  facet_grid(ntaxa~bl) +
                  geom_abline(slope=1, intercept=0, color="red", size=0.7) +
                  scale_x_continuous(limits=c(0,1), breaks=c(0, 0.5, 1), name = "dN, FEL1") +
                  scale_y_continuous(limits=c(0,1), breaks=c(0, 0.5, 1), name = "dN, FEL2") +
                  ggtitle("No codon bias")
bias.fel.dn   <- bias.compare.dn %>%
                  ggplot(aes(x = FEL1, y = FEL2)) +
                  geom_point() +
                  facet_grid(ntaxa~bl) +
                  geom_abline(slope=1, intercept=0, color="red", size=0.7) +
                  scale_x_continuous(limits=c(0,1), breaks=c(0, 0.5, 1), name = "dN, FEL1") +
                  scale_y_continuous(limits=c(0,1), breaks=c(0, 0.5, 1), name = "dN, FEL2") +
                  ggtitle("Codon bias")
dn.same.grid <- plot_grid(nobias.fel.dn, bias.fel.dn, nrow=1, labels=c("A", "B"), label_size=20)
save_plot(paste0(PLOTDIR, "samedn.pdf"), dn.same.grid, base_width = 10, base_height=5)


raw.nobias %>% dplyr::select(-dn, -ds) %>%
  na.omit() %>%
  group_by(ntaxa, bl, site, method) %>%
  summarize(mean.dnds = mean(dnds)) %>%
  spread(method, mean.dnds) -> nobias.compare.dnds
raw.bias  %>% dplyr::select(-dn, -ds) %>%
  na.omit() %>%
  group_by(ntaxa, bl, site, method) %>%
  summarize(mean.dnds = mean(dnds)) %>%
  spread(method, mean.dnds) -> bias.compare.dnds

nobias.fel.dnds <- nobias.compare.dnds %>% filter(bl >= 0.01) %>%
                  ggplot(aes(x = FEL1, y = FEL2)) +
                  geom_point() +
                  facet_grid(ntaxa~bl) +
                  geom_abline(slope=1, intercept=0, color="red", size=0.7) +
                  scale_x_continuous(limits=c(0,1), breaks=c(0, 0.5, 1), name = "dN/dS, FEL1") +
                  scale_y_continuous(limits=c(0,1), breaks=c(0, 0.5, 1), name = "dN/dS, FEL2") +
                  ggtitle("No codon bias")
bias.fel.dnds  <- bias.compare.dn %>% filter(bl >= 0.01) %>%
                  ggplot(aes(x = FEL1, y = FEL2)) +
                  geom_point() +
                  facet_grid(ntaxa~bl) +
                  geom_abline(slope=1, intercept=0, color="red", size=0.7) +
                  scale_x_continuous(limits=c(0,1), breaks=c(0, 0.5, 1), name = "dN/dS, FEL1") +
                  scale_y_continuous(limits=c(0,1), breaks=c(0, 0.5, 1), name = "dN/dS, FEL2") +
                  ggtitle("Codon bias")
dnds.same.grid <- plot_grid(nobias.fel.dnds, bias.fel.dnds, nrow=1, labels=c("A", "B"), label_size=20)
save_plot(paste0(PLOTDIR, "samednds.pdf"), dnds.same.grid, base_width = 10, base_height=5)



####################################################################################################
############################### Density of inferred dS estimates for FEL2  #########################
####################################################################################################

raw.bias %>% filter(method == "FEL2") %>% dplyr::select(biastype, ds, site, rep, ntaxa, bl)-> fel2.bias
raw.nobias %>% filter(method == "FEL2") %>% dplyr::select(biastype, ds, site, rep, ntaxa, bl) %>% rbind(fel2.bias) %>% filter(ds <= 3)-> fel2
fel2$biastype <- factor(fel2$biastype, levels=c("nobias", "bias"), labels=c("No codon bias", "Codon bias"))
fel2 %>% filter(bl >= 0.01, ntaxa >= 256) %>%
  ggplot(aes(x = ds, fill = biastype, group = biastype)) +
  geom_density(alpha=0.6) +
  facet_grid(ntaxa~bl) +
  scale_fill_manual(values = dn_ds_colors, name = "Simulation set") +
  xlab("dS estimate by FEL2") + ylab("Density") -> fel2.ds.density
save_plot(paste0(PLOTDIR, "fel2_ds_density.pdf"), fel2.ds.density, base_width = 8, base_height=5)




####################################################################################################
################################# Boxplots of counted substitutions ################################
####################################################################################################

unequal.nobias <- counted.repwise %>% filter(biastype == "nobias") %>%
                    ggplot(aes(x = as.factor(ntaxa), y = mean_ratio, fill = as.factor(bl))) +
                    geom_boxplot(outlier.size = 0.75) +
                    geom_hline(yintercept=1) +
                    scale_color_hue(l=40, name = "Branch Lengths") +
                    xlab("Number of Taxa") + ylab("Mean N/S substitutions") +
                    ggtitle("No codon bias") +
                    scale_y_continuous(limits=c(0,4)) +
                    background_grid("xy")+
                    theme(legend.position = "none")
unequal.bias <- counted.repwise %>% filter(biastype == "bias") %>%
                    ggplot(aes(x = as.factor(ntaxa), y = mean_ratio, fill = as.factor(bl))) +
                    geom_boxplot(outlier.size = 0.75) +
                    geom_hline(yintercept=1) +
                    scale_color_hue(l=40, name = "Branch Lengths") +
                    xlab("Number of Taxa") + ylab("Mean N/S substitutions") +
                    ggtitle("Codon bias") +
                    scale_y_continuous(limits=c(0,4)) +
                    background_grid("xy")+
                    theme(legend.position = "none")
unequal.count <- plot_grid(unequal.nobias, unequal.bias, nrow=1, labels=c("A", "B"), scale=0.95)
unequal.count.withlegend <- plot_grid(unequal.count, count.legend, nrow=2, rel_heights=c(1,0.1))
save_plot(paste0(PLOTDIR, "counted_ratio_boxplots.pdf"), unequal.count.withlegend, base_width = 8, base_height = 3.5)







####################################################################################################
############################# Counted substitutions vs true dN/dS ##################################
####################################################################################################

### Plot of *site-specific* N:S ratios. Small plot for maintext and larger for SI

counted.sitewise %>%
  filter(pitype == "unequalpi", biastype == "nobias") %>%
  ggplot(aes(x = truednds, group=truednds, y = mean_ratio)) + geom_point() +
  facet_grid(ntaxa~bl) +
  geom_hline(yintercept=1, color="red") +
  xlab("True dN/dS") + ylab("Mean N/S substitutions") + ggtitle("No codon bias") -> p.nobias
counted.sitewise %>%
  filter(pitype == "unequalpi", biastype == "bias") %>%
  ggplot(aes(x = truednds, group=truednds, y = mean_ratio)) + geom_point() +
  facet_grid(ntaxa~bl) +
  geom_hline(yintercept=1, color="red") +
  xlab("True dN/dS") + ylab("Mean N/S substitutions") + ggtitle("Codon bias") -> p.bias
p.grid <- plot_grid(p.nobias, p.bias, nrow=1)
save_plot(paste0(PLOTDIR,"sitewisecounts_vs_truednds_full.pdf"), p.grid, base_width = 12, base_height=5)


counted.sitewise %>% filter(pitype == "unequalpi", bl == 0.04, ntaxa == 512) -> sub.counted.sitewise
sub.counted.sitewise$biastype <- factor(sub.counted.sitewise$biastype, levels=c("nobias", "bias"), labels=c("No codon bias", "Codon bias"))
sub.counted.sitewise %>%
  ggplot(aes(x = truednds, group=truednds, y = mean_ratio)) + geom_point() +
  facet_grid(~biastype) +
  geom_hline(yintercept=1, color="red") +
  xlab("True dN/dS") + ylab("Mean N/S substitutions") -> p.sub
save_plot(paste0(PLOTDIR,"sitewisecounts_vs_truednds_subset.pdf"), p.sub, base_width = 6, base_height=3)




####################################################################################################
############################# Counted substitutions vs true dN/dS INTERCEPT ########################
####################################################################################################


counted.sitewise %>%
  filter(mean_ratio > 1e-10, pitype == "unequalpi") %>%
  group_by(biastype, ntaxa, bl) %>%
  do(tidy(lm(truednds ~ mean_ratio, data=.))) %>%
  dplyr::select(-std.error, -statistic, -p.value) %>%
  summarize(intercept = sum(estimate)) -> dnds.nsratio.intercept
dnds.nsratio.intercept$biastype <- factor(dnds.nsratio.intercept$biastype, levels=c("nobias", "bias"), labels=c("No codon bias", "Codon bias"))

dnds.nsratio.intercept %>%
  ggplot(aes(x = as.factor(ntaxa), y = intercept, group = as.factor(bl), color = as.factor(bl))) +
    geom_line(size=0.75, alpha=0.8) + geom_point(size=2.75, alpha=0.8) +
    facet_grid(~biastype) +
    scale_color_hue(l=40, name = "Branch Lengths     ") +
    background_grid("xy") +
    scale_y_continuous(limits=c(0.2, 1.1)) +
    xlab("Number of Taxa") + ylab("dN/dS at Intercept") +
    theme(legend.position = "bottom") -> intercept.plot
#    geom_text_repel(nudge_x=0.15, size=3.5, aes(label = round(intercept,2))) +

save_plot(paste0(PLOTDIR,"truednds_at_intercept.pdf"), intercept.plot, base_width = 8, base_height=3.5)






####################################################################################################
######################## Lineplot correlations for high and low dN/dS sites ########################
####################################################################################################

bias <- raw.bias %>%
  filter(method %in% c("FEL1", "FEL2")) %>%
  inner_join(counted) %>%
  dplyr::select(dnds, site, rep, ntaxa, bl, method, nscount.ratio, truednds) %>%
  na.omit() %>% filter(!is.infinite(nscount.ratio)) %>%
  process.raw.dnds()

bias.fel2 <- raw.bias %>%
  filter(method == "FEL2") %>%
  inner_join(counted) %>%
  dplyr::select(dn, ds, site, rep, ntaxa, bl, nscount.ratio, truedn, trueds) %>%
  na.omit() %>% filter(!is.infinite(nscount.ratio)) %>%
  process.raw.dn.ds()

nobias <- raw.nobias %>%
  filter(method %in% c("FEL1", "FEL2")) %>%
  inner_join(counted) %>%
  dplyr::select(dnds, site, rep, ntaxa, bl, method, nscount.ratio, truednds) %>%
  na.omit() %>% filter(!is.infinite(nscount.ratio)) %>%
  process.raw.dnds()
nobias %>% filter(method == "FEL2") -> nobias.fel2

bias %>%
  ggplot(aes(x = factor(ntaxa), y = meanr.dnds, group = interaction(sitetype, method), color=method, shape = sitetype)) +
  geom_line(size=0.75, alpha=0.7) + geom_point(size=2.75, alpha=0.7) +
  scale_color_manual(values = dn_ds_colors, name = "Inference Method") +
  scale_shape_manual(values=c(16,17), name = "Site N/S substitutions") +
  theme(legend.title = element_text(size=12),legend.text = element_text(size=11))-> raw
grobs <- ggplotGrob(raw)$grobs
low.high.legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]


nobias %>%
  ggplot(aes(x = factor(ntaxa), y = meanr.dnds, group = interaction(sitetype, method), color=method, shape = sitetype)) +
  geom_line(size=0.75, alpha=0.6) + geom_point(size=2.75, alpha=0.6) +
  facet_grid(~bl, scales = "free_x") +
  scale_y_continuous(limits=c(-0.1, 1), breaks = c(0.1, 0.3, 0.5, 0.7, 0.9), name ="Correlation") +
  xlab("Number of Taxa") +
  ggtitle("No codon bias") +
  background_grid("xy") +
  scale_color_manual(values = c("red", "dodgerblue"), name = "Inference Method") +
  scale_shape_manual(values=c(16,17), name = "Site N/S Ratio")+
  theme(axis.text.x = element_text(size = 10, angle=30), legend.position = "none") -> low.high.nobias

bias %>%
    ggplot(aes(x = factor(ntaxa), y = meanr.dnds, group = interaction(sitetype, method), color=method, shape = sitetype)) +
    geom_line(size=0.75, alpha=0.6) + geom_point(size=2.75, alpha=0.6) +
    facet_grid(~bl, scales = "free_x") +
    scale_y_continuous(limits=c(-0.1, 1), breaks = c(0.1, 0.3, 0.5, 0.7, 0.9), name ="Correlation") +
    xlab("Number of Taxa") +
    ggtitle("Codon bias") +
    background_grid("xy") +
    scale_color_manual(values = c("red", "dodgerblue"), name = "Inference Method") +
    scale_shape_manual(values=c(16,17), name = "Site N/S substitutions")+
    theme(axis.text.x = element_text(size = 10, angle=30), legend.position = "none") -> low.high.bias

ggdraw() +
      draw_plot(low.high.nobias, 0, .5, 0.8, .5) +
      draw_plot(low.high.bias, 0, 0, 0.8, .5) +
      draw_plot(low.high.legend, 0.78, 0.35, 0.25, 0.25) +
      draw_plot_label(c("A", "B"), c(0, 0), c(1, 0.5), size = 15) -> low.high.grid
save_plot(paste0(PLOTDIR, "low_high_lineplots.pdf"),low.high.grid, base_width=10, base_height=5)







####################################################################################################
######################### Correlations and RMSD for real tree simulations ##########################
####################################################################################################
real.dnds.sum$biastype <- factor(real.dnds.sum$biastype, levels=c("nobias", "bias"), labels=c("No codon bias", "Codon bias"))

r.dnds <- real.dnds.sum %>%
           ggplot(aes(x = dataset, fill = method, y = r)) +
           geom_boxplot(position=position_dodge(0.77), outlier.size = 0.75, size=0.4) +
           facet_wrap(~biastype) +
           scale_fill_manual(values=dn_ds_colors, name = "Inference Method") +
           xlab("Dataset") + ylab("Correlation") +
           scale_y_continuous(limits=c(0,1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.)) +
           theme(legend.position = "none")
rmsd.dnds <- real.dnds.sum %>%
              ggplot(aes(x = dataset, fill = method, y = rmsd)) +
              geom_boxplot(position=position_dodge(0.77), outlier.size = 0.75, size=0.4) +
              facet_wrap(~biastype) +
              scale_fill_manual(values=dn_ds_colors, name = "Inference Method") +
              xlab("Dataset") + ylab("RMSD") +
              scale_y_continuous(limits=c(0,2), breaks=c(0, 0.5, 1, 1.5, 2)) +
              theme(legend.position = "none")
real.plots <- plot_grid(r.dnds, rmsd.dnds, nrow=2, labels=c("A", "B"))
real.plots.withlegend <- plot_grid(real.plots, fubar.legend, nrow=2, rel_heights=c(1,0.08) )
save_plot(paste0(PLOTDIR, "real_r_rmsd.pdf"), real.plots.withlegend, base_width = 7.25, base_height = 4.75)

r.dn.ds <- real.dn.ds.sum %>%
            ggplot(aes(x = dataset, fill = parameter, y = r)) +
            geom_boxplot(outlier.size = 0.75, size=0.4) +
            scale_fill_manual(values = dn_ds_colors, name = "Parameter") +
            ylab("Correlation") + xlab("Dataset") +
            theme(legend.position = "none")
rmsd.dn.ds <- real.dn.ds.sum %>%
                ggplot(aes(x = dataset, fill = parameter, y = rmsd)) +
                geom_boxplot(outlier.size = 0.75, size=0.4) +
                scale_fill_manual(values = dn_ds_colors, name = "Parameter") +
                ylab("RMSD") + xlab("Dataset") +
                theme(legend.position = "none")
real.plots.dn.ds <- plot_grid(r.dn.ds, rmsd.dn.ds, nrow=1, labels=c("A", "B"))
real.plots.dn.ds.legend <- plot_grid(real.plots.dn.ds, dn.ds.legend, nrow=2, rel_heights=c(1,0.1) )
save_plot(paste0(PLOTDIR, "real_dn_ds_r_rmsd.pdf"), real.plots.dn.ds.legend, base_width = 9, base_height = 2.5)



####################################################################################################
############################## Divergence vs. number of taxa #######################################
####################################################################################################


ntaxa.bl.r <- treelen.dat %>%
                filter(pitype == "unequalpi") %>%
                ggplot(aes(x = as.factor(bl), y = r, fill = as.factor(biastype))) +
                geom_boxplot(position = position_dodge(0.75), lwd=0.3, outlier.size=0.8) +
                scale_x_discrete(labels = c("N = 2048\nB = 0.04","N = 512\nB = 0.16","N = 128\nB = 0.64"), name = "Parameterization") +
                scale_y_continuous(limits=c(0.45, 0.95), name = "Correlation") +
                scale_fill_manual(values = c("grey40", "grey80")) + theme(legend.position = "none")

ntaxa.bl.rmsd <- treelen.dat %>%
                  filter(pitype == "unequalpi") %>%
                  ggplot(aes(x = as.factor(bl), y = rmsd, fill = as.factor(biastype))) +
                  geom_boxplot(position = position_dodge(0.75), lwd=0.3, outlier.size=0.8) +
                  scale_x_discrete(labels = c("N = 2048\nB = 0.04","N = 512\nB = 0.16","N = 128\nB = 0.64"), name = "Parameterization") +
                  scale_y_continuous(limits=c(0, 0.75), name = "RMSD") +
                  scale_fill_manual(values = c("grey40", "grey80")) + theme(legend.position = "none")
ntaxa.bl.r.rmsd <- plot_grid(ntaxa.bl.r, ntaxa.bl.rmsd, nrow=1, labels=c("A", "B"))
ntaxa.bl.r.rmsd.legend <- plot_grid(ntaxa.bl.r.rmsd, treelen.legend, nrow=2, rel_heights=c(1,0.1))
save_plot(paste0(PLOTDIR, "ntaxa_bl_r_rmsd_unequalpi.pdf"), ntaxa.bl.r.rmsd.legend, base_width=8, base_height=3.25)


ntaxa.bl.r <- treelen.dat %>%
                filter(pitype == "equalpi") %>%
                ggplot(aes(x = as.factor(bl), y = r, fill = as.factor(biastype))) +
                geom_boxplot(position = position_dodge(0.75), lwd=0.3, outlier.size=0.8) +
                scale_x_discrete(labels = c("N = 2048\nB = 0.04","N = 512\nB = 0.16","N = 128\nB = 0.64"), name = "Parameterization") +
                scale_y_continuous(limits=c(0.45, 0.95), name = "Correlation") +
                scale_fill_manual(values = c("grey40", "grey80")) + theme(legend.position = "none")

ntaxa.bl.rmsd <- treelen.dat %>%
                  filter(pitype == "equalpi") %>%
                  ggplot(aes(x = as.factor(bl), y = rmsd, fill = as.factor(biastype))) +
                  geom_boxplot(position = position_dodge(0.75), lwd=0.3, outlier.size=0.8) +
                  scale_x_discrete(labels = c("N = 2048\nB = 0.04","N = 512\nB = 0.16","N = 128\nB = 0.64"), name = "Parameterization") +
                  scale_y_continuous(limits=c(0, 0.75), name = "RMSD") +
                  scale_fill_manual(values = c("grey40", "grey80")) + theme(legend.position = "none")
ntaxa.bl.r.rmsd <- plot_grid(ntaxa.bl.r, ntaxa.bl.rmsd, nrow=1, labels=c("A", "B"))
ntaxa.bl.r.rmsd.legend <- plot_grid(ntaxa.bl.r.rmsd, treelen.legend, nrow=2, rel_heights=c(1,0.1))
save_plot(paste0(PLOTDIR, "ntaxa_bl_r_rmsd_equalpi.pdf"), ntaxa.bl.r.rmsd.legend, base_width=8, base_height=3.25)
