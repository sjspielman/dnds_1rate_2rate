require(dplyr)
require(tidyr)
require(grid)
require(cowplot)
require(readr)
require(broom)

source("load_data_legends.R") # Script **fully dependent on this one** to load data and make most figure legends.





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

opt.bl <- ggplot(bl, aes(y = meanbl, x = factor(ntaxa))) +
           facet_grid(biastype~bl) +
           geom_violin(alpha=0.8, scale = "width", fill="grey70") +
           xlab("Number of Taxa") + ylab("Mean optimized branch length") +
           theme(axis.text.x = element_text(angle=30, size=10))
save_plot(paste0(PLOTDIR, "optimized_bl_violins.pdf"), opt.bl, base_width = 8, base_height=4.5)

####################################################################################################
########################### Scatterplots of dN vs. dN and dN/dS vs. dN/dS  #########################
####################################################################################################

raw.bias %>% dplyr::select(-dn, -ds) %>%
  na.omit() %>%
  filter(rep == 1, method %in% c("FEL1","FEL2")) %>%
  spread(method, dnds)-> bias.compare.dnds
raw.nobias %>% dplyr::select(-dn, -ds) %>%
  na.omit() %>%
  filter(rep == 1, method %in% c("FEL1","FEL2")) %>%
  spread(method, dnds)-> nobias.compare.dnds
raw.bias %>% dplyr::select(-dnds, -ds) %>%
  na.omit() %>%
  filter(rep == 1, method %in% c("FEL1","FEL2")) %>%
  spread(method, dn)-> bias.compare.dn
raw.nobias %>% dplyr::select(-dnds, -ds) %>%
  na.omit() %>%
  filter(rep == 1, method %in% c("FEL1","FEL2")) %>%
  spread(method, dn)-> nobias.compare.dn


nobias.fel.dnds <- nobias.compare.dnds %>%
                  ggplot(aes(x = FEL1, y = FEL2)) +
                  geom_point() +
                  facet_grid(ntaxa~bl) +
                  geom_abline(slope=1, intercept=0, color="red", size=0.7) +
                  scale_x_continuous(limits=c(0,1), breaks=c(0, 0.5, 1), name = "dN/dS, FEL1") +
                  scale_y_continuous(limits=c(0,1), breaks=c(0, 0.5, 1), name = "dN/dS, FEL2") +
                  ggtitle("No codon bias") +
                  theme(axis.text = element_text(size=10.5))
bias.fel.dnds  <- bias.compare.dnds %>% filter(bl >= 0.01) %>%
                  ggplot(aes(x = FEL1, y = FEL2)) +
                  geom_point() +
                  facet_grid(ntaxa~bl) +
                  geom_abline(slope=1, intercept=0, color="red", size=0.7) +
                  scale_x_continuous(limits=c(0,1), breaks=c(0, 0.5, 1), name = "dN/dS, FEL1") +
                  scale_y_continuous(limits=c(0,1), breaks=c(0, 0.5, 1), name = "dN/dS, FEL2") +
                  ggtitle("Codon bias") +
                  theme(axis.text = element_text(size=10.5))
nobias.fel.dn <- nobias.compare.dn %>%
                  ggplot(aes(x = FEL1, y = FEL2)) +
                  geom_point() +
                  facet_grid(ntaxa~bl) +
                  geom_abline(slope=1, intercept=0, color="red", size=0.7) +
                  scale_x_continuous(limits=c(0,1), breaks=c(0, 0.5, 1), name = "dN, FEL1") +
                  scale_y_continuous(limits=c(0,1), breaks=c(0, 0.5, 1), name = "dN, FEL2") +
                  ggtitle("No codon bias") +
                  theme(axis.text = element_text(size=10.5))
bias.fel.dn   <- bias.compare.dn %>%
                  ggplot(aes(x = FEL1, y = FEL2)) +
                  geom_point() +
                  facet_grid(ntaxa~bl) +
                  geom_abline(slope=1, intercept=0, color="red", size=0.7) +
                  scale_x_continuous(limits=c(0,1), breaks=c(0, 0.5, 1), name = "dN, FEL1") +
                  scale_y_continuous(limits=c(0,1), breaks=c(0, 0.5, 1), name = "dN, FEL2") +
                  ggtitle("Codon bias") +
                  theme(axis.text = element_text(size=10.5))
dnds.dn.grid <- plot_grid(nobias.fel.dnds, bias.fel.dnds, nobias.fel.dn, bias.fel.dn, labels="AUTO", nrow=2)
save_plot(paste0(PLOTDIR, "scatter_dn_dnds_rep1.pdf"), dnds.dn.grid, base_width = 12, base_height=10)



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
                    scale_fill_hue(l=50) +
                    xlab("Number of Taxa") + ylab("Mean N/S substitutions") +
                    ggtitle("No codon bias") +
                    scale_y_continuous(limits=c(0,4)) +
                    background_grid("xy")+
                    theme(legend.position = "none")
unequal.bias <- counted.repwise %>% filter(biastype == "bias") %>%
                    ggplot(aes(x = as.factor(ntaxa), y = mean_ratio, fill = as.factor(bl))) +
                    geom_boxplot(outlier.size = 0.75) +
                    geom_hline(yintercept=1) +
                    scale_fill_hue(l=50) +
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

counted.sitewise %>% filter(pitype == "unequalpi", bl == 0.04, ntaxa == 512) -> sub.counted.sitewise
counted.sitewise %>%
  filter(mean_ratio > 1e-10, pitype == "unequalpi") %>%
  group_by(biastype, ntaxa, bl) %>%
  do(tidy(lm(truednds ~ mean_ratio, data=.))) %>%
  dplyr::select(-std.error, -statistic, -p.value) %>%
  summarize(intercept = sum(estimate)) -> dnds.nsratio.intercept
dnds.nsratio.intercept$biastype <- factor(dnds.nsratio.intercept$biastype, levels=c("nobias", "bias"), labels=c("No codon bias", "Codon bias"))
sub.counted.sitewise$biastype <- factor(sub.counted.sitewise$biastype, levels=c("nobias", "bias"), labels=c("No codon bias", "Codon bias"))

sub.counted.sitewise %>%
  ggplot(aes(x = truednds, group=truednds, y = mean_ratio)) + geom_point() +
  geom_hline(yintercept=1, color="red") +
  facet_grid(~biastype) +
  scale_y_continuous(limits=c(0,6), breaks=c(0,1,2,3,4,5,6)) +
  scale_x_continuous(limits=c(0,1.1)) +  xlab("True dN/dS") + ylab("Mean N/S substitutions") -> true.ratio.scatter

dnds.nsratio.intercept %>%
  ggplot(aes(x = as.factor(ntaxa), y = intercept, group = as.factor(bl), color = as.factor(bl))) +
    geom_line(size=0.75, alpha=0.8) + geom_point(size=2.75, alpha=0.8) +
    scale_color_hue(l=50, name = "Branch Lengths     ") +
    facet_grid(~biastype) +
    background_grid("xy") +
    scale_y_continuous(limits=c(0.2, 1.1)) +
    xlab("Number of Taxa") + ylab("dN/dS at Intercept") +
    theme(legend.position = "bottom") -> true.intercept


nsratio.intercept.grid <- plot_grid(true.ratio.scatter, true.intercept, nrow=2, labels="AUTO")
save_plot(paste0(PLOTDIR,"truednds_counts_intercept.pdf"), nsratio.intercept.grid, base_width = 7, base_height=6.5)






####################################################################################################
######################## Lineplot correlations for high and low dN/dS sites ########################
####################################################################################################

bias <- raw.bias %>%
  filter(method %in% c("FEL1", "FEL2")) %>%
  inner_join(counted) %>%
  dplyr::select(dnds, site, rep, ntaxa, bl, method, nscount.ratio, truednds) %>%
  na.omit() %>% filter(!is.infinite(nscount.ratio)) %>%
  process.raw.dnds()

nobias <- raw.nobias %>%
  filter(method %in% c("FEL1", "FEL2")) %>%
  inner_join(counted) %>%
  dplyr::select(dnds, site, rep, ntaxa, bl, method, nscount.ratio, truednds) %>%
  na.omit() %>% filter(!is.infinite(nscount.ratio)) %>%
  process.raw.dnds()

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
