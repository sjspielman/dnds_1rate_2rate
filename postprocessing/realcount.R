unequal.nobias <- counted.real.repwise %>% filter(biastype == "nobias") %>%
                    ggplot(aes(x = as.factor(dataset), y = mean_ratio)) +
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
    scale_color_hue(l=40, name = "Branch Lengths     ") +
    facet_grid(~biastype) +
    background_grid("xy") +
    scale_y_continuous(limits=c(0.2, 1.1)) +
    xlab("Number of Taxa") + ylab("dN/dS at Intercept") +
    theme(legend.position = "bottom") -> true.intercept


nsratio.intercept.grid <- plot_grid(true.ratio.scatter, true.intercept, nrow=2, labels="AUTO")
save_plot(paste0(PLOTDIR,"truednds_counts_intercept.pdf"), nsratio.intercept.grid, base_width = 7, base_height=6.5)




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
