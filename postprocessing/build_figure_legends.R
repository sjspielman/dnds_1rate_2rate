
### Frequently used legends
reordered.dnds.sum.mean <- dnds.sum.mean
reordered.dnds.sum.mean$method <- factor(reordered.dnds.sum.mean$method, levels=c("SLAC1", "SLAC2", "FEL1", "FEL2", "FUBAR1", "FUBAR2"), labels=c("SLAC1  ", "SLAC2  ", "FEL1  ", "FEL2  ", "FUBAR1", "FUBAR2"))
reordered.method.colors <- c("deeppink3", "skyblue", "red", "dodgerblue", "darkred", "blue")
methods.legend.grob <- reordered.dnds.sum.mean %>% filter(pitype == "equalpi", biastype == "nobias") %>% ggplot(aes(x = factor(ntaxa), y = estbias, group = method, color = method)) + geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + scale_color_manual(values = reordered.method.colors, name = "Inference Method") + theme(legend.position = "bottom",legend.text = element_text(size=10), legend.title = element_text(size=11), legend.key.size = unit(.3, "cm"))
grobs <- ggplotGrob(methods.legend.grob)$grobs
methods.legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

reordered.dnds.sum.mean2 <- dnds.sum.mean %>% filter(method %in% c("FEL1", "FEL2", "FUBAR1", "FUBAR2"))
reordered.dnds.sum.mean2$method <- factor(reordered.dnds.sum.mean2$method, levels=c("FEL1", "FEL2", "FUBAR1", "FUBAR2"), labels=c("FEL1  ", "FEL2  ", "FUBAR1", "FUBAR2"))
reordered.method.colors2 <- c("red", "dodgerblue", "darkred", "blue")
methods.legend.grob2 <- reordered.dnds.sum.mean2 %>% filter(pitype == "equalpi", biastype == "nobias") %>% ggplot(aes(x = factor(ntaxa), y = estbias, group = method, color = method)) + geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) + scale_color_manual(values = reordered.method.colors2, name = "Inference Method") + theme(legend.position = "bottom",legend.text = element_text(size=10), legend.title = element_text(size=11), legend.key.size = unit(.3, "cm"))
grobs2 <- ggplotGrob(methods.legend.grob2)$grobs
methods.legend2 <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]


dn.ds.legend.grob <- dn.ds.sum %>% filter(pitype == "equalpi", biastype == "bias", method == "FUBAR2") %>% ggplot(aes(x = factor(ntaxa), y = r, fill = parameter)) + geom_violin(scale = "width") + scale_fill_manual(values = dn_ds_colors, name = "Parameter", labels = c("dN    ", "dS")) + theme(legend.position = "bottom", legend.key.size = unit(.3, "cm"))
grobs <- ggplotGrob(dn.ds.legend.grob)$grobs
dn.ds.legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]



fubar.legend.grob <- realdat.sum %>% ggplot(aes(x = dataset, fill = method, y = r)) + geom_boxplot(position=position_dodge(0.77), outlier.size = 0.75, size=0.4) + facet_wrap(~type) + scale_fill_manual(values=dn_ds_colors, name = "Inference Method   ", labels=c("FUBAR1 ", "FUBAR2")) + theme(legend.position = "bottom")
grobs <- ggplotGrob(fubar.legend.grob)$grobs
fubar.legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]


treelen.legend.grob <- treelen.dat %>% ggplot(aes(x = as.factor(bl), y = r, fill = as.factor(type))) + geom_boxplot(position = position_dodge(0.9)) + xlab("Parameterization") + scale_fill_manual(name = "", labels = c("No Codon bias     ", "Codon bias"), values = c("grey40", "grey80")) + theme(legend.position = "bottom")
treelen.grobs <- ggplotGrob(legend.grob)$grobs
treelen.legend <- grobs[[which(sapply(treelen.grobs, function(x) x$name) == "guide-box")]]



na.value <- "grey90"
heatmap.legend.raw <- linmodels.heatmap %>%
  ggplot(aes(x = as.factor(ntaxa), y = as.factor(bl), fill = r.diff)) +
    geom_tile() +
    scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "red", na.value = na.value, name = "R difference", limits=c(-0.05, 0.25), breaks=c(-0.05, 0, 0.05, 0.1, 0.15, 0.2, 0.25)) +
    theme(legend.position = "bottom", legend.key.width = unit(1.5, "cm"), legend.title = element_text(size = 16), legend.text = element_text(size = 14))
grobs <- ggplotGrob(heatmap.legend.raw)$grobs
heatmap.legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
