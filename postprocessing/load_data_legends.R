

PLOTDIR <- "figures/"
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

# process balanced counts
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



bl <- read.csv(paste0(DATADIR, "optimized_bl.csv"))
bl$biastype <- factor(bl$biastype, levels=c("nobias", "bias"), labels=c("No codon bias", "Codon bias"))



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
methods.legend2 <- grobs[[which(sapply(grobs2, function(x) x$name) == "guide-box")]]


dn.ds.legend.grob <- real.dn.ds.sum %>% filter(biastype == "bias", method == "FUBAR2") %>% ggplot(aes(x = factor(dataset), y = r, fill = parameter)) + geom_violin(scale = "width") + scale_fill_manual(values = dn_ds_colors, name = "Parameter", labels = c("dN    ", "dS")) + theme(legend.position = "bottom", legend.key.size = unit(.3, "cm"))
grobs <- ggplotGrob(dn.ds.legend.grob)$grobs
dn.ds.legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]



fubar.legend.grob <- real.dnds.sum %>% ggplot(aes(x = dataset, fill = method, y = r)) + geom_boxplot(position=position_dodge(0.77), outlier.size = 0.75, size=0.4) + facet_wrap(~biastype) + scale_fill_manual(values=dn_ds_colors, name = "Inference Method   ", labels=c("FUBAR1 ", "FUBAR2")) + theme(legend.position = "bottom")
grobs <- ggplotGrob(fubar.legend.grob)$grobs
fubar.legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]


treelen.legend.grob <- treelen.dat %>% ggplot(aes(x = as.factor(bl), y = r, fill = as.factor(biastype))) + geom_boxplot(position = position_dodge(0.9)) + xlab("Parameterization") + scale_fill_manual(name = "", labels = c("No Codon bias     ", "Codon bias"), values = c("grey40", "grey80")) + theme(legend.position = "bottom")
grobs <- ggplotGrob(treelen.legend.grob)$grobs
treelen.legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]



na.value <- "grey90"
heatmap.legend.raw <- linmodels.heatmap %>%
  ggplot(aes(x = as.factor(ntaxa), y = as.factor(bl), fill = r.diff)) +
    geom_tile() +
    scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "red", na.value = na.value, name = "R difference", limits=c(-0.05, 0.25), breaks=c(-0.05, 0, 0.05, 0.1, 0.15, 0.2, 0.25)) +
    theme(legend.position = "bottom", legend.key.width = unit(1.5, "cm"), legend.title = element_text(size = 16), legend.text = element_text(size = 14))
grobs <- ggplotGrob(heatmap.legend.raw)$grobs
heatmap.legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]


count.grob1 <- ggplot(counted.repwise, aes(x = as.factor(ntaxa), y = mean_ratio, fill = as.factor(bl))) +
                 geom_boxplot(outlier.size = 0.75) +
                 scale_fill_hue(l=50, name = "Branch Lengths   ") +
                 theme(legend.position = "bottom")
count.grob <- ggplotGrob(count.grob1)$grobs
count.legend <- count.grob[[which(sapply(count.grob, function(x) x$name) == "guide-box")]]

count.grob1 <- ggplot(counted.repwise, aes(x = as.factor(ntaxa), y = mean_ratio, color = as.factor(bl))) +
                 geom_line(size=0.75) + geom_point(size=2) +
                 scale_color_hue(l=50, name = "Branch Lengths   ") +
                 theme(legend.position = "bottom", legend.text = element_text(size=11), legend.title = element_text(size=12))
count.grob <- ggplotGrob(count.grob1)$grobs
count.legend.ptline <- count.grob[[which(sapply(count.grob, function(x) x$name) == "guide-box")]]


# dummy <- data.frame(x = 1:5, y = unique(counted.sitewise$bl))
# ggplot(dummy, aes(x=x, y=y,color=as.factor(y))) + geom_point() + geom_line() +
# scale_color_hue(l=50, name = "Branch Lengths") +theme(legend.position = "bottom")-> raw.bl.legend
# grobs <- ggplotGrob(raw.bl.legend)$grobs
# bl.legend.counted <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]


######### ALSO INCLUDE SOME FUNCTIONS USED IN PLOTTING !!!!!!!


############### LINEPLOTS OF CORRELATIONS DEPENDING ON N/S ratio ###################


process.raw.dnds <- function(dat){
  dat %>% na.omit() %>%
  mutate(sitetype = ifelse(nscount.ratio>=1, "N/S >= 1", "N/S < 1")) %>%
  group_by(ntaxa, bl, rep, sitetype, method) %>%
  do(rraw1 = cor(.$dnds, .$truednds)) %>%
  mutate(r.dnds = rraw1[1]) %>%
  ungroup() %>%
  group_by(ntaxa, bl, sitetype, method) %>%
  na.omit() %>%
  mutate(meanr.dnds = mean(r.dnds)) %>% dplyr::select(-r.dnds) %>% unique() -> dat2
  dat2
}
#
# process.raw.dn.ds <- function(dat){
#   dat %>% na.omit() %>%
#   mutate(sitetype = ifelse(nscount.ratio>=1, "N/S >= 1", "N/S < 1")) %>%
#   group_by(ntaxa, bl, rep, sitetype) %>%
#   do(rraw1 = cor(.$dn, .$truedn),
#      rraw2 = cor(.$ds, .$trueds)) %>%
#   mutate(r.dn   = rraw1[1],
#          r.ds   = rraw2[1]) %>%
#   dplyr::select(-rraw1, -rraw2) %>%
#   ungroup() %>%
#   na.omit() %>% unique() %>%
#   mutate(r.diff = r.dn - r.ds) -> dat2
#
#   #gather(parameter, r, r.dn, r.ds) -> dat2
#   dat2
# }
