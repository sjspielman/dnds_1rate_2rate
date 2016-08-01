


# Where does the switch happen?
dummy <- data.frame(x = 1:5, y = unique(counted.sitewise$bl))
ggplot(dummy, aes(x=x, y=y,color=as.factor(y))) + geom_point() + geom_line() +
scale_color_hue(l=40, name = "Branch Lengths") +theme(legend.position = "bottom")-> raw.bl.legend
grobs <- ggplotGrob(raw.bl.legend)$grobs
bl.legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

counted.sitewise %>%
  group_by(pitype, biastype, ntaxa, bl) %>%
  filter(mean_ratio > 1e-10) %>%
  do(tidy(lm(truednds ~ mean_ratio, data=.))) %>%
  select(-std.error, -statistic,-p.value) %>%
  summarize(intercept = sum(estimate)) -> dnds.nsratio.intercept

# Panel A: ntaxa=512, bl=0.04, nobias. N:S Ratio regressed on true dN/dS
counted.sitewise %>%
    filter(pitype == "unequalpi", biastype == "nobias", ntaxa == 512, bl == 0.04) %>%
    ggplot(aes(x = truednds, group=truednds, y = mean_ratio)) + geom_point() +
    geom_hline(yintercept=1, color="red") +
    scale_x_continuous(limits=c(0, 1.1), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    xlab("True dN/dS") + ylab("Mean N/S count") + ggtitle("No codon bias")-> panela

# Panel C: ntaxa=512, bl=0.04, bias. N:S Ratio regressed on true dN/dS
counted.sitewise %>%
  filter(pitype == "unequalpi", biastype == "bias", ntaxa == 512, bl == 0.04) %>%
  ggplot(aes(x = truednds, group=truednds, y = mean_ratio)) + geom_point() +
  geom_hline(yintercept=1, color="red") +
  scale_x_continuous(limits=c(0, 1.1), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  xlab("True dN/dS") + ylab("Mean N/S count") + ggtitle("Codon bias")-> panelb

g <- plot_grid(panela, panelb, nrow=2, labels="AUTO")
save_plot("nscounts_vs_truednds_subset.pdf", g, base_width = 4, base_height=6)

# Panel A: across conditions, what's the intercept? nobias
dnds.nsratio.intercept %>%
  filter(pitype == "unequalpi", biastype == "nobias") %>%
  ggplot(aes(x = as.factor(ntaxa), y = intercept, group = as.factor(bl), color = as.factor(bl))) +
    geom_line(size=0.75, alpha=0.7) + geom_point(size=2.75, alpha=0.7) +
    scale_color_hue(l=40, name = "Branch Lengths") +
    background_grid("xy") +
    scale_y_continuous(limits=c(0.2, 1.1)) +
    xlab("Number of Taxa") + ylab("dN/dS at Intercept") +
    theme(legend.position = "none") -> panela

# Panel D: across conditions, what's the intercept? nobias
dnds.nsratio.intercept %>%
  filter(pitype == "unequalpi", biastype == "bias") %>%
  ggplot(aes(x = as.factor(ntaxa), y = intercept, group = as.factor(bl), color = as.factor(bl))) +
    geom_line(size=0.75, alpha=0.7) + geom_point(size=2.75, alpha=0.7) +
    scale_color_hue(l=40, name = "Branch Lengths") +
    background_grid("xy") +
    scale_y_continuous(limits=c(0.2, 1.1)) +
    xlab("Number of Taxa") + ylab("dN/dS at Intercept") +
    theme(legend.position = "none") -> panelb

#    geom_text_repel(nudge_x=0.15, size=3.5, aes(label = round(intercept,2))) +

grid.raw <- plot_grid(panela, panelb, nrow=2, labels="AUTO")
intercept.grid <- plot_grid(grid.raw, bl.legend, nrow=2, rel_heights=c(1,0.1))
save_plot("dnds_at_intercept.pdf", intercept.grid, base_width = 5, base_height=7.5)











############### LINEPLOTS OF CORRELATIONS DEPENDING ON N/S ratio ###################


process.raw <- function(dat){
  dat %>% na.omit() %>%
  mutate(sitetype = ifelse(nscount.ratio>=1, "N/S >= 1", "N/S < 1")) %>%
  group_by(ntaxa, bl, rep, sitetype, method) %>%
  do(rraw1 = cor(.$dnds, .$truednds),
     rraw2 = cor(.$dn, .$truedn),
     rraw3 = cor(.$ds, .$trueds)) %>%
  mutate(r.dnds = rraw1[1],
         r.dn   = rraw2[1],
         r.ds   = rraw3[1]) %>%
  select(-rraw1, -rraw2, -rraw3) %>%
  ungroup() %>%
  group_by(ntaxa, bl, sitetype, method) %>%
  na.omit() %>%
  mutate(meanr.dnds = mean(r.dnds), meanr.dn = mean(r.dn), meanr.ds = mean(r.ds)) %>% select(-r.dnds, -r.dn, -r.ds) %>% unique() -> dat2
  dat2
}

process.raw.dn.ds <- function(dat){
  dat %>% na.omit() %>%
  mutate(sitetype = ifelse(nscount.ratio>=1, "N/S >= 1", "N/S < 1")) %>%
  group_by(ntaxa, bl, rep, sitetype) %>%
  do(rraw1 = cor(.$dn, .$truedn),
     rraw2 = cor(.$ds, .$trueds)) %>%
  mutate(r.dn   = rraw1[1],
         r.ds   = rraw2[1]) %>%
  select(-rraw1, -rraw2) %>%
  ungroup() %>%
  na.omit() %>% unique() %>%
  mutate(r.diff = r.dn - r.ds) -> dat2

  #gather(parameter, r, r.dn, r.ds) -> dat2
  dat2
}



raw.bias <- read_csv("dataframes/results_balancedtrees_bias_unequalpi.csv", col_types = list(
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
raw.nobias <- read_csv("dataframes/results_balancedtrees_nobias_unequalpi.csv", col_types = list(
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

bias <- raw.bias %>%
  filter(method %in% c("FEL1", "FEL2")) %>%
  inner_join(counted) %>%
  select(dnds, site, rep, ntaxa, bl, method, nscount.ratio, truednds) %>%
  na.omit() %>% filter(!is.infinite(nscount.ratio)) %>%
  process.raw()

bias.fel2 <- raw.bias %>%
  filter(method == "FEL2") %>%
  inner_join(counted) %>%
  select(dn, ds, site, rep, ntaxa, bl, nscount.ratio, truedn, trueds) %>%
  na.omit() %>% filter(!is.infinite(nscount.ratio)) %>%
  process.raw.dn.ds()

nobias <- raw.nobias %>%
  filter(method %in% c("FEL1", "FEL2")) %>%
  inner_join(counted) %>%
  select(dnds, site, rep, ntaxa, bl, method, nscount.ratio, truednds) %>%
  na.omit() %>% filter(!is.infinite(nscount.ratio)) %>%
  process.raw()

nobias %>% filter(method == "FEL2") -> nobias.fel2


bias %>%
  ggplot(aes(x = factor(ntaxa), y = meanr.dn, group = interaction(sitetype, method), color=method, shape = sitetype)) +
  geom_line(size=0.75, alpha=0.7) + geom_point(size=2.75, alpha=0.7) +
  facet_grid(~bl, scales = "free_x") +
  scale_y_continuous(limits=c(-0.1, 1), breaks = c(-0.1,0.1, 0.3, 0.5, 0.7, 0.9), name ="Correlation") +
  xlab("Number of Taxa") +
  ggtitle("Codon bias") +
  background_grid("xy") +
  scale_color_manual(values = c("red", "dodgerblue"), name = "Inference Method") +
  scale_shape_manual(values=c(16,17), name = "Site N/S Ratio")+
  theme(axis.text.x = element_text(size = 10, angle=30)) -> low.high.bias
  save_plot("dnds_lineplots_low_high_bias.pdf", low.high.bias, base_width = 10, base_height=3.5)


nobias %>%
  ggplot(aes(x = factor(ntaxa), y = meanr, group = interaction(sitetype, method), color=method, shape = sitetype)) +
  geom_line(size=0.75, alpha=0.7) + geom_point(size=2.75, alpha=0.7) +
  facet_grid(~bl, scales = "free_x") +
  scale_y_continuous(limits=c(0, 1), breaks = c(0.1, 0.3, 0.5, 0.7, 0.9), name ="Correlation") +
  xlab("Number of Taxa") +
  ggtitle("No codon bias") +
  background_grid("xy") +
  scale_color_manual(values = c("red", "dodgerblue"), name = "Inference Method") +
  scale_shape_manual(values=c(16,17), name = "Site N/S Ratio")+
  theme(axis.text.x = element_text(size = 10, angle=30)) -> low.high.nobias
save_plot("dnds_lineplots_low_high_nobias.pdf", low.high.nobias, base_width = 10, base_height=3.5)
