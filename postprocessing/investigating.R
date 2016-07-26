

raw %>% na.omit() %>%
left_join(truednds.piunequal) %>%
filter(biastype == "bias") %>%
select(-pitype, -biastype) %>%
mutate(sitetype = ifelse(site <=25, "lowdnds", "highdnds")) %>%
group_by(ntaxa, bl, rep, sitetype, method) %>%
do(rraw = cor(.$dnds, .$truednds)) %>%
mutate(r = rraw[1]) %>% select(-rraw) %>%
ungroup() %>%
group_by(ntaxa, bl, sitetype, method) %>%
na.omit() %>%
mutate(meanr = mean(r)) %>% select(-r) %>% unique()-> wompwomp
wompwomp$method <- factor(wompwomp$method, levels=c("SLAC1", "FEL1", "FUBAR1", "SLAC2", "FEL2", "FUBAR2"))
wompwomp %>%
ggplot(aes(x = factor(ntaxa), y = meanr, group = method, color=method)) +
              geom_line(size=0.75, alpha=0.7) + geom_point(size=2, alpha=0.7) +
              facet_grid(sitetype~bl, scales = "free_x") +
              scale_y_continuous(limits=c(-0.1, 1), breaks = c(-0.1,0.1, 0.3, 0.5, 0.7, 0.9), name = "Correlation") +
              xlab("Number of Taxa") +
              background_grid("xy") +
              scale_color_manual(values = method_colors, name = "Inference Method") +
              theme(axis.text.x = element_text(size = 9, angle=30)) -> low.high.dnds.lines
save_plot("dnds_lineplots_low_high.pdf", low.high.dnds.lines, base_width = 10, base_height=5)

nobias.low <- dnds.sum.mean %>%
                       filter(pitype == "unequalpi", biastype == "nobias", site >25) %>%
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











empdnds.balanced %>% select(-pitype) %>%
  left_join(truednds.piunequal) %>%
  filter(biastype == "bias") %>%
  select(-biastype) %>% na.omit() %>%
  group_by(ntaxa, bl, rep) %>%
  do(rraw1 = cor(.$empdn, .$truedn), rraw2 = cor(.$empds, .$trueds)) %>%
  mutate(r.dn = rraw1[1], r.ds = rraw2[1]) %>% select(-rraw1, -rraw2) %>%
  na.omit() %>%
  gather(rtype, r, r.dn, r.ds) %>%
  ggplot(aes(x = as.factor(ntaxa), y = r, fill=rtype))+
  geom_boxplot() + facet_grid(~bl) + coord_cartesian(ylim=c(-0.5, 1))


p <- counted %>% group_by(pitype, biastype, site, ntaxa, bl) %>% 
  mutate(mean_ratio = mean(ncount)/mean(scount)) %>%
  filter(pitype == "unequalpi", biastype == "bias") %>%
  ggplot(aes(x = site, group=site, y = mean_ratio)) + geom_point() + facet_grid(ntaxa~bl)+ geom_hline(yintercept=1, color="red") + xlab("Site (low to high dN/dS)") + ylab("Ratio of N to S substitutions counted") + scale_y_log10(limits=c(1e-3,10))
save_plot("N_over_S_counts.pdf", p, base_width = 10, base_height=5)
raw<- read_csv("dataframes/results_balancedtrees_bias_unequalpi.csv", col_types = list(
  pit
  ype = col_character(),
  biastype = col_character(),
  dn = col_double(),
  ds = col_double(),
  dnds = col_double(),
  site = col_integer(),
  rep = col_integer(),
  ntaxa = col_factor(c(128, 256, 512, 1024, 2048)),
  bl = col_factor(c(0.0025, 0.01, 0.04, 0.16, 0.64)),
  method = col_factor(c("FEL1", "FEL2", "SLAC1", "SLAC2", "FUBAR1", "FUBAR2"))))

raw %>% left_join(truednds.piunequal) %>% filter(biastype == "bias") %>% select(-pitype, -biastype, -dn, -truedn, -dnds, -truednds) %>%
  group_by(site, ntaxa, bl) -> raw2
raw2 %>% do(rmsd.raw = sqrt(mean(.$ds - .$trueds)^2)) %>% mutate(rmsd = rmsd.raw[1]) %>%
na.omit() %>% filter(!is.infinite(rmsd)) %>%
ggplot(aes(x = site, group = site, y = rmsd))+geom_point() + facet_grid(ntaxa~bl)+
coord_cartesian(ylim=c(0,5))
