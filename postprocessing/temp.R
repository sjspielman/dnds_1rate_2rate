empdnds.balanced %>%
  select(-pitype) %>%
  left_join(truednds.piunequal) %>%
  na.omit() %>%
  group_by(rep, ntaxa, bl, biastype) %>%
  do(rmsd.raw = sqrt(mean((.$empdnds - .$truednds)^2))) %>%
  mutate(rmsd = rmsd.raw[[1]]) %>% select(-rmsd.raw) -> emp.true.rmsd
emp.true.rmsd$biastype <- factor(emp.true.rmsd$biastype, levels=c("nobias", "bias"), labels=c("No codon bias", "Codon bias"))
emp.true.rmsd %>% ggplot(aes(x = as.factor(ntaxa), y = rmsd, fill = biastype)) +
                    geom_boxplot(outlier.size = 0.9) +
                    background_grid() +
                    facet_grid(~bl) +
                    scale_y_continuous(limits=c(0,8), expand=c(0,0)) +
                    scale_fill_manual(values=c("orange", "seagreen"), name ="Simulation type") +
                    scale_color_manual(values=c("orange", "seagreen")) +
                    xlab("Number of Taxa") + ylab("RMSD") +
                    theme(axis.text.x = element_text(angle=30))
