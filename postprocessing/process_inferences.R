# Process results into nice and clean dataframes, to be saved in dataframes/ directory.

require(readr)
require(dplyr)
require(tidyr)
require(purrr)
require(stringr)
source("summary_functions.R")

max_threshold = 9999     # Hyphy assigns this value (or greater, some decimal threshold lots of points out) to parameter upon failure to converge
zero_threshold = 1e-15   # My zero threshold
true.path <- "../simulation/"
DATADIR = "dataframes/"

# True dN/dS
files <- dir(path = true.path, pattern = "truednds*")
data_frame(filename = files) %>%
  mutate(file_contents = map(filename,
           ~ read_csv(file.path(true.path, .)))) %>%
  unnest() %>%
  separate(filename, c("deleteme", "pitype", "biastype1"), sep="_") %>%
  mutate(biastype = str_replace(biastype1, "(.+).csv", "\\1"),
         truedn = dn,
         trueds = ds,
         truednds = truedn/trueds) %>%
  dplyr::select(-deleteme, -biastype1, -dnds, -dn, -ds) -> truednds
truednds %>% filter(pitype == "unequalpi") %>% dplyr::select(-pitype) -> truednds.piunequal


# Counted substitutions, balanced trees
counted.path <- "../count_substitutions/balanced_trees_counts"
files <- dir(path = counted.path, pattern = "*counted.txt")
data_frame(filename = files) %>%
  mutate(file_contents = map(filename,
           ~ read_tsv(file.path(counted.path, .)))) %>%
  unnest() %>%
  separate(filename, c("rep1", "ntaxa1", "bl1", "pitype", "biastype", "deleteme"), sep="_") %>%
  mutate(rep = as.integer(str_replace(rep1, "rep([0-9]+)","\\1")),
         ntaxa = 2**(as.integer(str_replace(ntaxa1, "n([0-9]+)", "\\1"))),
         bl = as.numeric(str_replace(bl1, "bl(0.[0-9]+)", "\\1")),
         ncount = ns_changes,
         scount = s_changes) %>%
  dplyr::select(-rep1, -ntaxa1, -bl1, -deleteme, -ns_sites, -s_sites, -ns_changes, -s_changes) -> counted.dat

######### Process inferences for real tree simulations ##############
data.path <- "../results/realtrees_results/"
files <- dir(path = data.path, pattern = "*FUBAR[12].txt")
data_frame(filename = files) %>%
  mutate(file_contents = map(filename,
           ~ read_csv(file.path(data.path, .)))) %>%
  unnest() %>%
  mutate(dn = beta,
         ds = alpha,
         dnds = ifelse(dn/ds >= max_threshold | dn/ds == Inf, NA, dn/ds)) %>%
  dplyr::select(filename, dn, ds, dnds) %>%
  group_by(filename) %>%
  mutate(site = 1:n()) %>% ungroup() %>%
  separate(filename, c("rep1", "dataset", "pitype", "biastype", "method1"), sep="_") %>%
  mutate(rep = as.integer(str_replace(rep1, "rep([0-9]+)","\\1")),
         method = str_replace(method1, "(FUBAR[12]).txt", "\\1")) %>%
  dplyr::select(-rep1, -method1, -pitype) -> real.dat


######### Process and merge inferences for balanced tree simulations ###########
data.path <- "../results/balancedtrees_results/"
# FUBAR12
files <- dir(path = data.path, pattern = "*FUBAR[12].txt")
data_frame(filename = files) %>%
  mutate(file_contents = map(filename,
           ~ read_csv(file.path(data.path, .)))) %>%
  unnest() %>%
  mutate(dn = beta,
         ds = alpha,
         dnds = ifelse(dn/ds >= max_threshold | dn/ds == Inf, NA, dn/ds)) %>%
  dplyr::select(filename, dn, ds, dnds) -> fubar.raw

# SLAC12
files <- dir(path = data.path, pattern = "*SLAC.txt")
data_frame(filename = files) %>%
  mutate(file_contents = map(filename,
           ~ read_tsv(file.path(data.path, .)))) %>%
  unnest() -> slac.raw.data
  slac.raw.data %>%
  mutate(dn = dN,
         ds = dS,
         dnds = ifelse(dn/ds >= max_threshold | dn/ds == Inf, NA, dn/ds)) %>%
  dplyr::select(filename, dn, ds, dnds) %>%
  mutate(filename = str_replace(filename, "(.+)SLAC.txt", "\\1SLAC2.txt")) -> slac2.raw
  slac.raw.data %>% group_by(filename) %>%
  mutate(dn = dN,
         ds = mean(dS),
         dnds = ifelse(dn/ds >= max_threshold | dn/ds == Inf, NA, dn/ds)) %>%
  ungroup() %>%
  dplyr::select(filename, dn, ds, dnds) %>%
  mutate(filename = str_replace(filename, "(.+)SLAC.txt", "\\1SLAC1.txt")) -> slac1.raw



# FEL1
files <- dir(path = data.path, pattern = "*FEL1.txt")
data_frame(filename = files) %>%
  mutate(file_contents = map(filename,
           ~ read_csv(file.path(data.path, .)))) %>%
  unnest() %>%
  mutate(dn = `dN/dS`,
         ds = 1,
         dnds = ifelse( (LRT == 0 & `p-value` == 1 & (dn == 0 | dn == 1)) | dn >= max_threshold | dn == Inf, NA, dn)) %>%
  dplyr::select(filename, dn, ds, dnds) -> fel1.raw


# FEL2
files <- dir(path = data.path, pattern = "*FEL2.txt")
data_frame(filename = files) %>%
  mutate(file_contents = map(filename,
           ~ read_csv(file.path(data.path, .)))) %>%
  unnest() %>%
  mutate(dn = dN,
         ds = dS,
         dnds = ifelse( (LRT == 0 & `p-value` == 1) | dn >= max_threshold | ds >= max_threshold | dn == Inf | ds ==  Inf | ds <= zero_threshold, NA, dn/ds)) %>%
  dplyr::select(filename, dn, ds, dnds) -> fel2.raw


# Merge all of the balanced inferences
slac1.raw %>% rbind(slac2.raw) -> slac.raw
fel1.raw %>% rbind(fel2.raw) -> fel.raw
full.balanced.inf <- rbind(fel.raw, slac.raw)
full.balanced.inf <- rbind(full.balanced.inf, fubar.raw)



# Unpack file names into columns and merge with true dN/dS values
full.balanced.inf %>% group_by(filename) %>%
  mutate(site = 1:n()) %>% ungroup() %>%
  separate(filename, c("rep1", "ntaxa1", "bl1", "pitype", "biastype", "method1"), sep="_") %>%
  mutate(rep = as.integer(str_replace(rep1, "rep([0-9]+)","\\1")),
         ntaxa = 2**(as.integer(str_replace(ntaxa1, "n([0-9]+)", "\\1"))),
         bl = as.numeric(str_replace(bl1, "bl(0.[0-9]+)", "\\1")),
         method = str_replace(method1, "([A-Z]+[12]).txt", "\\1")) %>%
  dplyr::select(-rep1, -ntaxa1, -bl1, -method1) -> balanced.dat
balanced.dat.with.true <- left_join(balanced.dat, truednds)



########## Build summary dataframes ############
real.dat.with.true <- left_join(real.dat, truednds.piunequal)
balanced.dat.with.true <- left_join(balanced.dat, truednds)
counted.dat.with.true <- left_join(counted.dat, truednds)
counted.realdat.with.true <- left_join(counted.realdat, truednds.piunequal)

# Summary of balanced compared to true
balanced.sum.dnds <- summarize_dnds(balanced.dat.with.true)
balanced.sum.dn <- summarize_dn(balanced.dat.with.true)
balanced.sum.ds <- summarize_ds(balanced.dat.with.true)
balanced.sum.dn.ds <- rbind(balanced.sum.dn, balanced.sum.ds)


# Summary of real compared to true
real.sum.dnds <- summarize_dnds_real(real.dat.with.true)
real.sum.dn <- summarize_dn_real(real.dat.with.true)
real.sum.ds <- summarize_ds_real(real.dat.with.true)
real.sum.dn.ds <- rbind(real.sum.dn, real.sum.ds)

# Mean of summaries
mean.sum.real.true.dnds <- real.sum.dnds %>% group_by(method, dataset, biastype) %>% summarize_means()
mean.sum.real.true.dn.ds <- real.sum.dn.ds %>% group_by(method, dataset, biastype, parameter) %>% summarize_means()
mean.sum.balanced.true.dnds <- balanced.sum.dnds %>% group_by(method, ntaxa, bl, biastype, pitype) %>% summarize_means()
mean.sum.balanced.true.dn.ds <- balanced.sum.dn.ds %>% group_by(method, ntaxa, bl, biastype, pitype, parameter) %>% summarize_means()


### Save all the things!
write_csv(counted.realdat.with.true, paste0(DATADIR,"substitution_counts_realtrees.csv"))
write_csv(counted.dat.with.true, paste0(DATADIR,"substitution_counts.csv"))
write_csv(balanced.sum.dnds, paste0(DATADIR,"summary_balanced_dnds.csv"))
write_csv(balanced.sum.dn.ds, paste0(DATADIR,"summary_balanced_dn_ds.csv"))
write_csv(real.sum.dnds, paste0(DATADIR,"summary_real_dnds.csv"))
write_csv(real.sum.dn.ds, paste0(DATADIR,"summary_real_dn_ds.csv"))

write_csv(mean.sum.real.true.dnds, paste0(DATADIR,"mean_summary_real_dnds.csv"))
write_csv(mean.sum.real.true.dn.ds, paste0(DATADIR,"mean_summary_real_dn_ds.csv"))
write_csv(mean.sum.balanced.true.dnds, paste0(DATADIR,"mean_summary_balanced_dnds.csv"))
write_csv(mean.sum.balanced.true.dn.ds, paste0(DATADIR,"mean_summary_balanced_dn_ds.csv"))

write_csv(real.dat, paste0(DATADIR,"results_realtrees.csv"))
balanced.dat %>% filter(pitype == "equalpi", biastype == "nobias") %>% write_csv(paste0(DATADIR, "results_balancedtrees_nobias_equalpi.csv"))
balanced.dat %>% filter(pitype == "equalpi", biastype == "bias") %>% write_csv(paste0(DATADIR, "results_balancedtrees_bias_equalpi.csv"))
balanced.dat %>% filter(pitype == "unequalpi", biastype == "nobias") %>% write_csv(paste0(DATADIR, "results_balancedtrees_nobias_unequalpi.csv"))
balanced.dat %>% filter(pitype == "unequalpi", biastype == "bias") %>% write_csv(paste0(DATADIR, "results_balancedtrees_bias_unequalpi.csv"))
