This directory contains code for performing sequence simulation, as well as csv/txt files containing information about the simulation parameters. Contents:
   - [`simulate_stationary_frequencies_equalpi_nobias.py`](./simulate_stationary_frequencies_equalpi_nobias.py) simulates the base codon frequencies for simulations with equal mutation rates and no codon bias. Fitness values for codons are extracted as the log of frequencies within this file for use in simulations with unequal mutation rates and no codon bias. 

   - [`simulate_stationary_frequencies_equalpi_bias.py`](./simulate_stationary_frequencies_equalpi_bias.py) manipulates the codon frequencies produced by [`simulate_stationary_frequencies_equalpi_nobias.py`](./simulate_stationary_frequencies_equalpi_nobias.py) to incorporate codon bias. Fitness values for codons are extracted as the log of frequencies within this file for use in simulations with unequal mutation rates and codon bias. 

   - [`simulate_stationary_frequencies_unequalpi.py`](./simulate_stationary_frequencies_unequalpi.py) calculates codon frequencies for simulations with unequal mutation rates (both with and without codon bias) from the codon fitness values determined for simulations with equal mutation rates.

   - [`simulate_stationary_frequencies_equalpi_nobias.py`](./simulate_stationary_frequencies_equalpi_nobias.py) to incorporate codon bias. Fitness values for codons are extracted as the log of frequencies within this file for use in simulations with unequal mutation rates and codon bias. 

   - Simulated codon frequencies are stored in files named `codon_freq_lib_<equal/unequal>pi_<no>bias.txt`. 

   - The file [`lambda_term.csv`](./lambda_term.csv) indicates the parameter used for initial derivation of codon frequencies from a Boltzmann distribution.

   - The file [`bias_term.csv`](./bias_term.csv) indicates the extent of codon bias at each site, produced by [`simulate_stationary_frequencies_equalpi_bias.py`](./simulate_stationary_frequencies_equalpi_bias.py)
   
   - [`compute_true_dnds.py`](./compute_true_dnds.py) calculates true *dN/dS* for each parameterization. True *dN/dS* values are saved in the files named `truednds_<equal/unequal>pi_<no>bias.csv`.
   
   - [`simulate_alignments.py`](./simulate_alignments.py) performs the actual alignment simulation, using [Pyvolve](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0139047).
   
   - [`function_library.py`](./function_library.py) contains functions shared throughout scripts.
