#### Repository for "A comparison of one-rate and two-rate inference frameworks for site-specific *dN/dS* estimation" by SJS\*, SW, COW.

##### Contents of repository:

- [`count_substitutions/`](./count_substitutions/) contains code and results for substitutions counted along simulations using ancestral sequences.

- [`dnds_inference/`](./dnds_inference/) contains all HyPhy code and batchfiles for inferring *dN/dS*.

- [`results/`](./results/) contains all *dN/dS* inference files made by HyPhy, for both simulations along [`balanced trees`](./results/balancedtrees_results/) and [`empirical trees`](./results/realtrees_results/). The directory [`unequalpi_nucfits/`](./results/unequalpi_nucfits/) contains results for the nucleotide fits that HyPhy performs before rate calculation (branch lengths are optimized here). The files here are specifically from SLAC inference. Note that all inference methods, as coded in HyPhy, employ more or less the same procedure for nucleotide fits (slight FUBAR exception, although not terribly different), so this data applies across all methods examined here.

- [`trees/`](./trees) contains all trees used for simulation, as well as a simple script to calculate tree statistic for empirical phylogenies. You can grab all simulated alignments using [this link](https://www.dropbox.com/sh/6k92n9h7x2vgkoq/AABbV7qGDcNy3imRQ1D80JUfa?dl=0). If that link is broken, please contact me at stephanie.spielman@gmail.com.

- [`simulation/`](./simulation/) contains all scripts and data files used during sequence simulation, as well as true simulation parameters.

- [`postprocessing/`](./postprocessing/) contains R code for statistical analyses, processed results (into CSV files), and figures. Note that processed result files are in the subdirectory [`dataframes/`](./postprocessing/dataframes/), and all figures are in the subdirectory [figures/](./postprocessing/figures/) (further organized into [`maintext/`](./postprocessing/figures/maintext/) and [`SI/`](./postprocessing/figures/SI/) for figures in the main text and supplementary info, respectively).

Note that all scripts named ``*.sh`` and ``*.slurm`` were used for submitting jobs to TACC.


\*This repository is maintained by SJS. Please file any questions/comments in [Issues](https://github.com/sjspielman/dnds_1rate_2rate/issues/), or email stephanie.spielman@gmail.com
