##### Repository for "One-rate models outperform two-rate models for site-specific *dN/dS* estimation" by SJS\*, SW, COW. 

##### Contents of respository:

- [`simulation/`](./simulation/) contains all scripts and data files used during sequence simulation.
   - [compute_true_dnds.py](./scripts/compute_true_dnds.py) and [compute_dnds_from_mutsel.py](./scripts/compute_dnds_from_mutsel.py) are used to calculate *dN/dS* from equilibrium frequencies.

- [`dnds_inference/`](./dnds_inference) contains all HyPhy code and batchfiles for inferring *dN/dS*.

- [`count_substitutions/`](./count_substitutions) contains code and results for substitutions counted along simulations using ancestral sequences.

- [`postprocessing/`](./postprocessing/) contains R code for statistical analyses, processed results (into CSV files), and figures. The script [`run_all_postprocessing.R`](./postprocessing/run_all_postprocessing.R) will run all steps together (although you might want to comment out certain steps for time-saving considerations, as indicated within this file).


Note that all scripts named ``*.sh`` and ``*.slurm`` were used for submitting jobs to TACC.

\*This repository is maintained by SJS. Please file any questions/comments in [Issues](https://github.com/sjspielman/sitewise_dnds_mutsel/issues/).
