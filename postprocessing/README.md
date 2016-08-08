This directory processes, plots, and performs statistics on all inferences using mostly R, with some Python. Dependencies include the following R packages (basically the Hadley Wickham repertoire and a few more):

+ dplyr
+ tidyr
+ cowplot (this builds on ggplot2)
+ readr
+ purrr
+ broom
+ stringr
+ grid
+ lme4
+ multcomp


To reproduce analysis and figures, run in this order...
Rscript process_inferences.R      # Extracts dN/dS inferences from raw FEL, SLAC, and FUBAR inference files for simulations performed along balanced trees
Rscript build_linear_models.R     # Builds linear models to compare performance of 1/2 rate models and inference methods
Rscript plot_figures.R            # Performs all plotting
