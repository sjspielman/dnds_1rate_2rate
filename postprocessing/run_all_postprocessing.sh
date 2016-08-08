# SJS
# Post-processing full pipeline.
# NOTE: If you only want to re-create plots, you might want just run the "plot_results.R" script to save some time.

Rscript process_inferences.R      # Extracts dN/dS inferences from raw FEL, SLAC, and FUBAR inference files for simulations performed along balanced trees
Rscript build_linear_models.R     # Builds linear models to compare performance of 1/2 rate models and inference methods
Rscript plot_figures.R            # Performs all plotting
