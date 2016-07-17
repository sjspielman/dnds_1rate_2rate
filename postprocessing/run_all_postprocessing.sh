# SJS
# Post-processing full pipeline.
# NOTE: If you only want to re-create plots, you might want just run the "plot_results.R" script to save some time.

BALANCED_RESULT_DIRECTORY="" # Directory where results are stored for balanced tree inferences
REAL_RESULT_DIRECTORY="" # Directory where results are stored for empirical tree inferences
COUNTED_RESULT_DIRECTORY="" # Directory where results are stored for counted changes using ancestors


Rscript process_inferences.R      # Extracts dN/dS inferences from raw FEL, SLAC, and FUBAR inference files for simulations performed along balanced trees
#Rscript process_realtrees.R $REAL_RESULT_DIRECTORY       # Extracts dN/dS inferences from raw SLAC inference files for simulations performed along real trees
#Rscript process_counted.R $COUNTED_RESULT_DIRECTORY      # Process counted number of changes for bias simulations
#Rscript summarize_results.R                              # Summarize dN/dS, dN, and dS estimates to compute, for each quantity, correlation, estimator bias, and RMSD with true values. The resulting statistics are averaged across replicates.
Rscript build_linear_models.R                            # Builds linear models to compare performance of 1/2 rate models and inference methods
Rscript plot_figures.R                                   # Performs all plotting
