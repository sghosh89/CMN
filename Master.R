# This is the master file to generate all figures for the main manuscript

# All figures to be saved as pdfs in ARMN_Results/ folder

# The contour plots showing the regime for co-existence are plotted using the csv files saved in "ARMN_Results/" folder 

source("ARMN_plotter.R") # for Figs 1, 2, 3, 4

source("eigenval_analysis.R") # additional plot for eigen value analysis: eigen value vs. Ps, and eigen value vs. fidelity

source("get_fmax.R") # for Fig 5: regime of co-existence for different combination of parameters s, aM, aN,
                                  # and some additional plot for range of fidelity vs. s

source("ARMN_alternative.R") # to get Fig for appendix