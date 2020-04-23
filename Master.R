# This is the master file to generate all figures for the main manuscript

# All figures to be saved as pdfs in ARMN_Results/ folder

# The contour plot showing the regime for co-existence is plotted using the csv file "ARMN_Results/ps_f_PM.csv" 

source("ARMN_plotter.R") # for Figs 1, 2, 3(A,B,D,E), 5(C)

source("eigenval_analysis.R") # for Figs 3(C,F)

source("get_fmax.R") # for Figs 4, 5(A,B)