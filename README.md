# Repository of analyses: Symbiont coexistence may be common in mutualisms because of negative physiological feedback in preferential allocation and resource partitioning

Shyamolina Ghosh, University of Kansas 

Daniel C. Reuman, University of Kansas

James D. Bever, University of Kansas


## Introduction

This repository records the complete workflow that produced the paper. All analyses can be reproduced.

## Dependencies

   - R, R studio
   - latex, any tex compiler
   - FORTRAN 77, ifort compiler
   - Origin pro/ any other data drawing tool, e.g., Adobe illustrator (if necessary)
   
## How to compile
 
The main manuscript is written using standard microsoft word and Endnotes (for citations), Appendix can be compiled using any tex compiler (SG used Texmaker) from Appendix.tex. First, all numerical simulation from model analysis were carried out using FORTRAN 77 and ifort compiler (source code ARMN.f) and results were stored as .dat files in ARMN_Results/ARMN_dat/ folders. For convenience, we already provided the .dat files, so one can skip this step. Alternatively, for someone who is not used to FORTRAN, alternative R-script is provided as "ARMN.R". The file ARMN.R will execute the same task and can generate same data files as we stored as .dat files. The only difference is in the computation time, FORTRAN is much faster for long term simulation. Next, one needs to run the Master.R script to get all the plots to be saved as pdfs in ARMN_Results/ folder.

## Short summary for each file
   
   - **ARMN.R**: R script solving 4 variable ODE, arguments are explained in that script.
   - **ARMN.f**: equivalent FORTRAN code doing same things as the R script **ARMN.R** does, but in much faster way.
   - **ARMN.in**: input file for **ARMN.f**.
   - **ARMN_plotter.R**: plotter function related to output from **ARMN.f** and generates Figs (1 - 4) for the manuscript.
   - **eigenval_analysis.R**: R script with the function to compute eigen values at equilibrium. If maximum of real part of eigenvalues is negative, that means the equilibrium is stable.
   - **get_fmax.R**: R script that generates data plotted in Fig 5, i.e., regime of co-existence for different combination of parameters *s*, *aM*, *aN*. It makes some additional plots for range of fidelity vs. *s* given other parameters.
   - **get_MNAR_eqm_analytical.R**: R script with the function to get equilibrium of the ODE model used in ARMN.f using the analytical expression which is derived in the manuscript.
   - **LookAtStability_cleaned.r**: R script to check stability within the fidelity range (*f_min,f_max*) given all possible combination of other parameters used in the model.
   - **Appendix.tex**: supplementary latex file used for the manuscript.
   - **ARMN_alternative.R**: solver for 4 variable ODE with different functional form of *F(M,N)* used in the appendix.
   - **Master.R**: R script to generate all figures for the main manuscript, a master file.
   
## Acknowlegements 

SG and DCR are supported by a grant from the James S McDonnell Foundation, grant 1714195 from the National Science Foundation, and a grant from the California Department of Fish and Wildlife Delta Science Program. JB acknowledges support from NSF grants DEB 1556664, DEB 1738041 and OIA 1656006.
   
