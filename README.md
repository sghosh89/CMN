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
   - Origin pro/ any other data drawing tool (if necessary)
   
## How to compile
 
The main manuscript is written using standard microsoft word and Endnotes (for citations), Appendix can be compiled using any tex compiler (SG used Texmaker) from Appendix.tex. First, all numerical simulation from model analysis were carried out using FORTRAN 77 and ifort compiler (source code ARMN.f) and results were stored as .dat files in ARMN_Results/ARMN_dat/ folders. For convenience, we already provided the .dat files, so one can skip this step. Next, one needs to run the Master.R script to get all the plots to be saved as pdfs in ARMN_Results/ folder.
SG used R version 3.4.4 to run the R script.
   
## Acknowlegements 

SG and DCR are supported by a grant from the James S McDonnell Foundation, grant 1714195 from the National Science Foundation, and a grant from the California Department of Fish and Wildlife Delta Science Program. JB acknowledges support from NSF grants DEB 1556664, DEB 1738041 and OIA 1656006.
   
