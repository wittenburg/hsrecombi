AUTHORS: D. Wittenburg, N. Melzer
DATE: October 13, 2025

OBJECTIVE: 
Create a genetic map and interactively visualise the outcome

SOURCE: 
The R package CLARITY is available at https://github.com/nmelzer/CLARITY
An online version of CLARITY is available at https://nmelzer.shinyapps.io/clarity/

INSTRUCTIONS:
1. The pipeline for creating a male genetic map from genotypes of paternal half-siblings with the 
   hsrecombi approach (Hampel et al. 2018) is available in subfolder makefile_for_<real/sim>_data. 
   The pipeline requires make version >= 4.3; this version allows for grouped targeting.
   Ensure that required R packages have been installed (hsrecombi >=1.0.0, ggplot2; optional: 
   AlphaSimR), then call "make all". For each chromosome i, genetic and physical positions are 
   stored as a list in geneticpositions_chr<i>.Rdata. Further interim results are available.
2. Additionally, sex-specific and average genetic maps employing all genotypes can be estimated 
   with the LINKPHASE3 approach (Druet & George 2015). To this end, run "run_linkphase.R" in the
   same working directory as "make all" (requires R packages magrittr, hsrecombi, dplyr, data.table).
3. To prepare input data for the R-Shiny app CLARITY version >= 3.0.0 in the required format, run 
   "prepare_data_rshiny.R" in the same working directory as "make all" (requires R packages 
   magrittr, hsrecombi, dplyr). The output Rdata will then be content of folder 
   CLARITY/inst/extdata. Note that, unless working with nchr = 29 bovine autosomes, modifications
   to the app become necessary.

