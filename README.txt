AUTHORS: D. Wittenburg, N. Melzer
DATE: November 30, 2022

OBJECTIVE: 
Create a breed-specific genetic map from genotypes of half-siblings and interactively visualise the
outcome

SOURCE: 
The R package CLARITY is available at https://github.com/nmelzer/CLARITY
An online version of CLARITY is available at https://nmelzer.shinyapps.io/clarity/

INSTRUCTIONS:
1. The pipeline for creating a genetic map (available in subfolder makefile_for_<real/sim>_data) 
   requires make version >= 4.3; this version allows for grouped targeting.
   Ensure that required R packages have been installed (hsrecombi, ggplot2; optional: AlphaSimR), 
   then call "make all". For each chromosome i, marker names, genetic and physical positions are 
   stored as a list in geneticpositions_chr<i>.RData. Further interim results are available.
2. To prepare input data for the R-Shiny app CLARITY version >= 1.0.0 in the required format, run 
   "prepare_data_rshiny.R" in the same working directory as "make all" (requires R packages 
   magrittr, hsrecombi). The output Rdata will then be content of folder 
   CLARITY/inst/extdata. Note that, unless working with nchr = 29 bovine autosomes, modifications
   to the app become necessary.
(3.) You may check the Rdata format by running "description_data_rshiny.Rmd" and comparing with
     https://github.com/nmelzer/CLARITY/blob/main/description_data_rshiny.html.