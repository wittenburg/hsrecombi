AUTHORS: D. Wittenburg, N. Melzer
DATE: October 19, 2022

OBJECTIVE: 
Create breed-specific genetic map from genotypes of half-siblings and interactively visualise the outcome

SOURCE: 
The R package CLARITY is available at https://github.com/nmelzer/CLARITY
An online version of CLARITY is available at https://nmelzer.shinyapps.io/clarity/

INSTRUCTIONS:
1. Pipeline for creating the genetic map (available in subfolder makefile_for_<real/sim>_data) requires make version >= 4.3; this version allows for grouped targeting.
   Ensure that required R packages have been installed, then call "make all".
2. Run "prepare_data_rshiny.R" in the same working directory as "make all" to prepare input data for R-Shiny app CLARITY v0.2.0 in the required format. 
   The data will be content of folder CLARITY/inst/extdata.
(3.) Unless working with nchr = 29 autosomes of the cattle genome, modifications to the app become necessary.
