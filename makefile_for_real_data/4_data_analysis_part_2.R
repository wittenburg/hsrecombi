library(hsrecombi)

args <- commandArgs(trailingOnly = TRUE)
path <- read.table("directory.tmp")[1, 1]

resfile <- args[1]
candfile <- args[2]
genofile <- args[3]
chr <- as.numeric(args[4])

## ----estimation of genetic map positions--------------------------------------
load(resfile)
excl <- read.table(candfile)[, 1]

out <- geneticPosition(final, exclude = excl)
pos <- list(pos.cM = c(out, rep(NA, length(map$locus_Mb) - length(out))), pos.Mb = map$locus_Mb, name = map$Name)

save(list = 'pos', file =  file.path(path, paste0('geneticpositions_chr', chr, '.Rdata')), compress = 'xz')


## ----deterministic approach revisited with misplaced markers excluded---------
if(any(!is.na(excl))){
  # 2: Genotype matrix
  genomatrix <- data.table::fread(genofile)
  X <- as.matrix(genomatrix[, -c(1:6)])
  
  # 3: Assign daughters to sire IDs
  daughterSire <- genomatrix$PAT
  
  # 4: Estimate sire haplotypes and format data (exclude data when sire is unknown)
  hap <- makehappm(setdiff(unique(daughterSire), "0"), daughterSire, X, exclude = excl)
  save('hap', file = file.path(path, paste0('hsphase_output_chr', chr, '.Rdata')), compress = 'xz')
}

Sys.info()
