library(hsrecombi)

args <- commandArgs(trailingOnly = TRUE)
path <- read.table("directory.tmp")[1, 1]

resfile <- args[1]
candfile <- args[2]
chr <- as.numeric(args[3])

## ----estimation of genetic map positions--------------------------------------
load(resfile)
excl <- read.table(candfile)[, 1]

out <- geneticPosition(final, exclude = excl)
pos <- list(pos.cM = c(out, rep(NA, length(map$locus_Mb) - length(out))), pos.Mb = map$locus_Mb, name = map$Name)

save(list = 'pos', file =  file.path(path, paste0('geneticpositions_chr', chr, '.RData')))

Sys.info()
