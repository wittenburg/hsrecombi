library(AlphaSimR)


## ----global variables---------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
path <- read.table("directory.tmp")[1, 2]

# Number of chromosomes
nchr <- as.numeric(args[1])
# Number of progeny in each half-sib family 
n <-  as.numeric(args[2])
# Number of SNPs (e.g., p = 1500 resembles an average bovine chromosome)
p <-  as.numeric(args[3])


## ----alphasim-----------------------------------------------------------------
founderPop <- runMacs2(nInd = 1000, nChr = nchr, segSites = p)
SP <- SimParam$new(founderPop)
SP$setSexes("yes_sys")
# Enable tracing location of recombination events
SP$setTrackRec(TRUE)  

pop <- newPop(founderPop)
N <- 10; ntotal <- N * n
my_pop <- selectCross(pop = pop, nFemale = 500, nMale = N, use = "rand", nCrosses = ntotal)

probRec <- list()
for(chr in 1:nchr){
  co.pat <- matrix(0, ncol = p, nrow = ntotal)
  for(i in 1:ntotal){
    if(nrow(SP$recHist[[1000 + i]][[chr]][[2]]) > 1){
      # 1. line contains 1 1 by default
      loci <- SP$recHist[[1000 + i]][[chr]][[2]][-1, 2] 
      co.pat[i, loci] <- 1 
    }
  }
  probRec[[chr]] <- colMeans(co.pat)
}

save(list = c('SP', 'founderPop', 'pop', 'my_pop', 'ntotal', 'probRec'), file = file.path(path, 'pop.RData'))

## ----genetic data-------------------------------------------------------------
PAT <- my_pop@father
rown <- paste(rep(unique(PAT), each = 2), c(1, 2), sep = '_')
H.pat <- pullSegSiteHaplo(pop)[rown, ]
X <- pullSegSiteGeno(my_pop)
# Physical position of markers in Mbp
map.snp <- lapply(founderPop@genMap, function(z) z * 100)

## ----plink-format-------------------------------------------------------------
map <- data.frame(Chr = rep(1:length(map.snp), unlist(lapply(map.snp, length))), 
                  Name = paste0('SNP', 1:length(unlist(map.snp))),
                  locus_Mb = unlist(map.snp), 
                  locus_bp = unlist(map.snp) * 1e+6)

colnames(X) <- map$Name
FID <- 'FAM001'
IID <- my_pop@id
MAT <- my_pop@mother
SEX <- 2
PHENOTYPE <- -9

for(chr in 1:nchr){
  write.table(map[map$Chr == chr, ], file.path(path, paste0('map_chr', chr, '.map')), 
              col.names = F, row.names = F, quote = F)
  write.table(cbind(FID, IID, PAT, MAT, SEX, PHENOTYPE, X[, map$Chr == chr]), 
              file.path(path, paste0('hsphase_input_chr', chr, '.raw')), 
              col.names = T, row.names = F, quote = F) 
}

Sys.info()
