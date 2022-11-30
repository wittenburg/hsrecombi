library(hsrecombi)

args <- commandArgs(trailingOnly = TRUE)
path <- read.table("directory.tmp")[1, 1]

mapfile <- args[1]
genofile <- args[2]
chr <- as.numeric(args[3])


## ----data needs to be provided in plink-format--------------------------------
## ----recombination rate estimation--------------------------------------------

# 1: Physical  map
map <- read.table(mapfile, col.names = c('Chr', 'Name', 'locus_Mb', 'locus_bp'), colClasses = c('integer', 'character', 'numeric', 'integer'))
map$SNP <- 1:nrow(map)
locus_Mb <- map$locus_Mb

# 2: Genotype matrix
genomatrix <- data.table::fread(genofile)
X <- as.matrix(genomatrix[, -c(1:6)])

# 3: Assign daughters to sire IDs
daughterSire <- genomatrix$PAT

# 4: Estimate sire haplotypes and format data (exclude data when sire is unknown)
hap <- makehappm(setdiff(unique(daughterSire), "0"), daughterSire, X)
save('hap', file = file.path(path, paste0('hsphase_output_chr', chr, '.Rdata')), compress = 'xz')

# Check order and dimension
io <- sapply(1:nrow(map), function(z){grepl(x = colnames(X)[z], pattern = map$Name[z])})
if(sum(io) != nrow(map)) stop("ERROR in dimension")

# 5: Estimate recombination rates
res <- hsrecombi(hap, X, map$SNP)
final <- editraw(res, map)

if(nrow(final) == 0) message(paste('no result on chr', chr))


## ----candidates of misplaced SNPs---------------------------------------------
excl <- checkCandidates(final)


save(list = c('final', 'map'), file = file.path(path, paste0("Results_chr", chr, ".Rdata")), compress = 'xz')
write.table(excl, file = file.path(path, paste0('candidates_chr', chr, '.txt')), row.names = F, col.names = F)

Sys.info()
