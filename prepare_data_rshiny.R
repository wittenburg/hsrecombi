# Authors: D. Wittenburg, N. Melzer
# Date: May 6, 2022
#
# With this script, output from pipeline "hsrecombi" is prepared as input for R-Shiny app "CLARITY"
# R-Shiny app is available at https://github.com/nmelzer/CLARITY
#

library(magrittr)

nchr <- 2

path.input <- "input_rshiny"
dir.create(path.input, showWarnings = F)

path.results <- read.table("directory.tmp")[1]


### 1: table of genetic map with all markers
geneticMap <- c()
for(chr in 1:nchr){
  cat('Chr', chr, '\n')
  load(file.path(path.results, paste0('geneticpositions_chr', chr, '.RData')))

  gen.em <- pos$pos.cM
  id <- gen.em < 0
  if(sum(id, na.rm = T) > 0){
    cat('set to zero', sum(id, na.rm = T), '\n')
    gen.em[id] <- 0 # due to numerics in optimization approach (even though restrictions have been set properly)
  }
  
  load(file.path(path.results, paste0('hsphase_output_chr', chr, '.Rdata')))
  gen.hs <- c(0, cumsum(hap$probRec)) * 100
  rec.hs <- c(0, hap$probRec)
  
  dis <- c(NA, diff(pos$pos.Mb))
  geneticMap <- rbind(geneticMap, cbind(chr, 1:length(pos$pos.cM), pos$pos.Mb, pos$pos.Mb * 1e+6, gen.em, gen.hs, rec.hs, dis))
}

colnames(geneticMap) <- c('Chr', 'Name', 'Mbp_position', 'bp_position', 'cM_likelihood', 'cM_deterministic', 
                       'recrate_adjacent_deterministic', 'Mbp_inter_marker_distance')

geneticMap <- as.data.frame(geneticMap)
sum(is.na(geneticMap$cM_likelihood)) # no analysis at these SNPs because of overall homozygosity (with few exceptions)

save(geneticMap, file = file.path(path.input, "geneticMap.Rdata"))



### 2: table of genetic map summary
tab <- c()
for(Chr in 1:nchr){
  load(file.path(path.results, paste0('hsphase_output_chr', Chr, '.Rdata'))) # estimated crossover events
  tab.chr <- geneticMap[geneticMap$Chr == Chr, ]
  nSNP <- nrow(tab.chr)
  max_bp <- max(tab.chr$bp_position)
  Gap_bp <- max(diff(tab.chr$bp_position))
  Space_kb <- mean(diff(tab.chr$bp_position)) * 1e-3
  nRec <- lapply(hap$numberRec, function(z){sum(z, na.rm = T)}) %>% unlist %>% sum
  D_M <- max(tab.chr$cM_deterministic) / 100
  cMMb_D <- max(tab.chr$cM_deterministic) / (max_bp * 1e-6)
  L_M <- max(tab.chr$cM_likelihood, na.rm = T) / 100
  cMMb_L <- max(tab.chr$cM_likelihood, na.rm = T) / (max_bp * 1e-6)
  tab <- rbind(tab, c(Chr, nSNP, max_bp, Gap_bp, Space_kb, nRec, D_M, cMMb_D, L_M, cMMb_L))
}
colnames(tab) <- c("Chr", "nSNP", "max_bp", "Gap_bp", "Space_kb", "nRec", "D_M", "cMMb_D", "L_M", "cMMb_L")

end <- "#"
for(i in 2:ncol(tab)){
  if(i < 5) {
    end <- c(end, round(sum(tab[, i]), 0))
    tab[, i] <- round(tab[, i], 0)
  }
  if(i == 6 || i == 7 || i == 9){
    end <- c(end, round(sum(tab[, i]), 3))
    tab[, i] <- round(tab[, i], 3)
  }
  if(i==5 || i==8 || i==10){
    end <- c(end, round(mean(tab[, i]), 3))
    tab[, i] <- round(tab[, i], 3)
  }
}

genetic_map_summary <- rbind(tab, end)
save(genetic_map_summary, file = file.path(path.input, "genetic_map_summary.Rdata"))



### 3: table of recombination rates between adjacent SNPs for hotspot analysis 
adjacentRecRate <- data.frame(Chr = geneticMap$Chr, SNP = geneticMap$Name, cM = geneticMap$cM_deterministic, 
                              BP = geneticMap$bp_position, Theta = geneticMap$recrate_adjacent_deterministic, 
                              Dis = c(NA, diff(geneticMap$bp_position)))
save(adjacentRecRate, file = file.path(path.input, "adjacentRecRate.Rdata"))



### 4: table of misplaced markers
dirs0 <- list.dirs(recursive = F)
curl::curl_download(url ="https://www.radar-service.eu/radar-backend/archives/TKLpmZWNBmENepNd/versions/1/content",  "tmp.tar")
untar('tmp.tar')
file.remove('tmp.tar')

dat.tmp <- list.files(pattern = "MisplacedSNPs_problem_regions_paper", recursive = T)
misplacedMarkers <- readxl::read_xlsx(dat.tmp)
save(misplacedMarkers, file = file.path(path.input, "misplacedMarkers.Rdata"))

dirs1 <- list.dirs(recursive = F) 
filesstrings::dir.remove(dirs1[!(dirs1 %in% dirs0)])



### 5: table of general problematic regions
curl::curl_download(url = "https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fage.13205&file=age13205-sup-0004-TableS2.xlsx", "tmp.xlsx")
generalProblematicRegions <- readxl::read_xlsx("tmp.xlsx")
file.remove('tmp.xlsx')
save(generalProblematicRegions, file = file.path(path.input, "generalProblematicRegions.Rdata"))



### 6: table of parameters for genetic-mapping functions and data for curves
set.seed(10)
source('mapping_functions.R')

out <- c()
store2 <- list()
for (chr in 1:nchr){
  store <- list()
  
  load(file.path(path.results, paste0('geneticpositions_chr', chr, '.RData')))
  load(file.path(path.results, paste0('Results_chr', chr, '.RData')))
  
  genmap <- pos$pos.cM
  final$dist_M <- (genmap[final$SNP2] - genmap[final$SNP1]) / 100
  
  idx <- !is.na(final$dist_M)
  final <- final[idx, ]
  
  # minimization Haldane's function with scaling (Morgan -> theta)
  targetfun <- function(a) {
    mean((haldane(a * final$dist_M, inverse = T) - final$theta) ^ 2)
  }
  sln.haldane <- optim(0.1, targetfun, method = "Brent", lower = 0, upper = 10) 
  
  #  minimizing Rao's function (theta -> Morgan)
  targetfun <- function(a) {
    mean((rao(a, final$theta, inverse = F) - final$dist_M) ^ 2)
  }
  sln.rao <- optim(0.1, targetfun, method = "Brent", lower = 0, upper = 1) 
  # computationally slow if all lines included, thus:
  if(nrow(final) <= 100000) {
    idx <- 1:nrow(final)
  } else{
    idx <- sample(1:nrow(final), 100000) 
  }
  sln.rao$value <- mean((rao.inv(sln.rao$par, final$dist_M[idx]) - final$theta[idx]) ^ 2) 
  # -> MSE is now on the same scale as other map functions
  
  # minimization Felstenstein's function (Morgan -> theta)
  targetfun <- function(a) {
    mean((felsenstein(a, final$dist_M, inverse = T) - final$theta) ^ 2)
  }
  sln.felsenstein <- optim(0.1, targetfun, method = "Brent", lower = -10, upper = 10) 
  # NOTE: a < 0 phenomenon of "map expansion" which should not occur
  
  # minimization Karlin's function (Morgan -> theta)
  s <- c()
  for(N in 2:5){
    s[N] <- mean((karlin(N, final$dist_M, inverse = T) - final$theta) ^ 2)
  } 
  
  out <- rbind(out, c(chr, sln.haldane$value, sln.haldane$par, 
                      sln.rao$value, sln.rao$par, sln.felsenstein$value, sln.felsenstein$par,
                      min(s, na.rm = T), which.min(s)))
  
  r <- seq(0, max(final$dist_M), length = 101)
  y.hal <- haldane(sln.haldane$par * r, T)
  y.rao <- rao.inv(sln.rao$par, r)
  y.fel <- felsenstein(sln.felsenstein$par, r, T)
  y.kar <- karlin(which.min(s), r, T)
  
  xs <- cbind(r, r, r, r)
  ys <- cbind(y.hal, y.rao, y.fel, y.kar)
  
  red.final <- cbind(final$dist_M, final$theta)
  colnames(red.final) <- c("dist_M", "theta")
  store[[1]] <- red.final
  store[[2]] <- xs
  store[[3]] <- ys
  store2[[chr]] <- store
  save(store, file = file.path(path.input, paste0("curve-",chr,".Rdata")), compress = 'xz')
}

colnames(out) <- c('Chr', 'Haldane_scaled_mse', 'Haldane_scaled_par',
                   'Rao_mse', 'Rao_par', 'Felsenstein_mse', 'Felsenstein_par', 
                   'Karlin_mse', 'Karlin_par')

save(list = 'out', file = file.path(path.input, 'bestmapfun.RData'))
save(store2, file = file.path(path.input, "curve-all.Rdata"), compress = 'xz')

