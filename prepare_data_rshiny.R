### Output from pipeline "hsrecombi" is prepared as input for R-Shiny app "CLARITY"
library(hsrecombi)
library(magrittr)

nchr <- 2

path.input <- "input_rshiny"
dir.create(path.input, showWarnings = F)

path.results <- read.table("directory.tmp")[1, 1]


### 1: table of genetic map with all markers
geneticMap <- list()
for(chr in 1:nchr){
  cat('Chr', chr, '\n')
  
  load(file.path(path.results, paste0('geneticpositions_chr', chr, '.Rdata')))
  dis <- c(0, diff(pos$pos.Mb))
  load(file.path(path.results, paste0('hsphase_output_chr', chr, '.Rdata')))
  
  geneticMap[[chr]] <- data.frame(Chr = chr, Name = pos$name, Mbp_position = pos$pos.Mb, 
                                  bp_position = pos$pos.Mb * 1e+6, cM_likelihood = round(pos$pos.cM, 8),
                                  cM_deterministic = round(hap$gen, 8), 
                                  recrate_adjacent_deterministic = round(c(0, hap$probRec), 8), 
                                  Mbp_inter_marker_distanc = round(dis, 8), stringsAsFactors = FALSE)
}

geneticMap <- rlist::list.rbind(geneticMap)
sum(is.na(geneticMap$cM_likelihood)) # no analysis at these SNPs because of overall homozygosity (with few exceptions)

# list of misplaced markers
geneticMap$candidate_misplacement <- 0
for(chr in 1:29) {
  candfile <- file.path(path.results, paste0('candidates_chr', chr, '_verified.txt'))
  if(file.exists(candfile)) {
    list.excl <- read.table(candfile)[, 1]
    geneticMap$candidate_misplacement[geneticMap$Chr == chr][list.excl] <- 1
  }
}

save(geneticMap, file = file.path(path.input, "geneticMap.Rdata"), compress = 'xz')
#WriteXLS::WriteXLS(geneticMap, 'GeneticMap.xlsx', na = "NA", SheetNames = Sys.Date())



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
  D_M <- max(tab.chr$cM_deterministic, na.rm = T) / 100
  cMMb_D <- max(tab.chr$cM_deterministic, na.rm = T) / (max_bp * 1e-6)
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

genetic_map_summary <- rbind(tab, end) %>% as.data.frame(., stringsAsFactors = FALSE)
save(genetic_map_summary, file = file.path(path.input, "genetic_map_summary.Rdata"), compress = 'xz')



### 3: table of recombination rates between adjacent SNPs for hotspot analysis 
adjacentRecRate <- data.frame(Chr = geneticMap$Chr, SNP = geneticMap$Name, cM = round(geneticMap$cM_deterministic, 8),
                              BP = geneticMap$bp_position, Theta = round(geneticMap$recrate_adjacent_deterministic, 8),
                              Dis = c(0, diff(geneticMap$bp_position)), stringsAsFactors = FALSE)
save(adjacentRecRate, file = file.path(path.input, "adjacentRecRate.Rdata"), compress = 'xz')



### 4 and 5: table of parameters for genetic-mapping functions and reduced data for scatter plots
set.seed(10)

out <- c()
for (chr in 1:nchr){
  store <- list()
  
  load(file.path(path.results, paste0('geneticpositions_chr', chr, '.Rdata')))
  load(file.path(path.results, paste0('Results_chr', chr, '.Rdata')))
  
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
                      min(s, na.rm = T), which.min(s))) %>% round(., 8)
  
  r <- seq(0, max(final$dist_M), length = 101)
  y.hal <- haldane(sln.haldane$par * r, T)
  y.rao <- rao.inv(sln.rao$par, r)
  y.fel <- felsenstein(sln.felsenstein$par, r, T)
  y.kar <- karlin(which.min(s), r, T)
  
  xs <- cbind(r, r, r, r)
  ys <- cbind(y.hal, y.rao, y.fel, y.kar)
  
  red.final <- cbind(final$dist_M, final$theta)
  colnames(red.final) <- c("dist_M", "theta")
  
  max.plot <- 200000
  dim.red <- nrow(red.final)
  df2 <- red.final
  
  if(dim.red >= max.plot){
    gg <- sqrt((red.final[2:dim.red, 1] - red.final[1:(dim.red - 1), 1])^2 + 
                 (red.final[2:dim.red, 2] - red.final[1:(dim.red - 1), 2])^2)
    prozent <- max.plot / dim.red
    posit <- which(gg > quantile(gg, 1 - prozent))
    df2 <- red.final[posit, ]
  }
  
  prozent <- nrow(df2) / dim.red * 100
  
  store[[1]] <- df2
  store[[2]] <- xs
  store[[3]] <- ys
  store[[4]] <- prozent 
  
  save(store, file = file.path(path.input, paste0("curve-short-", chr, ".Rdata")), compress = 'xz')
}

colnames(out) <- c('Chr', 'Haldane_scaled_mse', 'Haldane_scaled_par',
                   'Rao_mse', 'Rao_par', 'Felsenstein_mse', 'Felsenstein_par', 
                   'Karlin_mse', 'Karlin_par')

save(list = 'out', file = file.path(path.input, 'bestmapfun.Rdata'), compress = 'xz')

