### Run of LINKPHASE3
### Output from linkphase is prepared as input for R-Shiny app "CLARITY"
library(magrittr)
library(dplyr)
library(data.table)

gw <- getwd()
path.link <- "linkphase_tmp"
dir.create(path.link, showWarnings = F)

path.results <- read.table("directory.tmp")[1, 1]
path.data <- read.table("directory.tmp")[1, 2]

# genotype data
ls <- list.files(pattern = "chr[0-9]+.raw", path = path.data)
# extract chromosome number within file name
chroms <- as.numeric(gsub("\\D", "", ls))

# full path to LINKPHASE3 executable
file.exe <- system('find ~+ . -type f -name "LINKPHASE3"', intern = T)[1]
if(!file.exists(file.exe)){
  system('git clone https://github.com/tdruet/LINKPHASE3/')
  setwd('LINKPHASE3')
  system('gfortran LINKPHASE3.f90 -o LINKPHASE3')
  file.exe <- paste0(getwd(), '/LINKPHASE3')
  setwd('..')
}

# write content of parameter file linkin.txt
param <- c('#PEDIGREE_FILE', 'input_ped', '#GENOTYPE_FILE', 'input_geno', '#MARKER_FILE', 'input_map', '#HALFSIB_PHASING',
           'yes', '#HMM_PHASING', 'yes', '#N_TEMPLATES', '50', '#CHECK_PREPHASING', 'yes', '#ITERATIONS', '5')
write.table(param, file.path(path.link, 'linkin.txt'), col.names = F, row.names = F, quote = F)


# LINKPHASE 1. run
nrec <- c()
for(chr in chroms){
  # write input files for LINKPHASE3
  mapfile <- list.files(path = path.data, pattern = paste0('chr', chr, '.map'), full.names = T)
  system(paste0('awk \'{ print NR, $2, $3 }\' ', mapfile, ' > ', path.link, '/input_map'))
  genofile <- list.files(path = path.data, pattern = paste0('chr', chr, '.raw'), full.names = T)
  system(paste0('awk \'FNR>1 {print $2, $3, $4}\' ', genofile, ' > ', path.link, '/input_ped'))
  system(paste0('awk -f ../split_genotypes.awk ', genofile, ' > ', path.link, '/input_geno'))  
  
  setwd(path.link)
  system(file.exe)
  file.rename('MCS.txt', paste0('MCS_chr', chr, '.txt'))
  dat <-  read.table('nrec_hmm.txt', col.names = c('offspring', 'parent', 'parent.sex',  'prog.number', 'parent.phase.info', 'mate.geno.info', 'p.het', 'p.hom', 'p.inf' ,'co.count', 'unused'))
  tmp <- dat %>% select(offspring, parent, co.count) %>% group_by(offspring) %>% summarise(nrec_chr = sum(co.count))
  if(chr == 1) nrec <- tmp else nrec <- merge(nrec, tmp, by = 'offspring', all = T)
  setwd(gw)
}  

# filter SNPs with MCS > 0.986 & animals with nrec <= 2 * nchr -> second run and process results 
nrec <- nrec %>% mutate(total.nco = rowSums(.[-1]))
animal.ex <- nrec$offspring[nrec$total.nco > 2 * max(chroms)]

# LINKPHASE 2. run
for(chr in chroms){  
  # physical map
  mapfile <- list.files(path = path.data, pattern = paste0('chr', chr, '.map'), full.names = T)
  map <- read.table(file = mapfile, 
                    col.names = c('Chr', 'SNP', 'locus_Mb', 'locus_bp'), 
                    colClasses = c('integer', 'character', 'numeric', 'integer'))
 
  # filtering for run 2 
  mcs <- read.table(file.path(path.link, paste0('MCS_chr', chr, '.txt')), 
                    col.names = c('ID', 'Pos_M1', 'Pos_M2', 'recrate_adjacent_average', 'n.parents', 'error.rates', 'n.progeny', 'entropy', 'mcs'))
  marker.ex.id <- mcs$ID[mcs$mcs < 0.986]
  marker.ex <- map$SNP[marker.ex.id]
  
  genofile <- list.files(path = path.data, pattern = paste0('chr', chr, '.raw'), full.names = T)
  geno <- fread(genofile) %>% 
    filter(!(IID %in% animal.ex)) %>% 
    select(!matches(paste0(marker.ex, '$')))
  
  setwd(path.link)
  write.table(geno, 'geno.tmp', col.names = T, row.names = F, quote = F)
  mapR <- map[-marker.ex.id, ]
  write.table(cbind(1:nrow(mapR), mapR$SNP, mapR$locus_Mb), 'input_map', col.names = F, row.names = F, quote = F)
  system('awk \'FNR>1 {print $2, $3, $4}\' geno.tmp > input_ped')
  system(paste('awk -f', file.path(gw, '../split_genotypes.awk'), 'geno.tmp > input_geno'))
 
  system(file.exe)
  setwd(gw)
  
  ## processing the final outcome
  # recombination rates between neighbours and linkage map
  linkmap <- read.table(file.path(path.link, file = 'emap.txt'), 
                        col.names = c('ID1', 'Loc1_Mb100', 'Loc2_Mb100', 'Hm_recrate_adjacent', 'Hf_recrate_adjacent'))
  # MCS
  mcs <- read.table(file.path(path.link, file = 'MCS.txt'), 
                    col.names = c('ID', 'Pos_M1', 'Pos_M2', 'recrate_adjacent_average', 'n.parents', 'error.rates', 'n.progeny', 'entropy', 'mcs'))
  
  # emap: here RR between marker and next -> shift by -1 index to get RR between marker and previous as in other approaches
  linkmap <- rbind(c(0, NA, NA, 0, 0), linkmap)
  linkmap$ID <- linkmap$ID + 1
  linkmap$Name <- map$SNP[linkmap$ID]
  linkmap$bp <- map$locus_bp[linkmap$ID]
  linkmap$MCS <- mcs$mcs
  linkmap$Ha_recrate_adjacent <- c(0, mcs$recrate_adjacent_average[-nrow(mcs)]) # recreate between snp and previous snp
  
  linkmap <- linkmap %>% 
    mutate(Chr = chr, Hm_cM = cumsum(Hm_recrate_adjacent) * 100, 
           Hf_cM = cumsum(Hf_recrate_adjacent) * 100, 
           Ha_cM = cumsum(Ha_recrate_adjacent) * 100) %>% 
    select(Chr, Name, bp, Hm_cM, Hm_recrate_adjacent, Hf_cM, Hf_recrate_adjacent, Ha_cM, Ha_recrate_adjacent, MCS) %>% 
    mutate(across(where(is.numeric), \(x) round(x, 8)))
  
  # crossover counts
  dat <-  read.table(file.path(path.link, file = 'nrec_hmm.txt'), 
                     col.names = c('offspring', 'parent', 'parent.sex',  'prog.number', 'parent.phase.info', 'mate.geno.info', 'p.het', 'p.hom', 'p.inf' ,'co.count', 'unused'))
  dat1 <-  read.table(file.path(path.link, file = 'nrec.txt'), 
                      col.names = c('offspring', 'parent', 'parent.sex',  'prog.number', 'parent.phase.info', 'mate.geno.info', 'p.het', 'p.hom', 'p.inf' ,'co.count'))
  nrec <- list(male = dat[dat$parent.sex == 1, c('parent', 'co.count')], female = dat[dat$parent.sex == 2, c('parent', 'co.count')],
               average = dat1[, c('parent', 'co.count')])
  # all data
  save(list = c('linkmap', 'nrec'), file = paste0(path.results, '/linkphase_output_chr', chr, '.Rdata'), compress = 'xz')
}

