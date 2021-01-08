library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
path <- read.table("directory.tmp")[1, 1]

resfile <- args[1]
candfile <- args[2]
chr <- as.numeric(args[3])

## ----visual verification of preselected SNPs; keep candidates with unusually--
## ----high recombination rate to the succeeding SNPs---------------------------

load(resfile)
excl <- read.table(candfile)[, 1]

for(cand in excl){
## ----Heatmap plot of recombination rates for visual verification--------------

  win <- cand + (-100:100)
  win <- win[(win >= 1) & (win <= max(final$SNP2))]
  
  target <- final[(final$SNP1 %in% win) & (final$SNP2 %in% win), ]
  
  gg1 <- ggplot(data = target, aes(SNP2, SNP1, fill = theta)) + 
    geom_tile() +
    xlab("Locus 2") +
    ylab("Locus 1") +
    coord_equal() + 
    scale_y_continuous(trans = "reverse") +
    theme(panel.background = element_blank(), 
          panel.grid.major = element_line(colour = "grey", size = 0.1),
          panel.grid.minor = element_line(colour = "grey")) +
    theme(text = element_text(size = 18)) +
    scale_fill_gradientn(colours = c('yellow', 'red'), limits = c(0, 1+1e-10), na.value = 'white')
# -> conspicious only if band of increased values in orange

  ggsave(file.path(path, paste0('chr_', chr, '_cand_', cand, '.pdf')), gg1)
}

## ----keep index of SNPs which have been inspected and confirmed---------------
## ----you may alter the "_verified.txt" file manually afterwards---------------
excl <- NA
write.table(excl, file = file.path(path, paste0('candidates_chr', chr, '_verified.txt')), 
            row.names = F, col.names = F)

Sys.info()
