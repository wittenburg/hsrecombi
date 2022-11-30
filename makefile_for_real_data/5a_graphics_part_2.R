library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
path <- read.table("directory.tmp")[1, 1]

mapfile <- args[1]
hsoutfile <- args[2]
chr <- as.numeric(args[3])

## ----genetic map positions per chromosome-------------------------------------
load(mapfile)
data <- data.frame(Chr = rep(chr, length(pos$pos.cM)), 
                   Mbp = pos$pos.Mb, cM = pos$pos.cM)
gg2 <- ggplot(data, aes(x = Mbp, y = cM)) + geom_point(na.rm = T)

ggsave(file.path(path, paste0('genetic_map_chr', chr, '.png')), plot = gg2, device = 'png')

## ----check outcome of deterministic approach----------------------------------
load(hsoutfile)
hsphase.cM <- rep(NA, length(hap$probRec))
hsphase.cM[!is.na(hap$probRec)] <- cumsum(hap$probRec[!is.na(hap$probRec)]) * 100
hsphase.cM <- c(0, hsphase.cM)

## ----comparison of methods----------------------------------------------------
comp <- data.frame(pos = pos$pos.Mb, lik.cM = pos$pos.cM, hsphase.cM = hsphase.cM)
gg3 <- ggplot(comp) +
  geom_point(aes(x = pos, y = lik.cM, color = '1'), na.rm = T) +    
  geom_point(aes(x = pos, y = hsphase.cM, color = '2')) +
  labs(x = "physical position (Mbp)", y = "genetic position (cM)", color = "") +
  scale_color_manual(values = c("black", "grey"),
                     labels = c("likelihood-based", "deterministic"))

ggsave(file.path(path, paste0('genetic_map_all_methods_chr', chr, '.png')), plot = gg3, device = 'png')

Sys.info()
