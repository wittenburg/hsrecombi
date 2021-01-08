# Authors: P. O. Westermark, D. Wittenburg
# Date: Jan 8, 2021

# for debugging
print-%:
	@echo '$* = $($*)'

base_dir := .
# directory for results 
res_dir := $(base_dir)/results

# data are stored in "dat_dir"
# genotype data files (in PLINK format) are named as sample_data_chr<chr>.raw
# map data files are denoted as map_chr<chr>.map
dat_dir := $(base_dir)/data

# number of chromosomes
nchr := 2

# Commands to run R, uses logfile
rbatch = Rscript --verbose $< &>> $(res_dir)/log.txt
rbatch_chr = Rscript --verbose $^ $* &>> $(res_dir)/log.txt

# Save all intermediate files
.SECONDARY:

# Phony targets (will always execute independently whether there are files with
# the target name in the directory. Also skips implicit rule searching)
.PHONY: all clean

# 2_data_analysis_part1.R accepts an argument chr with chromosome number 
# (the matched %, substituted in as $* in rbatch_chr)
$(res_dir)/hsphase_output_chr%.RData $(res_dir)/Results_chr%.RData \
$(res_dir)/candidates_chr%.txt: 2_data_analysis_part_1.R $(dat_dir)/map_chr%.map \
$(dat_dir)/sample_data_chr%.raw 
	mkdir -p $(res_dir)
	echo $(res_dir) $(dat_dir) > directory.tmp
	printf 'Estimation of genetic map: %s\n' "$$(date)" >$(res_dir)/log.txt	
	$(rbatch_chr)

# NOTE: necessary to check all graphic pdf's and alter candidates_chr%_verified
# manually, otherwise it contains NA for no SNP to be discarded!
$(res_dir)/candidates_chr%_verified.txt: 3_graphics_part_1.R \
$(res_dir)/Results_chr%.RData $(res_dir)/candidates_chr%.txt
	$(rbatch_chr)

$(res_dir)/geneticpositions_chr%.RData: 4_data_analysis_part_2.R \
$(res_dir)/Results_chr%.RData $(res_dir)/candidates_chr%_verified.txt
	$(rbatch_chr)

$(res_dir)/genetic_map_chr%.pdf: 5a_graphics_part_2.R \
$(res_dir)/geneticpositions_chr%.RData $(res_dir)/hsphase_output_chr%.RData
	$(rbatch_chr)

all: $(addprefix $(res_dir)/geneticpositions_chr, \
$(addsuffix .RData, $(shell seq 1 $(nchr)))) \
$(addprefix $(res_dir)/genetic_map_chr, \
$(addsuffix .pdf, $(shell seq 1 $(nchr))))

clean:
	rm -rf $(res_dir) *.tmp

