# Authors: P. O. Westermark, D. Wittenburg
# Date: Nov 30, 2022

# for debugging
print-%:
	@echo '$* = $($*)'

# parameters for simulation study
nchr := 2
nanimals := 1000
nsnps := 150

base_dir := .
res_dir := $(base_dir)/results
dat_dir := $(base_dir)/data


# Commands to run R, uses logfile
rbatch = Rscript --verbose $< &>> $(res_dir)/log.txt
rbatch_chr = Rscript --verbose $^ $* &>> $(res_dir)/log.txt
rbatch_sim = Rscript --verbose $< $(nchr) $(nanimals) $(nsnps) &>> $(res_dir)/log.txt

# Save all intermediate files
.SECONDARY:

# Phony targets (will always execute independently whether there are files with
# the target name in the directory. Also skips implicit rule searching)
.PHONY: all clean

# Constructs the file names of the simulated files. Automatized, since nchr could
# well be large
# NOTE, you probably want to give the nchr variable as argument to make
# (e.g., make all nchr=23)
sim_files_map := $(addprefix $(dat_dir)/map_chr, \
  $(addsuffix .map, $(shell seq 1 $(nchr))))
sim_files_raw := $(addprefix $(dat_dir)/hsphase_input_chr, \
  $(addsuffix .raw, $(shell seq 1 $(nchr))))

# &: stands for grouped target. Only works in GNU make >= 4.3. This is needed,
# in order to avoid making 2_data_simulation.R run every time one of its targets
# is a prerequisite! (it tells make that the R script generates all these files)
$(dat_dir)/pop.Rdata $(sim_files_map) $(sim_files_raw) &: 1_data_simulation.R
	mkdir -p $(res_dir)
	echo $(res_dir) $(dat_dir) > directory.tmp
	printf 'Estimation of genetic map: %s\n' "$$(date)" >$(res_dir)/log.txt	
	mkdir -p $(dat_dir)
	$(rbatch_sim)

# 2_data_analysis_part1.R accepts an argument chr with chromosome number 
# (the matched %, substituted in as $* in rbatch_chr)
$(res_dir)/hsphase_output_chr%.Rdata $(res_dir)/Results_chr%.Rdata \
$(res_dir)/candidates_chr%.txt: 2_data_analysis_part_1.R $(dat_dir)/map_chr%.map \
$(dat_dir)/hsphase_input_chr%.raw
	$(rbatch_chr)

# NOTE: necessary to check all graphic pdf's and alter candidates_chr%_verified
# manually, otherwise it contains NA for no SNP to be discarded!
$(res_dir)/candidates_chr%_verified.txt: 3_graphics_part_1.R \
$(res_dir)/Results_chr%.Rdata $(res_dir)/candidates_chr%.txt
	$(rbatch_chr)

$(res_dir)/geneticpositions_chr%.Rdata $(res_dir)/hsphase_output_chr%.Rdata: \
4_data_analysis_part_2.R $(res_dir)/Results_chr%.Rdata \
$(res_dir)/candidates_chr%_verified.txt $(dat_dir)/hsphase_input_chr%.raw
	$(rbatch_chr)

$(res_dir)/genetic_map_chr%.pdf: 5_graphics_part_2.R \
$(res_dir)/geneticpositions_chr%.Rdata $(res_dir)/hsphase_output_chr%.Rdata \
$(dat_dir)/pop.Rdata
	$(rbatch_chr)

all: $(addprefix $(res_dir)/geneticpositions_chr, \
$(addsuffix .Rdata, $(shell seq 1 $(nchr)))) \
$(addprefix $(res_dir)/genetic_map_chr, \
$(addsuffix .pdf, $(shell seq 1 $(nchr))))

clean:
	rm -rf $(res_dir) $(dat_dir) *.tmp

