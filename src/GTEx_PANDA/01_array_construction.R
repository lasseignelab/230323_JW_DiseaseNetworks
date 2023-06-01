# NOTE: be sure to run this script to set up list of files BEFORE constructing networks with '02_PANDA_array.sh' bash script. 
# This script was run inside of JW docker container jordanwhitlock/setbp1_manuscript:1.0.5 (setbp1_manuscript_1.0.5.sif) using Singularity/3.5.2-GCC-5.4.0-2.26 
set.seed(2178)

##Load libraries
library(here)
library(utils)

##Make list of .Rdata expression files needed for array job (generated using `MouseSetbp1_PANDA_expression_05.Rmd`):
files <- list.files("/data/project/lasseigne_lab/DATASET_dir/TCGA_GTEx_CCLE_HPA_220807/gtex/tpm/individual", full.names = TRUE) # ensure that ALL exppression inputs for PANDA are in the same directory and that nothing else is in the directory.
files #save this as a text file to be the array file

write.table(files, file = here("results/array_inputs/GTEx_exp_files.txt"), sep = "\t", row.names = FALSE, col.names = FALSE) # for other data user may want to change file name

##Remove " " around each file path

# Specify the input and output file paths
input_file <- here("results/array_inputs/GTEx_exp_files.txt")
output_file <- here("results/array_inputs/GTEx_exp_files_array.txt")

# Read the input file line by line
lines <- readLines(input_file)

# Remove quotes from each line
modified_lines <- gsub("\"", "", lines)

# Write the modified lines to the output file
writeLines(modified_lines, output_file)