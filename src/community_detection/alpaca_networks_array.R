#!/usr/bin/env Rscript
set.seed(2178)
# Load necessary libraries
library(dplyr)
library(tidyr)
library(netZooR)


args <- commandArgs(trailingOnly = TRUE)
file1 <- basename(args[1])
file2 <- basename(args[2])

cat("file1", file1, "\n")
cat("file2", file2, "\n")

# Initialize empty vectors to store the extracted information
gtex_tissue1 <- character(length(file1))
condition1 <- character(length(file1))
tissue1 <- character(length(file1))

gtex_tissue2 <- character(length(file2))
condition2 <- character(length(file2))
tissue2 <- character(length(file2))
  
# Get the filename without extension
file_name1 <- tools::file_path_sans_ext(basename(file1))
file_name2 <- tools::file_path_sans_ext(basename(file2))
  
# Split the filename based on underscores
name_parts1 <- strsplit(file_name1, "_")[[1]]
name_parts2 <- strsplit(file_name2, "_")[[1]]
  
# Extract the relevant information
gtex_tissue1 <- paste(name_parts1[1:(length(name_parts1)-2)], collapse = " ")
tissue_type1 <- paste(name_parts1[(length(name_parts1)-1)], collapse = " ") 

gtex_tissue2 <- paste(name_parts2[1:(length(name_parts2)-2)], collapse = " ")
tissue_type2 <- paste(name_parts2[(length(name_parts2)-1)], collapse = " ")


cat("GTEx Tissue:", gtex_tissue1, " ", gtex_tissue2, "\n")
#cat("Condition:", condition1, " ", condition2, "\n")
cat("Tissue tType:", tissue_type1, " ", tissue_type2, "\n")

# combining file name with full folder path  
folder_path <- "/data/user/jbarham3/230323_JW_DiseaseNetworks/results/PANDA"

full_path1 <- file.path(folder_path, file1)
full_path2 <- file.path(folder_path, file2)


cat("File path:", full_path1, full_path2, "\n")
  
# Load the Rdata file
panda_obj1 <- load(full_path1)
reg_net1 <- pandaResults@regNet
rm(pandaResults)
input1 <- as.data.frame(reg_net1) %>% mutate(TF = rownames(reg_net1)) %>% pivot_longer(cols = -TF, names_to = "gene", values_to = "weight1")

panda_obj2 <- load(full_path2)
reg_net2 <- pandaResults@regNet
rm(pandaResults)
input2 <- as.data.frame(reg_net2) %>% mutate(TF = rownames(reg_net2)) %>% pivot_longer(cols = -TF, names_to = "gene", values_to = "weight2")

alpaca_input <- merge(input1, input2, by = c("TF", "gene"))

# Specify output path
alpaca_result <- alpaca(alpaca_input, file.stem = paste("/data/user/jbarham3/230323_JW_DiseaseNetworks/results/alpaca/main/", gtex_tissue1, "_", tissue_type1, ".RData", sep="")
, verbose = TRUE)
alpaca_mat <- alpaca_result[[1]] 
alpaca_membership <- as.vector(alpaca_mat)
names(alpaca_membership) <- names(alpaca_mat)
  
save(alpaca_membership, file = paste("/data/user/jbarham3/230323_JW_DiseaseNetworks/results/alpaca/membership/", gtex_tissue1, "_", tissue_type1, "_membership.RData", sep=""))


# Print the extracted information
cat("Alpaca input dimensions:", dim(alpaca_input), "\n")
