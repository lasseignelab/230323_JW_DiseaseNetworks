#!/usr/bin/env Rscript
set.seed(2178)
# Load necessary libraries
print("Loading the necessary libraries and setting seed...")
library(dplyr)
library(tidyr)
library(netZooR)


args <- commandArgs(trailingOnly = TRUE)
file1 <- basename(args[1])


cat("file1", file1, "\n")

print("Extracting tissue/sample information...")
# Initialize empty vectors to store the extracted information
gtex_tissue1 <- character(length(file1))
condition1 <- character(length(file1))
tissue1 <- character(length(file1))
  
# Get the filename without extension
file_name1 <- tools::file_path_sans_ext(basename(file1))
  
# Split the filename based on underscores
name_parts1 <- strsplit(file_name1, "_")[[1]]
  
# Extract the relevant information
gtex_tissue1 <- paste(name_parts1[1:(length(name_parts1)-2)], collapse = " ")
tissue_type1 <- paste(name_parts1[(length(name_parts1)-1)], collapse = " ") 

cat("GTEx Tissue:", gtex_tissue1, " ", "\n")
#cat("Condition:", condition1, " ", "\n")
cat("Tissue tType:", tissue_type1, " ", "\n")

# combining file name with full folder path  
folder_path <- "/data/user/jbarham3/230323_JW_DiseaseNetworks/results/PANDA"

full_path1 <- file.path(folder_path, file1)

cat("File path:", full_path1, "\n")


print("Loading in the PANDA network and formatting for CONDOR...")
# Load the Rdata file
panda_obj1 <- load(full_path1)
reg_net1 <- pandaResults@regNet
rm(pandaResults)
condor_input <- as.data.frame(reg_net1) %>% mutate(TF = rownames(reg_net1)) %>% pivot_longer(cols = -TF, names_to = "gene", values_to = "weight1")

# Print the extracted information
cat("Condor input dimensions:", dim(condor_input), "\n")

# formatting input 
condor_input[,1] <- paste(condor_input[,1],"A",sep="_")
condor_input[,2] <- paste(condor_input[,2],"B",sep="_")
  
print("Detecting communities in network...")
net_pos <- condor_input[condor_input[,3]>=0,seq_len(3)]
  
elist <- data.frame(red=net_pos[,2],blue=net_pos[,1],weights=net_pos[,3])

condor_result <- createCondorObject(elist)
condor_result <- condorCluster(condor_result,project=FALSE) # condorCluster will cluster the nodes and produce the overall modularity along with two community membership data.frames

print("Detecting modularity...")
condor_result <- condorQscore(condor_result)

print("Saving Condor object")
save(condor_result, file = paste("/data/user/jbarham3/230323_JW_DiseaseNetworks/results/condor/main/", gtex_tissue1, "_", tissue_type1, "_condor.RData", sep="")) 

print("Extracting and saving membership...")
condor_memb <- c(condor_result$red.memb[,2],condor_result$blue.memb[,2])
names(condor_memb) <- c(as.character(condor_result$red.memb[,1]),as.character(condor_result$blue.memb[,1]))
save(condor_memb, file = paste("/data/user/jbarham3/230323_JW_DiseaseNetworks/results/condor/membership/", gtex_tissue1, "_", tissue_type1, "_membership.RData", sep=""))
