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
gtex_tissue <- character(length(file1))

  
# Get the filename without extension
file_name <- tools::file_path_sans_ext(file1)
  
# Extract Tissue Name
tissue_type <- stringr::str_match(file_name, "gtex_(.*)_PANDA")[,2]

cat("GTEx Tissue:", tissue_type, " ", "\n")
#cat("Condition:", condition1, " ", "\n")


# combining file name with full folder path  
folder_path <- "/data/project/lasseigne_lab/JordanWhitlock/230323_JW_DiseaseNetworks/results/PANDA"
full_path <- file.path(folder_path, file1)

cat("File path:", full_path, "\n")


print("Loading in the PANDA network and formatting for CONDOR...")
# Load the Rdata file
load(full_path)
reg_net <- pandaResults@regNet
rm(pandaResults)
condor_input <- as.data.frame(reg_net) %>% mutate(TF = rownames(reg_net)) %>%
  pivot_longer(cols = -TF, names_to = "gene", values_to = "weight1")

# Print the extracted information
cat("Condor input dimensions:", dim(condor_input), "\n")

# formatting input 
condor_input$TF <- paste(condor_input$TF,"A",sep="_")
condor_input$gene <- paste(condor_input$gene,"B",sep="_")
  
print("Detecting communities in network...")
net_pos <- condor_input[condor_input[,3]>=0,seq_len(3)]
  
elist <- data.frame(red=net_pos$gene, blue=net_pos$TF, weights=net_pos$weight1)

condor_result <- createCondorObject(elist)
condor_result <- condorCluster(condor_result,project=FALSE) # condorCluster will cluster the nodes and produce the overall modularity along with two community membership data.frames

print("Detecting modularity...")
condor_result <- condorQscore(condor_result)

print("Saving Condor object")
save(condor_result,
     file = paste0("/data/user/clarkad/230323_JW_DiseaseNetworks/results/condor/main/",
                   "gtex_",
                   tissue_type,
                   "_condor.RData",
                   sep="")) 

print("Extracting and saving membership...")
condor_memb <- c(condor_result$red.memb[,2],condor_result$blue.memb[,2])
names(condor_memb) <- c(as.character(condor_result$red.memb[,1]),as.character(condor_result$blue.memb[,1]))
save(condor_memb,
     file = paste0("/data/user/clarkad/230323_JW_DiseaseNetworks/results/condor/membership/",
                   "gtex_",
                   tissue_type,
                   "_membership.RData",
                   sep=""))
write.csv(condor_memb,
          file = paste0("/data/user/clarkad/230323_JW_DiseaseNetworks/results/condor/membership/",
                        "gtex_",
                        tissue_type,
                        "_membership.csv",
                        sep=""))

print(paste("Analysis Complete for", gtex_tissue, tissue_type))
      
