#!/usr/bin/env Rscript
set.seed(2178)

# Load necessary libraries
library(dplyr)
library(decoupleR)
library(CoSIA)
source("./src/functions.R")

# Fetch command line arguments
args <- commandArgs(trailingOnly = TRUE)

gtex_file <- args[1]
prior_net <- args[2]
min_n <- as.numeric(args[3])

print(gtex_file)
print(prior_net)
print(min_n)

#extracting tissue name
tissue <- basename(gtex_file)
pattern <- "^(.*?)_tpm.csv"
tissue <- sub(pattern, "\\1", tissue)
print(tissue)

# formatting expression input data
print("creating matrix for expression data")
mat <- read.csv(gtex_file)

mat <- mat[!grepl("PAR", mat$X), ] # removing PAR genes
mat$X <- gsub("\\.\\d+$", "", mat$X) # stripping ENSG versions
print("ENSG version IDs and PARs removed")

collapsed_mat <- aggregate(. ~X, data = mat, FUN = sum) # collapsing counts based off of ENSG
rownames(collapsed_mat) <- collapsed_mat$X # moving ENSG to rownames
collapsed_mat <- collapsed_mat[, -which(names(collapsed_mat) == "X")]  # dropping old ENSG column 'X'
dim(collapsed_mat)
print("Data collapsed based off of ENSG ID and ENSG moved to rownames")

## use CoSIA to convert Bgee Ensembl IDs to species gene symbols
species_rowid_vec <- rownames(collapsed_mat)
species_rowid_df <- as.data.frame(rownames(collapsed_mat)) |> dplyr::rename(X1 = `rownames(collapsed_mat)`)

converted_id <- CoSIA_ensembltosymbol("Homo_sapiens", species_rowid_df)
print("Check first converted ID")
print(head(converted_id, 1))

## match names between Ensembl id rownames and converted gene symbols
matched_names <- converted_id[,2][match(species_rowid_vec, converted_id[,1])]
print("Check first matched ID")
print(head(matched_names, 1))
collapsed_mat$HGNC <- matched_names #adding matched names back to collapsed column

## re-collapsing again based off of HGNC
collapsed_mat <- aggregate(. ~HGNC, data = collapsed_mat, FUN = sum) 
rownames(collapsed_mat) <- collapsed_mat$HGNC
collapsed_mat <- collapsed_mat[, -which(names(collapsed_mat) == "HGNC")] 
print("ENSG converted to HGNC for human expression input")

## converting to matrix
mat <- as.matrix(collapsed_mat)

# formatting CollecTRI input data
print("grabbing prior network")
net <- data.frame(read.csv(prior_net), header = TRUE)

# Define function
tf_activity <- function(mat, net, min_n){
  acts <- run_mlm(mat, net, .source = "source", .target = "target", .mor = "mor", minsize = min_n)
  return(acts)
}

# Run function and save results
print("calculating tf activity")
acts <- tf_activity(mat, net, min_n)

# Specify output path
output_path <- paste("/data/user/jbarham3/230323_JW_DiseaseNetworks/results/decoupleR/", tissue, "_acts.RData", sep = "") 
# Save the output
save(acts, file = output_path)
