# NOTE: please create a '.Rprofile' file in your project directory in the same location as your '.Rproj' file. 
# Include a single line with the following code: `R_PROFILE_USER=""`
# This will address any normalizePath warnings looking for the home directory which is not bound within the container
# This script is caled upon and run within the 02_PANDA_array.sh script

#set seed
set.seed(2178)
print("seed set")

#load in package libraries
library(netZooR)
library(data.table)
library(CoSIA)

#load in functions
source("/data/user/jbarham3/230323_JW_DiseaseNetworks/src/functions.R")
print("libraries and functions loaded")

getwd() #output wd
setwd("/data/user/jbarham3/230323_JW_DiseaseNetworks/")
getwd()
.libPaths() #output libPath

#enable args
args <- R.utils::commandArgs(trailingOnly = TRUE)

#load in the input data needed
motif <- read.table(file = "/data/user/jbarham3/230323_JW_DiseaseNetworks/data/hu_motif_all.txt", sep = "\t") #load in motif data
print("motif loaded")

ppi <- read.table(file = "/data/user/jbarham3/230323_JW_DiseaseNetworks/data/hg38_ppi.txt", sep = "\t") #load in ppi data
print("ppi loaded")

expression <- read.csv(args[1]) #load expression data from .Rdata in here from $SAMPLE_LIST
print("expression loaded")

# formatting expression data
expression <- expression[!grepl("PAR", expression$X), ] # removing PAR genes
expression$X <- gsub("\\.\\d+$", "", expression$X) # stripping ENSG versions
print("ENSG version IDs and PARs removed")

collapsed_expression <- aggregate(. ~X, data = expression, FUN = sum) # collapsing counts based off of ENSG
rownames(collapsed_expression) <- collapsed_expression$X # moving ENSG to rownames
collapsed_expression <- collapsed_expression[, -which(names(collapsed_expression) == "X")]  # dropping old ENSG column 'X'
dim(collapsed_expression)
print("Data collapsed based off of ENSG ID and ENSG moved to rownames")

# use CoSIA to convert Bgee Ensembl IDs to species gene symbols
species_rowid_vec <- rownames(collapsed_expression)
species_rowid_df <- as.data.frame(rownames(collapsed_expression)) |> dplyr::rename(X1 = `rownames(collapsed_expression)`)

converted_id <- CoSIA_ensembltosymbol("Homo_sapiens", species_rowid_df)
print("Check first converted ID")
print(head(converted_id, 1))

## match names between Ensembl id rownames and converted gene symbols
matched_names <- converted_id[,2][match(species_rowid_vec, converted_id[,1])]
print("Check first matched ID")
print(head(matched_names, 1))
collapsed_expression$HGNC <- matched_names #adding matched names back to collapsed column

## re-collapsing again based off of HGNC
collapsed_expression <- aggregate(. ~HGNC, data = collapsed_expression, FUN = sum) 
rownames(collapsed_expression) <- collapsed_expression$HGNC
collapsed_expression <- collapsed_expression[, -which(names(collapsed_expression) == "HGNC")] 
print("ENSG converted to HGNC for human expression input")

expression <- collapsed_expression #renaming for input for PANDA

#run panda on multi-omic inputs
pandaResults <- makePanda(motif, ppi, expression)
name <- sub(".Rdata", "", basename(args[1]))
save(pandaResults, file = paste0("/data/user/jbarham3/230323_JW_DiseaseNetworks/results/PANDA/", name, "_PANDA.Rdata"))
rm(pandaResults)
print(paste0(name, "_PANDA network made and saved."))