# Script for subsetting GTEx TPM counts for SETBP1 and its known targets (aggregated across sources in 01_setbp1_combine_targets.R script)

# load in package libraries
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(data.table)
})


# load in the input data needed
setbp1_targets <- read.csv(file = here("data/SETBP1_Targets/setbp1_targets_geneconversions.csv"), row.names = 1)
projects <- readLines("./results/array_inputs/GTEx_exp_files_array.txt")

startTime <- Sys.time()
for(i in 1:length(projects)){
  expression <- read.csv(projects[i]) 
  expression$X <- gsub("\\.\\d+$", "", expression$X) # stripping ENSG versions
  
  # filter for expression of just the setbp1 targets
  target_expression <- setbp1_targets %>% dplyr::select(., name, target, description) %>% inner_join(., expression, by = c("target" = "X"))

  # finding tissue name
  tissue <- str_extract(projects[i], "(?<=gtex_).*?(?=_tpm)")
  
  # saving results
  if (!dir.exists("./results/SETBP1_Expression/")) dir.create("./results/SETBP1_Expression/")
  write.csv(target_expression, file = paste0("./results/SETBP1_Expression/", tissue, "_SETBP1_targets_tpm.csv"))
}
endTime <- Sys.time()
print(endTime - startTime)
