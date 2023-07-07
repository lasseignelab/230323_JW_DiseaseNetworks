# script for plotting median TPM expression of SETBP1 and its targets across tissues

# load libraries
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  source(here("src/expression_functions.R"))
})

# list all subsetted files
setbp1_files <- list.files(path = here("results/SETBP1_Expression/"), pattern = "_targets_tpm.csv", all.files = TRUE, full.names = TRUE)
#TEST SET
#setbp1_files <- setbp1_files[1:2]

# function to calculate median TPM values across samples for each file
# calc_median_tpm <- function(file){
#   tpm <- read.csv(file, row.name = 1)
#   median_tpm <- tpm %>% dplyr::select(., starts_with("GTEX")) %>% apply(., 1, median)
#   median_tpm <- tpm %>% mutate(., Median = median_tpm, .after = "Tissue")
#   
# return(median_tpm)  
# }

#test out fxn
test <- calc_median_tpm(setbp1_files[1])
test2 <- pull_median_tpm(test)
test3 <- get_gtex_median_tpm(setbp1_files[1])

# calculate median tpm of each gene for each tissue and compile list
tissues_median_tpm <- map_dfc(setbp1_files, get_gtex_median_tpm)
# remove any columns with NAs
tissues_median_tpm <- tissues_median_tpm %>% select_if(~ !any(is.na(.)))

p <- tissues_median_tpm %>% 
  pivot_longer(cols = everything(), names_to = "Tissue", values_to = "Median_TPM") %>% 
  mutate(Median_TPM = log1p(Median_TPM), .keep = "unused") %>%
  mutate(Tissue = str_extract(Tissue, ".*(?=Median_TPM)"), .keep = "unused") %>%
  ggplot(aes(x = Tissue, y = Median_TPM)) + 
  geom_violin() + 
  ylab("Median TPM Scaled (1+log)") +
  ggtitle("Expression of SETBP1 and Targets Across GTEx Tissues") + 
  theme(axis.text.x = element_text(size = 5, angle = 80, vjust = 0.60, hjust = 0.2)) 
  #scale_y_continuous(trans = "pseudo_log") + 
#add boxplot 
p <- p + geom_boxplot(width = 0.1)
ggsave("./results/testviolin.png")

png(here("results/test.png"), width = 25, height = 30, units = "cm", res = 300)
set.seed(1)
ComplexHeatmap::Heatmap(as.matrix(tissues_median_tpm), row_labels = rownames(tissues_median_tpm)) 

dev.off()