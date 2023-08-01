# set seed and libraries
set.seed(2178)
library(here)
library(tidyr)
library(dplyr)

# load in data
load(here("results/tissue_specific_edges/merged_edges_object.Rdata")) # load merged_df


##We wanted to identify tissue-specific edges as previously defined in [Sonawane et al Cell Reports 2017](https://www.cell.com/cell-reports/fulltext/S2211-1247(17)31418-3#secsectitle0080); formulas 2 & 3 to calculate tissue-specific-edges
#turn any negative edge weights to zero
merged_df[merged_df < 0 ] <- 0

#remove edges that are all 0
merged_df <- merged_df %>%
  filter(rowSums(.[, -1] != 0) > 0) # skipping edges column

# make df for storing median and iqr
med_iqr_df <- data.frame(edges = merged_df$edges)

# calculate median and IQR for each edge (row) across all tissues 
med_iqr_df$median_all <- apply(merged_df[, -1], 1, median) # median for each row
med_iqr_df$IQR_all <- apply(merged_df[, -1], 1, IQR) # IQR for each row

# compare edge weights across tissues using Sonawane et al Equation 2 to calculate specificity scores
merged_df[, -1] <- t(apply(merged_df[, -1], 1, function(x) (x - med_iqr_df$median_all) / med_iqr_df$IQR_all)) #transpose because apply returns rows as rtissues and columns as edges

# calculate multiplicity
merged_df$multiplicity <- apply(merged_df[-1], 1, function(x) sum(x > 2))
multiplicity <- merged_df[,c(1,33)]

save(merged_df, file = here("results/ts_edges/specificity_multiplicity.Rdata"))

multiplicity <- multiplicity %>%
  separate(edges, into = c("TF", "gene"), sep = "_", remove = FALSE) # split edges

save(multiplicity, file = here("results/ts_edges/multiplicity.Rdata"))

sessionInfo()
