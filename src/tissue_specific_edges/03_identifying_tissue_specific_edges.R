# set seed and load libraries
set.seed(2178)
library(here)
library(tidyr)
library(dplyr)

# load in data
load(here("results/ts_edges/specificity_multiplicity.Rdata"))

# define tissue specific edges as those with a specificity scores > 2
merged_df <- merged_df[,-33]
merged_df[, -1]<- lapply(merged_df[, -1], function(x) ifelse(x > 2, 1, 0))

# melt tissue-specific edges and grab only those that are tissue-specific
tissue_specific <- melt(merged_df, id.vars = "edges", variable.name = "Tissue", value.name = "Tissue Specific")
tissue_specific <- tissue_specific[tissue_specific$`Tissue Specific` == 1, ]

# split edges
tissue_specific <- tissue_specific %>%
  separate(edges, into = c("TF", "gene"), sep = "_", remove = FALSE)

save(tissue_specific, file = here("results/ts_edges/tissue_specific_edges.Rdata"))

sessionInfo()