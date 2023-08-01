set.seed(2178)
library(here)
library(tidyr)
library(dplyr)

# generating edges from the regNet
files <- list.files(here("results/PANDA"), full.names = TRUE)

names <- sub("expression_PANDA.Rdata", "", basename(files)) %>% sub("gtex_", "", .) %>% sub("_PANDA.Rdata", "", .)

regNet_edges <- list()

for(i in files){
  #load in the pandaResults
  load(i)
  #extracting-networks
  regNet <- pandaResults@regNet
  regNet <- reshape2::melt(regNet, varnames = c("TF", "gene"), value.name = "edge_weight")
  
  #creating edges
  regNet$edges <- paste(regNet$TF, regNet$gene, sep = "_")
  regNet <- subset(regNet, select = -c(TF, gene))
  
  #store in list
  regNet_edges[[i]] <- regNet
  
  #print
  print(paste0(i, "melted and edges created"))
}

names(regNet_edges) <- names

save(regNet_edges, file = here("results/tissue_specific_edges/edges_list_object.Rdata"))


# combine all regNets for all tissues 
merged_df <- regNet_edges[[1]]
for (i in 2:length(regNet_edges)) {
  merged_df <- left_join(merged_df, regNet_edges[[i]], by = "edges")
}

# reorder the columns
merged_df <- merged_df %>%
  select(edges, edge_weight, everything())

# rename the columns 
names(merged_df)[-1] <- names(regNet_edges)

# save the object
save(merged_df, file = here("results/tissue_specific_edges/merged_edges_object.Rdata"))



##We wanted to identify tissue-specific edges as previously defined in [Sonawane et al Cell Reports 2017](https://www.cell.com/cell-reports/fulltext/S2211-1247(17)31418-3#secsectitle0080); formulas 2 & 3
# calculate tissue-specific-edges
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
