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

sessionInfo()