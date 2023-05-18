# Funtions for 01_EdgeSpecificity.Rmd

calculate_proportions <- function (gene_list,
                          data_table,
                          tissue_types) {
  # Subset data table for genes of interest
  gene_subsets <- list()
  for (i in 1:length(gene_list)) {
    gene_subset <- data_table[
      data_table$TargetGene == gene_list[[i]], 
    ]
    gene <- gene_list[[i]]
    gene_subsets[[gene]] <- gene_subset
  }
  gene_edge_counts <- list()
  # Calculate total number of edges for each gene of interest
  for (i in 1:length(names(gene_subsets))) {
    gene <- names(gene_subsets)[[i]]
    gene_edge_counts[[gene]] <- nrow(gene_subsets[[i]])
  }
  # Create a list of dataframes for each gene that is appropriately sized
  gene_tissue_dfs <- list()
  for (i in 1:length(names(gene_subsets))) {
    gene <- names(gene_subsets)[[i]]
    gene_tissue_dfs[[gene]] <- data.frame(
      matrix(
        ncol = 38, 
        nrow = gene_edge_counts[[gene]]
      )
    )
  }
  # Format the dataframes with the proper row and col names
  for (i in 1:length(names(gene_tissue_dfs))) {
    gene <- names(gene_tissue_dfs)[[i]]
    rownames(gene_tissue_dfs[[gene]]) <- gene_subsets[[gene]]$TF
    colnames(gene_tissue_dfs[[gene]]) <- tissue_types
  }
  # Assigne 1 for each edge that is specific per tissue
  for (i in 1:length(names(gene_tissue_dfs))) {
    gene <- names(gene_tissue_dfs)[[i]]
    for (j in 1:nrow(gene_tissue_dfs[[gene]])) {
      present <- unlist(strsplit(gene_subsets[[gene]]$Tissues[j], ","))
      for (k in 1:length(present)) {
        tissueForEdge <- present[k]
        gene_tissue_dfs[[gene]][j,which(colnames(gene_tissue_dfs[[gene]])==tissueForEdge)]<-1
      }
    }
  }
  
  # Assign 0 for TF-pkd gene pairs that are not present in each tissue
  for (i in 1:length(names(gene_tissue_dfs))) {
    gene_dataframe <- gene_tissue_dfs[[i]]
    gene_dataframe[is.na(gene_dataframe)] <- 0
    gene_tissue_dfs[[i]] <- gene_dataframe
  }
  # Sum the number of specific edges for each gene in each tissue
  edges_gene_tissues <- list()
  for (i in 1:length(names(gene_tissue_dfs))) {
    gene_dataframe <- gene_tissue_dfs[[i]]
    edges_gene_tissues[[i]] <- colSums(gene_dataframe)
  }
  # Calculate proportion of specific edges vs total edges per tissue
  proportions_per_gene <- list()
  for (i in 1:length(names(gene_tissue_dfs))) {
    gene <- names(gene_tissue_dfs)[[i]]
    gene_dataframe <- gene_tissue_dfs[[i]]
    gene_proportion <- edges_gene_tissues[[i]]/gene_edge_counts[[i]]
    proportions_per_gene[[gene]] <- gene_proportion
  }
  return(proportions_per_gene)
}


# Plotting proportions
plot_proportions <- function (proportions_df,
                              gene_list) {
  for (i in 1:length(names(formatted_proportions))) {
    gene <- names(formatted_proportions)[[i]]
    plot <- ggplot(
      formatted_proportions[[gene]],
      aes(
        x = Proportions,
        y = Tissues,
      )
    ) +
      geom_point(
        col = "tomato2",
        size = 3
      ) +
      xlim(0, 1) +
      labs(
        title = paste0(
          "Proportions of Specific Edges Per Tissue - ",
          gene_list[[i]]),
      )
    print(plot)
  }
}

permutation_test <- function(data_table,
                             tissue_types,
                             gene_oi,
                             nperm) {
  counter <- 0
  results <- data.frame(
    matrix(
      ncol = 1, 
      nrow = length(tissue_types)
    )
  )
  
  rownames(results) <- tissue_types
  colnames(results) <- gene_oi
  unique_targets <- unique(STable1_GTEx_TSInfo_Edges$TargetGene)
  gene_list <- c(
    gene_oi, 
    sample(
      unique_targets, 
      size = 999, 
      replace = FALSE
    )
  )
  print(length(gene_list))
  # Subset the dataframe based on the sampled gene list.
  permutation <- calculate_proportions(
    gene_list = gene_list,
    data_table = data_table,
    tissue_types = tissue_types
  )
  for (i in 1:length(tissue_types)) {
    for (j in 2:length(gene_list)){
      if (permutation[[1]][[tissue_types[i]]] <= permutation[[j]][[tissue_types[i]]]) {
        counter <- counter + 1
      }
    }
    p_val <- ((counter+1)/(nperm+1))
    results[i, 1] <- p_val
    counter <- 0
  }
  return(results)
}




# 
# formatted_proportions <- list()
# for (i in 1:length(names(proportions_df))) {
#   gene <- names(proportions_df)[[i]]
#   gene_dataframe <- as.data.frame(proportions_df[[gene]])
#   colnames(gene_dataframe) <- "Proportions"
#   gene_dataframe$Tissues <- rownames(gene_dataframe)
#   formatted_proportions[[gene]] <- gene_dataframe
# }