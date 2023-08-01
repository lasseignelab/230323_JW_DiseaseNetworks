# function to pull GTEx TPM values for one gene across tissues
pull_gene_gtex <- function(file, gene){
  tpm <- read.csv(file, row.name = 1)
  gene_tpm <- tpm %>% 
    dplyr::filter(., name == "SETBP1") %>% # filter for just gene of interest
    dplyr::select(., "Tissue", starts_with("GTEX"), starts_with("K.562")) %>% 
    pivot_longer(., cols = c(starts_with("GTEX"), starts_with("K.562")), values_to = unique(tpm$Tissue), values_drop_na = TRUE) #pivot for easier combining across tissues #pivot for easier combining across tissues
  
return(gene_tpm[,3]) #select just new column
}


# function to calculate median TPM values across GTEx samples for an input file
# file = filepath of GTEx TPM csv file
# cols_tpm = beginning identifier of desired columns' colnames to calculate median across
# identifier = desired new column location -- indicate which column to place the new column after
calc_median_tpm <- function(file, cols_tpm = "GTEX", identifier = "Tissue"){
  tpm <- read.csv(file, row.name = 1)
  #calculate median from 5:length(tpm) , mutate column ... across columns
  #median_tpm <- dplyr::mutate(rowwise(tpm), median = median(c(5:length(tpm))), .after = "Tissue")
  
  #use cols_tpm to select columns to calculate tpm rowwise
  median_tpm <- tpm %>% dplyr::select(., starts_with(cols_tpm)) %>% apply(., 1, median)
  #add calculated median TPM values to gtex df
  median_tpm <- tpm %>% mutate(., Median = median_tpm, .after = identifier)
  
return(median_tpm)  
}

# data wrangling specific to format from these analyses for easier plotting
# median_tpm_df = output from calc_median_tpm
# cols_tpm = beginning identifier of TPM counts to remove for minimal 
# identifier = column with same string throughout -- used as colname for Median values to make aggregation across results easier
pull_median_tpm <- function(median_tpm_df, identifier = "Tissue"){
  #print(paste(unique(median_tpm_df[identifier])))
  #id <- unique(median_tpm_df[identifier])
  #colnames(median_tpm_df)["Median"] <- paste(unique(median_tpm_df[identifier]))
  #rename the Median column to indicate the tissue
  colnames(median_tpm_df)[5] <- paste(unique(median_tpm_df[identifier]), "Median_TPM", collapse = "_")
  #select minimally needed columns for plotting
  median_tpm_simple <- median_tpm_df %>% 
    dplyr::select(., c(target) | ends_with("Median_TPM")) %>% 
    column_to_rownames(., var = "target")
return(median_tpm_simple)
}

# wrapper for calc_median_tpm and pull_median_tpm
get_gtex_median_tpm <- function(file, cols_tpm = "GTEX", identifier = "Tissue"){
  median_tpm_df <- calc_median_tpm(file, cols_tpm = cols_tpm, identifier = identifier)
  median_tpm_simple <- pull_median_tpm(median_tpm_df, identifier = identifier)
  
return(median_tpm_simple)
}

#function to calculate 1 + log(2)
#1 + log base 2 fxn
onelog2 <- function(x) {
  log(1 + x, base = 2)
}

