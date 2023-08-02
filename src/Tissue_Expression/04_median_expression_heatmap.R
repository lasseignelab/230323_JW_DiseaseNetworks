# script for plotting median TPM expression of SETBP1 and its targets across
# tissues as a heatmap

# load libraries
suppressPackageStartupMessages({
  set.seed(2178)
  library(here)
  library(tidyverse)
  library(ComplexHeatmap)
  source(here("src/expression_functions.R"))
})

ptm <- proc.time() # start timer

# list all subsetted files
setbp1_files <- list.files(
  path = here(
    "results/SETBP1_Expression/"
  ), pattern = "_targets_tpm.csv", all.files = TRUE,
  full.names = TRUE
)

# read in median TPM counts
tissues_median_tpm <- read.csv(here(
  "results/SETBP1_Expression/aggregated_tissues_tpm_median.csv"
), row.names = 1)

# affected tissues
affected <- c(
  "BLOOD",
  "BRAIN",
  "HEART",
  "KIDNEY",
  "BLADDER",
  "LUNG",
  "MUSCLE",
  "SMALL_INTESTINE",
  "STOMACH",
  "ESOPHAGUS"
)

# load in gene labels
setbp1_targets <- read.csv(file = here(
  "data/SETBP1_Targets/setbp1_targets_geneconversions.csv"
), row.names = 1)

# data wrangling for complexheatmap input
tissues_med_log <- tissues_median_tpm %>%
  rownames_to_column(., var = "target") %>%
  purrr::modify_if(., is.numeric, onelog2) %>% # 1 + log2 transformed scaling
  left_join(., setbp1_targets, by = "target") %>% # combine TPM w/ gene anno
  column_to_rownames(., var = "target") %>%
  dplyr::select(., name | ends_with("Median_TPM")) %>%
  rename_all(~ gsub(".Median_TPM", "", .)) %>%
  rename_all(~ str_to_title(gsub("_", " ", .)))


## Heatmap version 1:
# everything in one heatmap (not split)
png(
  here(
    "results/SETBP1_Expression/plots/median_tpm_scaled_heatmap.png"
  ),
  width = 25,
  height = 30, units = "cm", res = 300
)
set.seed(1)

ComplexHeatmap::Heatmap(
  as.matrix(tissues_med_log[2:length(tissues_med_log)])
) +
  rowAnnotation(link = anno_mark(
    at = which(tissues_med_log$name == "SETBP1"),
    labels = tissues_med_log$name[tissues_med_log$name == "SETBP1"]
  ))

dev.off()

## Heatmap version 2:
# setbp1 separated from targets in heatmap
targets_tissue <- tissues_med_log %>% dplyr::filter(., Name != "SETBP1")


# set color breaks
break_quart <- quantile(as.matrix(tissues_med_log[, 2:length(tissues_med_log)]))
col_fun <- circlize::colorRamp2(
  break_quart,
  c(
    "white",
    "#D3B1C5",
    "#A468A9",
    "#533180",
    "#2C1D6C"
  )
)

# Decide number of k means clusters, adapted code from https://www.r-bloggers.com/2020/05/how-to-determine-the-number-of-clusters-for-k-means-in-r/#:~:text=We%20can%20determine%20the%20number,K%20which%20maximizes%20that%20score.
# Use map_dbl to run many models with varying value of k (centers)
tot_withinss <- map_dbl(1:10,  function(k){
  model <- kmeans(x = tissues_med_log[, 2:length(tissues_med_log)], centers = k)
  model$tot.withinss
})

# Generate a data frame containing both k and tot_withinss
elbow_df <- data.frame(
  k = 1:10,
  tot_withinss = tot_withinss
)
# Plot the elbow plot
ggplot(elbow_df, aes(x = k, y = tot_withinss)) +
  geom_line() + geom_point()+
  scale_x_continuous(breaks = 1:10)


# pulling first 4 columns to wrangle setbp1 (4th row) to use as annotation
tempdf <- t(as.matrix(tissues_med_log[1:4, 2:length(tissues_med_log)]))
colnames(tempdf) <- tissues_med_log$Name[1:4]

# plot and save
png(
  here(
    "results/SETBP1_Expression/plots/median_tpm_scaled_heatmap_3clusters.png"
  ),
  width = 25, height = 30, units = "cm", res = 300 # height = 50 ensures legible y-axis
)
set.seed(1)
ht_opt(heatmap_column_names_gp = gpar(fontface = "bold"))
top_ha <- ComplexHeatmap::HeatmapAnnotation(
  SETBP1 = tempdf[, 4],
  name = "SETBP1",
  col = list(SETBP1 = col_fun),
  show_legend = FALSE,
  #gp = gpar(fontface = "bold"),
  annotation_label = "SETBP1"
)
hm_sb1_tar <- ComplexHeatmap::Heatmap(
  as.matrix(targets_tissue[, 2:length(targets_tissue)]),
  km = 3, 
  col = col_fun,
  #row_labels = targets_tissue[, 1],
  #row_names_gp = gpar(fontsize = 6),
  show_row_names = FALSE,
  heatmap_legend_param = list(title = "Scaled TPM"),
  top_annotation = top_ha
)

hm_sb1_tar
dev.off()

# adapted code for pulling clusters from here
# https://github.com/jokergoo/ComplexHeatmap/issues/136#issuecomment-671931967
dend <- row_dend(draw(hm_sb1_tar))
clust_list <- row_order(draw(hm_sb1_tar))

clu_df <- lapply(names(clust_list), function(i) {
  out <- data.frame(
    GeneID = targets_tissue[clust_list[[i]], 1],
    Cluster = paste0("cluster", i),
    stringsAsFactors = FALSE
  )
  return(out)
}) %>% # pipe (forward) the output 'out' to the function rbind
  do.call(rbind, .)

# FEA of each cluster
cluster_fea <- enrich_clusters(clu_df, custom_bg = clu_df$GeneID)

p <- cluster_fea %>% 
  dplyr::filter(., p_value < 0.05) %>%
  ggplot(., aes(x = Cluster, y = reorder(term_name, -p_value), size = recall, fill =
                                     p_value)) +
  geom_point(alpha = 0.7, shape = 21) +
  scale_size(range = c(2, 10), name = "Recall") + 
  scale_fill_distiller(palette = "Purples") + 
  labs(x = "Recall", y = "Terms") +
  theme(
    panel.grid.major = element_line("gray95"),
    #panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    text = element_text(family = "Helvetica"),
    #axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_text(color = "black"),
    axis.text.x = element_text(color = "black"),
    axis.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    title = element_text(face = "bold"),
    plot.title = element_text(size = 14, hjust = 0.5)
  )
p
ggsave(here("results/SETBP1_Expression/plots/fea_3clusters_bubbleplot.png"), width = 8, height = 9, bg = "white")

# end timer
fptm <- proc.time() - ptm
fptm[3] / 60 # script runtime in minutes: 0.19175

# save session info
saveRDS(
  sessionInfo(),
  here("results/SETBP1_Expression/versions_packages/versions_04.rds")
)
