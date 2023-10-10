# script for plotting median TPM expression of SETBP1 and its targets across
# tissues as a heatmap

# load libraries
suppressPackageStartupMessages({
  set.seed(2178)
  library(here)
  library(tidyverse)
  library(ComplexHeatmap)
  library(stats)
  library(NbClust)
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
  "ESOPHAGUS",
  "BONE_MARROW"
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

# set color breaks for heatmap
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

# setbp1 separated from targets in heatmap
targets_tissue <- tissues_med_log %>% dplyr::filter(., Name != "SETBP1")

# pulling first 4 columns to wrangle setbp1 (4th row) to use as heatmap anno
tempdf <- t(as.matrix(tissues_med_log[1:4, 2:length(tissues_med_log)]))
colnames(tempdf) <- tissues_med_log$Name[1:4]


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

top_ha <- ComplexHeatmap::HeatmapAnnotation(
  SETBP1 = tempdf[, 4],
  name = "SETBP1",
  col = list(SETBP1 = col_fun),
  show_legend = FALSE,
  annotation_label = "SETBP1"
)
hm_sb1_tar <- ComplexHeatmap::Heatmap(
  as.matrix(targets_tissue[, 2:length(targets_tissue)]),
  col = col_fun,
  show_row_names = FALSE,
  heatmap_legend_param = list(title = "Scaled TPM"),
  top_annotation = top_ha
)
ht_opt(heatmap_column_names_gp = gpar(fontface = "bold"))
hm_sb1_tar
dev.off()



## Heatmap version 2:

### Decide number of k means clusters
set.seed(2178)
opt_clust <- NbClust(tissues_med_log[, -1], min.nc = 2, max.nc = 15, method = "complete", index = "tracew")
opt_clust$Best.nc # optimal cluster and index value


tracew_p <- data.frame(opt_clust$All.index) %>%
  #rownames_to_column(var = "cluster")# %>%
  ggplot(aes(x = 2:15, y = opt_clust.All.index)) +
  geom_line() + 
  geom_point() +
  #scale_color_manual(values = c("Max_Second_Diff" = "steelblue")) +
  geom_vline(xintercept = 3, linetype = "dashed", color = "steelblue", size = 0.8) + # optimal cluster from "Best.nc"
  scale_x_continuous(breaks = 1:20) +
  xlab("Cluster") +
  ylab("Trace(W)")
tracew_p
ggsave("./results/SETBP1_Expression/plots/kmeans_traceW.png", width = 10, height = 8, units = "cm")

 # Elbow Method
# adapted code from https://www.r-bloggers.com/2020/05/how-to-determine-the-number-of-clusters-for-k-means-in-r/#:~:text=We%20can%20determine%20the%20number,K%20which%20maximizes%20that%20score.
# Use map_dbl to run many models with varying value of k (centers)
tot_withinss <- map_dbl(1:15,  function(k){
  set.seed(2178)
  model <- kmeans(x = tissues_med_log[, -1], #counts w/o name symbol anno col
                  centers = k)
  model$tot.withinss
})

# Generate a data frame containing both k and tot_withinss
elbow_df <- data.frame(
  k = 1:15,
  tot_withinss = tot_withinss
)
# Plot the elbow plot
ggplot(elbow_df, aes(x = k, y = tot_withinss)) +
  geom_line() + 
  geom_point() +
  scale_x_continuous(breaks = 1:20) +
  geom_vline(xintercept = 3, linetype = "dashed", color = "steelblue", size = 0.8) + # optimal cluster from bend
  xlab("Cluster") +
  ylab("Within SS")

ggsave(here("results/SETBP1_Expression/plots/elbowplot.png"), width = 12, height = 8, units = "cm")


# plot and save
png(
  here(
    "results/SETBP1_Expression/plots/median_tpm_scaled_heatmap_3clusters.png"
  ),
  width = 25, height = 30, units = "cm", res = 300 # height = 50 ensures legible y-axis
)
set.seed(1)

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
ht_opt(heatmap_column_names_gp = gpar(fontface = "bold"))
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
ggsave(
  here("results/SETBP1_Expression/plots/fea_3clusters_bubbleplot.png"),
  width = 8, height = 9, bg = "white")

# correlations
setbp1 <- t(as.matrix(tissues_med_log[4, 2:length(tissues_med_log)]))
targets <- t(as.matrix(tissues_med_log[-4, 2:length(tissues_med_log)]))

#cor_res <- list()
cor_res <- NULL

for(i in 1:ncol(targets)){
  res <- cor.test(setbp1, targets[,i])
  temp <- data.frame(cor = res$estimate, p_value = res$p.value)
  cor_res <- rbind(cor_res, temp)
  #cor_res[[colnames(targets)[i]]] <- temp
  rownames(cor_res)[i] <- colnames(targets)[i]
}

#t <- t(as.data.frame(cor_res))
sp1_target_cor <- cor_res %>%
  rownames_to_column(., var = "target") %>%
  inner_join(., setbp1_targets, by = "target")

sb1_cor_cluster <- left_join(clu_df, sp1_target_cor,
                             by = c("GeneID" = "name")) %>%
  select(GeneID, Cluster, SETBP1_Cor = cor, Cor_P_Value = p_value, description)
write.csv(sb1_cor_cluster, here(
  "results/SETBP1_Expression/setbp1_target_cor_clust.csv"))

# end timer
fptm <- proc.time() - ptm
fptm[3] / 60 # script runtime in minutes: 0.19175

# save session info
saveRDS(
  sessionInfo(),
  here("results/SETBP1_Expression/versions_packages/versions_04.rds")
)