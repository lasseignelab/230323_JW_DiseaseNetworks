# script for plotting median TPM expression of SETBP1 and its targets across tissues

# load libraries
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(ComplexHeatmap)
  source(here("src/expression_functions.R"))
})

# list all subsetted files
setbp1_files <- list.files(path = here("results/SETBP1_Expression/"), pattern = "_targets_tpm.csv", all.files = TRUE, full.names = TRUE)

# calculate median tpm of each gene for each tissue and compile list
tissues_median_tpm <- map_dfc(setbp1_files, get_gtex_median_tpm)
# remove any columns with NAs
tissues_median_tpm <- tissues_median_tpm %>% select_if(~ !any(is.na(.)))

#affected tissues 
## affected tissues list
affected <- c(
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

## VIOLIN PLOT
p <- tissues_median_tpm %>% 
  pivot_longer(cols = everything(), names_to = "Tissue", values_to = "Median_TPM") %>% 
  mutate(Median_TPM = log(1 + Median_TPM, base = 2), .keep = "unused") %>% #scaling
  mutate(Tissue = str_extract(Tissue, ".*(?=Median_TPM)"), .keep = "unused") %>% #remove excess from filenames
  mutate(Tissue = str_replace_all(Tissue, " ", "")) %>% #remove whitespace from names
  mutate(Affected = ifelse(Tissue %in% affected, "TRUE", "FALSE")) %>%
  mutate(Tissue = str_to_title(str_replace_all(Tissue, "_", " "))) %>%
  ggplot(aes(x = reorder(Tissue, Median_TPM, FUN = median), y = Median_TPM, color = Affected)) + 
  geom_violin(trim = FALSE) + 
  scale_color_manual(values = c("TRUE" = alpha("#2C1D6C", 0.75), "FALSE" = alpha("#D3B1C5", 0.5))) +
  ylab("Median TPM Scaled (1+log2)") +
  ggtitle("Expression of SETBP1 and Targets Across GTEx Tissues") +
  #theme_minimal() +
  labs(x = "GTEx Tissue") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        text = element_text(family = "Helvetica"),
        #axis.text.x = element_text(color = "black", size = 10, angle = 80, vjust = 0.4, hjust = 0.2),
        axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black"),
        axis.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        title = element_text(face = "bold"),
        plot.title = element_text(size = 14, hjust = 0.5))
#add boxplot 
p <- p + geom_boxplot(width = 0.1) + coord_flip()
ggsave(here("results/SETBP1_Expression/plots/median_tpm_scaled_violin_flip.png"), height = 8, width = 6)


# gene labels and data wrangling
setbp1_targets <- read.csv(file = here("data/SETBP1_Targets/setbp1_targets_geneconversions.csv"), row.names = 1)

tissues_med_log <- tissues_median_tpm %>% 
  rownames_to_column(., var = "target") %>% 
  purrr::modify_if(., is.numeric, onelog2) %>% # 1 + log2 transformed for scaling
#  mutate(Tissue = str_extract(Tissue, ".*(?=Median_TPM)"), .keep = "unused") %>%
  left_join(., setbp1_targets, by = "target") %>% #combine TPM with gene annotations
  column_to_rownames(., var = "target") %>%
  dplyr::select(., name | ends_with("Median_TPM")) %>%
  rename_all( ~ gsub("Median_TPM", "", .)) %>% 
  rename_all( ~ str_to_title(gsub("_", " ", .)))
#  mutate(Affected = ifelse(Tissue %in% affected, "TRUE", "FALSE")) %>%
#  mutate(Tissue = str_to_title(str_replace_all(Tissue, "_", " "))) %>%

## HEATMAP
# everything in one heatmap
png(here("results/SETBP1_Expression/plots/median_tpm_scaled_heatmap.png"), width = 25, height = 30, units = "cm", res = 300)
set.seed(1)
#ComplexHeatmap::Heatmap(as.matrix(tissues_med_log[2:length(tissues_med_log)]), row_labels = tissues_med_log$name) 

ComplexHeatmap::Heatmap(as.matrix(tissues_med_log[2:length(tissues_med_log)])) + 
                        #rowAnnotation(link = anno_mark(tissues_med_log$name[tissues_med_log$name == "SETBP1"],#Error: `at` should be numeric row index corresponding to the matrix. 
rowAnnotation(link = anno_mark(at = which(tissues_med_log$name == "SETBP1"), 
                            labels = tissues_med_log$name[tissues_med_log$name == "SETBP1"] ))#, 
#                            labels_gp = gpar(fontsize = 10), padding = unit(1, "mm"))))
                        #row_labels = tissues_med_log$name[tissues_med_log$name == "SETBP1"]) 
#rowAnnotation(link = anno_mark(at = which(ccl & base_mean > quantile(base_mean, 0.25)), 
#labels = rownames(mat)[ccl & base_mean > quantile(base_mean, 0.25)], 
#labels_gp = gpar(fontsize = 10), padding = unit(1, "mm"))) +
dev.off()

# setbp1 as separate heatmap
#setbp1_tissue <- tissues_med_log %>% dplyr::filter(., Name == "SETBP1") #%>% t() #weird class change issues with numbers..
#t <- tissues_med_log[4:5,]# %>% dplyr::filter(., name == "SETBP1") #%>% t()
targets_tissue <- tissues_med_log %>% dplyr::filter(., Name != "SETBP1")# %>% t()


#set color breaks
break_quart <- quantile(as.matrix(tissues_med_log[,2:length(tissues_med_log)]))

col_fun = circlize::colorRamp2(break_quart, c("white", "#D3B1C5", "#A468A9", "#533180", "#2C1D6C"))
#tf_setbp1_cortex_corr <- corrplot(tf_jaccard_setbp1_matrix, method = "circle", is.corr = FALSE, type = "upper", diag = FALSE, col = my_palette(100))
tempdf <- t(as.matrix(tissues_med_log[1:4,2:length(tissues_med_log)]))
colnames(tempdf) <- tissues_med_log$Name[1:4]

#plot and save
png(here("results/SETBP1_Expression/plots/median_tpm_scaled_heatmap_clusters.png"), width = 25, height = 30, units = "cm", res = 300)
set.seed(1)
#top_ha <- ComplexHeatmap::HeatmapAnnotation(SETBP1 = t(as.matrix(tissues_med_log[4,2:length(tissues_med_log)])),
top_ha <- ComplexHeatmap::HeatmapAnnotation(SETBP1 = tempdf[,4],
                                            name = "SETBP1",
                                            col = list(SETBP1 = col_fun),
                                            #annotation_name_gp = gpar(names = "SETBP1"),
                                            #show_annotation_name = FALSE,
                                            show_legend = FALSE,
                                            annotation_label = "SETBP1")
hm_sb1_tar <- ComplexHeatmap::Heatmap(as.matrix(targets_tissue[,2:length(targets_tissue)]), 
                                      km = 3,
                                      col = col_fun,
                                      row_labels = targets_tissue[,1],
                                      row_names_gp = gpar(fontsize = 6),
                                      heatmap_legend_param = list(title = "Scaled  TPM"),
                                      top_annotation = top_ha)

#hm_sb1 + hm_sb1_tar
hm_sb1_tar
dev.off()

#adapted code for pulling clusters from here https://github.com/jokergoo/ComplexHeatmap/issues/136#issuecomment-671931967
dend <- row_dend(draw(hm_sb1_tar))
clust_list <- row_order(draw(hm_sb1_tar))

clu_df <- lapply(names(clust_list), function(i){
   out <- data.frame(GeneID = targets_tissue[clust_list[[i]],1],#rownames(mat[clust_list[[i]],]),
                     Cluster = paste0("cluster", i),
                     stringsAsFactors = FALSE)
   return(out)
 }) %>%  #pipe (forward) the output 'out' to the function rbind to create 'clu_df'
  do.call(rbind, .)

# full gene annotations for each cluster
clu_df_anno <- clu_df %>% inner_join(., setbp1_targets, by = c("GeneID" = "name"), multiple = "all")

#FEA of each cluster
cl1_fea <- clu_df %>% dplyr::filter(., Cluster == "cluster1") #%>% 
cl1_fea <- gprofiler2::gost(cl1_fea$GeneID, evcodes = TRUE)

cl2_fea <- clu_df %>% dplyr::filter(., Cluster == "cluster2") #%>% 
cl2_fea <- gprofiler2::gost(cl2_fea$GeneID, evcodes = TRUE)


cl3_fea <- clu_df %>% dplyr::filter(., Cluster == "cluster3") #%>% 
cl3_fea <- gprofiler2::gost(cl3_fea$GeneID, evcodes = TRUE)

#cl4_fea <- clu_df %>% dplyr::filter(., Cluster == "cluster4") #%>% 
#cl4_fea <- gprofiler2::gost(cl4_fea$GeneID, evcodes = TRUE)
