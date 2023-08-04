# script for plotting median TPM expression of SETBP1 and its targets
# across tissues as a violin plot

# load libraries
suppressPackageStartupMessages({
  set.seed(2178)
  library(here)
  library(tidyverse)
  source(here("src/expression_functions.R"))
})

ptm <- proc.time() # start timer

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

# list all subsetted files generated in `02_Setbp1_expression.R`
setbp1_files <- list.files(
  path = here(
    "results/SETBP1_Expression/"
  ),
  pattern = "_targets_tpm.csv", all.files = TRUE, full.names = TRUE
)


### SETBP1 EXPRESSION
# pull SETBP1 expression across samples
setbp1_exp <- map_dfr(setbp1_files, pull_gene_gtex)

# data wrangling for plotting
setbp1_med_tpm_vp <- setbp1_exp %>%
  pivot_longer(
    cols = everything(),
    names_to = "Tissue", values_to = "TPM"
  ) %>%
  filter(!is.na(TPM)) %>%
  mutate(
    TPM = log(1 + TPM, base = 2),
    .keep = "unused"
  ) %>% # scaling
  mutate(Tissue = str_replace_all(Tissue, " ", "")) %>% # remove whitespace
  mutate(Affected = ifelse(Tissue %in% affected, "TRUE", "FALSE")) %>%
  mutate(Tissue = str_to_title(str_replace_all(Tissue, "_", " ")))  %>%
  group_by(Tissue) %>%
  mutate(nSamples = sum(!is.na(TPM))) %>%
  mutate(Max = max(TPM, na.rm = TRUE)) #%>%
  #mutate(Median = median(TPM, na.rm = TRUE))
  
  # VIOLIN PLOT - SETBP1 TISSUES
  p <- setbp1_med_tpm_vp %>% ggplot(aes(
    x = reorder(Tissue, TPM, FUN = median),
    y = TPM, color = Affected
  )) +
  geom_violin(trim = FALSE) +
  scale_color_manual(values = c(
    "TRUE" = alpha("#2C1D6C", 0.75),
    "FALSE" = alpha("#D3B1C5", 0.5)
  )) +
  #geom_text(aes(y = Max, label = nSamples, vjust = 0.5, hjust = -1, fontface = "bold")) +
  geom_text(aes(y = Max, label = nSamples, vjust = 0.5, fontface = "bold"), nudge_y = 1.2) +
  scale_y_continuous(limits = c(-0.7, 7)) +
  ylab("TPM Scaled (1+log2)") +
  labs(x = "GTEx Tissue") +
  theme(
    #plot.margin = margin(0.1, 2.5, 0.1, 0.1, "cm"), #top, right, bottom, left
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    text = element_text(family = "Helvetica"),
    axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_text(color = "black"),
    axis.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    title = element_text(face = "bold"),
    #plot.title = element_text(size = 14, hjust = 0.5)
  )
# add boxplot
p <- p + geom_boxplot(width = 0.1) + coord_flip()
ggsave(
  here(
    "results/SETBP1_Expression/plots/setbp1_tpm_scaled_violin.png"
  ),
  height = 8, width = 6
)


# RIDGELINE PLOT - SETBP1 TISSUES
sb1_tissue_plot <- ggplot(setbp1_med_tpm_vp, aes(x = TPM, y = reorder(Tissue, TPM, FUN = median), fill = Affected)) +
  ggridges:::stat_density_ridges(quantile_lines = TRUE, quantiles = 2,
                                 alpha = 0.8) +
  scale_fill_manual(values = c(
    "TRUE" = alpha("#2C1D6C", 0.75),
    "FALSE" = alpha("#D3B1C5", 0.5)
  )) +
  geom_text(aes(x = Max, label = nSamples, fontface = "bold"), nudge_x = 1, nudge_y = 0.22) +
  #scale_y_continuous(limits = c(-0.7, 7)) +
  #geom_text(aes(y = nSamples, label = Max, vjust = 0.5, fontface = "bold"), nudge_y = 1.2) +
  scale_x_continuous(limits = c(-0.5, 7)) +
  labs(x = "TPM Scaled (1+log2)", y = "GTEx Tissue") +
  theme(
    #plot.margin = margin(0.1, 2.5, 0.1, 0.1, "cm"), #top, right, bottom, left
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    text = element_text(family = "Helvetica"),
    axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_text(color = "black"),
    axis.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    title = element_text(face = "bold"),
    #plot.title = element_text(size = 14, hjust = 0.5)
  )
ggsave(
  here(
    "results/SETBP1_Expression/plots/setbp1_tissues_ridgeplot.png"
  ),
  height = 8, width = 6
)


### SETBP1 EXPRESSION BY BRAIN REGION
# gtex metadata for brain regions
gtex_meta <- read.delim(here("data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"))
brain_meta <- filter(gtex_meta, SMTS == "Brain")
brain_tpm <- read.csv(here("results/SETBP1_Expression/BRAIN_SETBP1_targets_tpm.csv"), row.names = 1)

# find brain regions
#brain_regions <- t(brain_tpm) %>% rownames_to_column(., var = "")
brain_regions <- brain_tpm %>% 
  pivot_longer(., cols = starts_with("GTEX"), names_to = "SAMPID", values_to = "TPM") %>%
  mutate(SAMPID = str_replace(SAMPID, "\\.1$", "")) %>% # remove trailing .1's
  mutate(SAMPID = str_replace_all(SAMPID, "\\.", "-")) %>%  # sub . for _
  inner_join(., brain_meta, by = "SAMPID") %>% #combine with gtex metadata
  mutate(Region = str_remove(SMTSD, "Brain - "), .before = TPM)


# RIDGELINE PLOT - SETBP1 BRAIN REGIONS
br_plot <- brain_regions %>% mutate(TPM = log(1 + TPM, base = 2)) %>%
  ggplot(., aes(x = TPM, y = reorder(Region, TPM, FUN = median), fill = Tissue)) +
  #x = reorder(Tissue, Median_TPM, FUN = median),
  ggridges:::stat_density_ridges(quantile_lines = TRUE, quantiles = 2) +
   scale_fill_manual(values = c(alpha("#2C1D6C", 0.65))) +
  labs(x = "TPM Scaled (1+log2)", y = "Brain Region") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    text = element_text(family = "Helvetica"),
    axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_text(color = "black"),
    axis.text = element_text(face = "bold"),
    legend.position = "none",
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    title = element_text(face = "bold"),
  )
ggsave(
  here(
    "results/SETBP1_Expression/plots/setbp1_brainregions_ridgeplot.png"
  ),
  height = 8, width = 6
)

### SETBP1 TARGET EXPRESSION ACROSS TISSUES
# calculate median tpm of each gene for each tissue and compile list
tissues_median_tpm <- map_dfc(setbp1_files, get_gtex_median_tpm)

# remove any columns with NAs
tissues_median_tpm <- tissues_median_tpm %>% select_if(~ !any(is.na(.)))
# save median TPM counts
write.csv(tissues_median_tpm, here(
  "results/SETBP1_Expression/aggregated_tissues_tpm_median.csv"
))

# remove setbp1 (ENSG00000152217 = setbp1)
tissues_median_tpm <- tissues_median_tpm[!(
  row.names(tissues_median_tpm) %in% "ENSG00000152217"), ]


# data wrangling for plotting
tissues_med_tpm_vp <- tissues_median_tpm %>%
  pivot_longer(
    cols = everything(),
    names_to = "Tissue", values_to = "Median_TPM"
  ) %>%
  mutate(
    Median_TPM = log(1 + Median_TPM, base = 2),
    .keep = "unused"
  ) %>% # scaling
  mutate(
    Tissue = str_extract(
      Tissue,
      ".*(?=Median_TPM)"
    ),
    .keep = "unused"
  ) %>% # remove excess from filenames
  mutate(Tissue = str_replace_all(Tissue, " ", "")) %>% # remove whitespace
  mutate(Affected = ifelse(Tissue %in% affected, "TRUE", "FALSE")) %>%
  mutate(Tissue = str_to_title(str_replace_all(Tissue, "_", " ")))  %>%
  group_by(Tissue) %>%
  mutate(Max = max(Median_TPM, na.rm = TRUE)) 
  
  
  # VIOLIN PLOT - SETBP1 TARGETS ACROSS TISSUES
  p <- tissues_med_tpm_vp %>% ggplot(aes(
    x = reorder(Tissue, Median_TPM, FUN = median),
    y = Median_TPM, color = Affected
  )) +
  geom_violin(trim = FALSE) +
  scale_color_manual(values = c(
    "TRUE" = alpha("#2C1D6C", 0.75),
    "FALSE" = alpha("#D3B1C5", 0.5)
  )) +
  ylab("Median TPM Scaled (1+log2)") +
  #ggtitle("Expression of SETBP1 Targets Across GTEx Tissues") +
  labs(x = "GTEx Tissue") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    text = element_text(family = "Helvetica"),
    axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_text(color = "black"),
    axis.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    title = element_text(face = "bold"),
    #plot.title = element_text(size = 14, hjust = 0.5)
  )
# add boxplot
p <- p + geom_boxplot(width = 0.1) + coord_flip()
ggsave(
  here(
    "results/SETBP1_Expression/plots/setbp1targets_median_tpm_scaled_violin.png"
  ),
  height = 8, width = 6
)

# end timer
fptm <- proc.time() - ptm
fptm[3] / 60 # script runtime in minutes: 0.3915333

# save session info
saveRDS(
  sessionInfo(),
  here("results/SETBP1_Expression/versions_packages/versions_03.rds")
)