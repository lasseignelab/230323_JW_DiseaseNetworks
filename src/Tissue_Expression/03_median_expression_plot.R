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

# data wrangling for violin plot
setbp1_med_tpm_vp <- setbp1_exp %>%
  pivot_longer(
    cols = everything(),
    names_to = "Tissue", values_to = "TPM"
  ) %>%
  mutate(
    TPM = log(1 + TPM, base = 2),
    .keep = "unused"
  ) %>% # scaling
  mutate(Tissue = str_replace_all(Tissue, " ", "")) %>% # remove whitespace
  mutate(Affected = ifelse(Tissue %in% affected, "TRUE", "FALSE")) %>%
  mutate(Tissue = str_to_title(str_replace_all(Tissue, "_", " "))) # %>%

# plot
p <- setbp1_med_tpm_vp %>% ggplot(aes(
  x = reorder(Tissue, TPM, FUN = median),
  y = TPM, color = Affected
)) +
  geom_violin(trim = FALSE) +
  scale_color_manual(values = c(
    "TRUE" = alpha("#2C1D6C", 0.75),
    "FALSE" = alpha("#D3B1C5", 0.5)
  )) +
  ylab("TPM Scaled (1+log2)") +
  ggtitle("Expression of SETBP1 Across GTEx Tissues") +
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
    plot.title = element_text(size = 14, hjust = 0.5)
  )
# add boxplot
p <- p + geom_boxplot(width = 0.1) + coord_flip()
ggsave(
  here(
    "results/SETBP1_Expression/plots/setbp1_tpm_scaled_violin.png"
  ),
  height = 8, width = 6
)

### SETBP1 TARGET EXPRESSION
# calculate median tpm of each gene for each tissue and compile list
tissues_median_tpm <- map_dfc(setbp1_files, get_gtex_median_tpm)

# remove any columns with NAs
tissues_median_tpm <- tissues_median_tpm %>% select_if(~ !any(is.na(.)))
# save median TPM counts
write.csv(tissues_median_tpm, here(
  "results/SETBP1_Expression/aggregated_tissues_tpm_median.csv"
))


## VIOLIN PLOT
# remove setbp1 (ENSG00000152217 = setbp1)
tissues_median_tpm <- tissues_median_tpm[!(
  row.names(tissues_median_tpm) %in% "ENSG00000152217"), ]

# data wrangling for violin plot
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
  mutate(Tissue = str_to_title(str_replace_all(Tissue, "_", " "))) # %>%
# plot
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
  ggtitle("Expression of SETBP1 Targets Across GTEx Tissues") +
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
    plot.title = element_text(size = 14, hjust = 0.5)
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
