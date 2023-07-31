# script for plotting median TPM expression of SETBP1 and its targets
# across tissues as a violin plot

# load libraries
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  source(here("src/expression_functions.R"))
})

ptm <- proc.time() # start timer

# list all subsetted files
setbp1_files <- list.files(
  path = here(
    "results/SETBP1_Expression/"
  ),
  pattern = "_targets_tpm.csv", all.files = TRUE, full.names = TRUE
)

# calculate median tpm of each gene for each tissue and compile list
tissues_median_tpm <- map_dfc(setbp1_files, get_gtex_median_tpm)

# remove any columns with NAs
tissues_median_tpm <- tissues_median_tpm %>% select_if(~ !any(is.na(.)))
# save median TPM counts
write.csv(tissues_median_tpm, here(
  "results/SETBP1_Expression/aggregated_tissues_tpm_median.csv"
))

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

## VIOLIN PLOT
p <- tissues_median_tpm %>%
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
  mutate(Tissue = str_to_title(str_replace_all(Tissue, "_", " "))) %>%
  ggplot(aes(
    x = reorder(Tissue, Median_TPM, FUN = median),
    y = Median_TPM, color = Affected
  )) +
  geom_violin(trim = FALSE) +
  scale_color_manual(values = c(
    "TRUE" = alpha("#2C1D6C", 0.75),
    "FALSE" = alpha("#D3B1C5", 0.5)
  )) +
  ylab("Median TPM Scaled (1+log2)") +
  ggtitle("Expression of SETBP1 and Targets Across GTEx Tissues") +
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
    "results/SETBP1_Expression/plots/median_tpm_scaled_violin_flip.png"
  ),
  height = 8, width = 6
)

# end timer
fptm <- proc.time() - ptm
fptm[3] / 60 # script runtime in minutes

# save session info
saveRDS(
  sessionInfo(),
  here("results/SETBP1_Expression/versions_packages/versions_03.rds")
)
