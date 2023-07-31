# pull together human SETBP1 targets across data sets that JW found
############# START JW's notes#############
#* Obtained gene list of targets from msigdb SETBP1 target list
#* (https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?
#* geneSetName=SETBP1_TARGET_GENES&fileType=TSV). The target list stored
#* in cheaha at '221018_Setbp1_msigdb_targets.csv' at /data/user/
#* jbarham3/2302227_JW_Setbp1Manuscript/results/seurat/, then converted
#* from mouse to human using bioDBnet Ortho GUI. Results are in file
#* '221019_bioDBnet_dbOrtho_msigdb.txt'
#* Additional genes were compiled from the below resources into '221019_
#* bioDBnet_dbOrtho_lit_signor.txt':
#  1. the literature (papers here): HCF1, KMT2A, PHF8, PHF6, SET, HOXA9,
#  HOXA10, PP2A, HMGA1, PHF6, BMP5, PDE4D, ERPP4, RUNX1, TCF4, BCL11A,
#  DNTT, MYB, MYC, CEBPB, PARP1, TAF1A, ANP32A, AKT, HMG2, TREX1, NME1, APEX1
# 2. Signor:
#  * Human: SET, HOXA9, HOXA10


# The final target list is composed of 209 Human genes mapped to 682
# mouse orthologs/aliases.
#############END JW's notes#############

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

ptm <- proc.time() # start timer

#setbp1 targets from msigdb
msig_setbp1 <- msigdbr::msigdbr(species = "Homo sapiens",
                                category = "C3",
                                subcategory = "GTRD")
msig_setbp1 <- msig_setbp1 %>% filter(., gs_name == "SETBP1_TARGET_GENES")

#targets compiled from literature by JW + SETBP1 itself
lit_setbp1 <- c("SETBP1",
                "HCF1",
                "KMT2A",
                "PHF8",
                "PHF6",
                "SET",
                "HOXA9",
                "HOXA10",
                "PP2A",
                "HMGA1",
                "PHF6",
                "BMP5",
                "PDE4D",
                "ERPP4",
                "RUNX1",
                "TCF4",
                "BCL11A",
                "DNTT",
                "MYB",
                "MYC",
                "CEBPB",
                "PARP1",
                "TAF1A",
                "ANP32A",
                "AKT",
                "HMG2",
                "TREX1",
                "NME1",
                "APEX1")


# signor
signor_setbp1 <- c("SET", "HOXA9", "HOXA10")
# any in signor that are no in the literature list of targets?
setdiff(signor_setbp1, lit_setbp1) # none, just using lit list moving forward

# any in the literature list of targets that are not in the msigdb one?
setdiff(lit_setbp1, msig_setbp1$gene_symbol)
length(setdiff(lit_setbp1, msig_setbp1$gene_symbol)) # 28


# convert to ensembl ID's
setbp1_targets_geneconversions <- gprofiler2::gconvert(
  unique(c(signor_setbp1, lit_setbp1, msig_setbp1$gene_symbol)),
  target = "ENSG"
)

# save setbp1 aggregated targets and gene conversions
write.csv(setbp1_targets_geneconversions,
  file = here("data/SETBP1_Targets/setbp1_targets_geneconversions.csv")
)

# end timer
fptm <- proc.time() - ptm
fptm[3] / 60 # script runtime in minutes

# save session info
saveRDS(
  sessionInfo(),
  here("results/SETBP1_Expression/versions_packages/versions_01.rds")
)

