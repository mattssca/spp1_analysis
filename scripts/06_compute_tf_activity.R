####################################################################################################
## compute_tf_activity
## Compute TF activity scores using decoupleR + CollecTRI
## Input:  expr_data_getmm.Rdata, metadata.Rdata
## Output: data/acts_mat.Rdata
## Note:   Run once; load acts_mat.Rdata in downstream scripts
####################################################################################################

library(decoupleR)
library(tidyverse)

####################################################################################################
## Load Data
####################################################################################################

load("C:/Users/matts/Desktop/GIT_REPOS/spp1_analysis/data/expr_data_getmm.Rdata")
load("C:/Users/matts/Desktop/GIT_REPOS/spp1_analysis/data/metadata.Rdata")

metadata$sample_base <- sub("_[0-9]+$", "", metadata$sample_id)

####################################################################################################
## Get CollecTRI Regulon Network
####################################################################################################

net <- get_collectri(organism = "human", split_complexes = FALSE)

cat("CollecTRI: ", length(unique(net$source)), "TFs,",
    length(unique(net$target)), "target genes\n")

####################################################################################################
## Run TF Activity Inference (ULM)
####################################################################################################

mat <- as.matrix(expr_data_getmm)

cat("Running ULM on", nrow(mat), "genes x", ncol(mat), "samples...\n")

tf_acts <- run_ulm(
  mat     = mat,
  net     = net,
  .source = "source",
  .target = "target",
  .mor    = "mor",
  minsize = 5
)

####################################################################################################
## Reshape to TF x Sample Matrix
####################################################################################################

acts_mat <- tf_acts %>%
  filter(statistic == "ulm") %>%
  select(source, condition, score) %>%
  pivot_wider(names_from = condition, values_from = score) %>%
  column_to_rownames("source") %>%
  as.matrix()

# Align column order to metadata
acts_mat <- acts_mat[, metadata$sample_id, drop = FALSE]

cat("acts_mat dimensions:", nrow(acts_mat), "TFs x", ncol(acts_mat), "samples\n")

####################################################################################################
## Save
####################################################################################################

save(acts_mat, net, file = "C:/Users/matts/Desktop/GIT_REPOS/spp1_analysis/data/acts_mat.Rdata")

cat("Saved -> data/acts_mat.Rdata (contains: acts_mat, net)\n")

####################################################################################################
## compute_tf_genesets
## Derive TF gene sets from DEG results using CollecTRI regulon network
## Input:  DEG_results.xlsx, acts_mat.Rdata (contains net)
## Output: data/tf_genesets.Rdata
## Note:   Run after 05b_compute_tf_activity.R
####################################################################################################

library(readxl)
library(dplyr)

####################################################################################################
## Load
####################################################################################################

load("C:/Users/matts/Desktop/GIT_REPOS/spp1_analysis/data/acts_mat.Rdata")  # loads: acts_mat, net

sheet_map <- c(
  hm1 = "VMCUB1_CISPL_vs_VMCUB1_SISPP1_C",
  hm2 = "VMCUB1_CTRL_vs_VMCUB1_CISPL",
  hm3 = "VMCUB1_CTRL_vs_VMCUB1_SISPP1",
  hm4 = "SCABER_CISPL_vs_SCABER_OESPP1_C",
  hm5 = "SCABER_CTRL_vs_SCABER_CISPL",
  hm6 = "SCABER_CTRL_vs_SCABER_OESPP1",
  hm7 = "SCABER_CTRL_vs_SCABER_RSPP1",
  hm8 = "SCABER_CISPL_vs_SCABER_RSPP1_CI"
)

deg_list <- lapply(sheet_map, function(sheet) {
  read_excel(
    "C:/Users/matts/Desktop/GIT_REPOS/spp1_analysis/DEG_results.xlsx",
    sheet = sheet
  )
})

####################################################################################################
## Build Gene Sets
####################################################################################################
# All 8 contrasts
all_sig_genes <- bind_rows(deg_list) %>%
  filter(adj.P.Val < 0.05) %>%
  pull(gene) %>%
  unique()

####################################################################################################
## Map to TFs via CollecTRI
####################################################################################################

tfs_all_deg    <- net %>% filter(target %in% all_sig_genes)    %>% pull(source) %>% unique()

cat("Gene sets:\n")
cat("  all_sig_genes:   ", length(all_sig_genes),    "genes\n")
cat("TF sets:\n")
cat("  tfs_all_deg:   ", length(tfs_all_deg),    "TFs\n")

####################################################################################################
## Save
####################################################################################################

save(
all_sig_genes,
tfs_all_deg,
  file = "C:/Users/matts/Desktop/GIT_REPOS/spp1_analysis/data/tf_genesets.Rdata"
)

cat("Saved -> data/tf_genesets.Rdata\n")
cat("  Objects: cispl_spp1_genes, all_sig_genes, tfs_cispl_spp1, tfs_all_deg\n")
