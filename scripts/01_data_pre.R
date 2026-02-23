#libraries
library(dplyr)
library(tibble)
library(tidyverse)


#read in expression data
load(file = "C:/Users/matts/Desktop/GIT_REPOS/spp1_analysis/data/RAW/SPP1_TPMstar_filt_HGNC.Rdata")

#read in GeTMM data
load(file = "C:/Users/matts/Desktop/GIT_REPOS/spp1_analysis/data/RAW/geTMM_star_HGNC.Rdata")

#read metadata
spp1_seq_meta_raw = read.table(file = "../DATA/SPP1/spp1_seq_metadata.txt", header = TRUE, sep = '\t')

spp1_seq_meta = spp1_seq_meta_raw %>% 
  rename(sample_id = Your.Sample.Name, comment = Comment, cell_line = Cell.line, spp1_profile = SPP1.profile, cisplatine = Cisplatine) %>% 
  select(id, sample_id, cell_line, spp1_profile, cisplatine, comment)

# Add replicate information and create sorting variables
metadata <- spp1_seq_meta %>%
  mutate(
    replicate = str_extract(sample_id, "\\d+$"),
    comment_order = factor(comment, levels = c("Untreated cells", "treated cells")),
    cisplatine_order = factor(cisplatine, levels = c("No", "Yes")),
    spp1_profile_order = factor(spp1_profile, 
                                levels = c("No", "Express by the cell line", 
                                           "Stable overexpression", "Recombinant protein", 
                                           "Inhibtion with siRNA"))
  ) %>%
  arrange(replicate, comment_order, cisplatine_order, spp1_profile_order)

#log2 transform expression data
spp1_salmon_tpm_log2 = log2(SPP1_TPMstar_filt_HGNC + 1)
geTMM_star_HGNC_log2 = log2(geTMM_star_HGNC + 1)

#export data
expr_data = as.data.frame(spp1_salmon_tpm_log2)
expr_data_getmm = as.data.frame(geTMM_star_HGNC_log2)
save(expr_data, file = "C:/Users/matts/Desktop/GIT_REPOS/spp1_analysis/data/expr_data.Rdata")
save(expr_data_getmm, file = "C:/Users/matts/Desktop/GIT_REPOS/spp1_analysis/data/expr_data_getmm.Rdata")
save(metadata, file = "C:/Users/matts/Desktop/GIT_REPOS/spp1_analysis/data/metadata.Rdata")

