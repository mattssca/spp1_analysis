#libraries
library(dplyr)
library(tibble)

#read in expression data
spp1_salmon_tpm_raw = read.table("../../../DATA/SPP1/RNA-seq/salmon.merged.gene_tpm.tsv", header = TRUE)

#read metadata
spp1_seq_meta_raw = read.table(file = "../../../DATA/SPP1/spp1_seq_metadata.txt", header = TRUE, sep = '\t')

#format data
spp1_salmon_tpm = spp1_salmon_tpm_raw %>% 
  select(-gene_id) %>% 
  rename(hgnc_symbol = gene_name)

#aggregate duplicates by summing expression values
spp1_salmon_tpm_aggregated <- spp1_salmon_tpm %>%
  summarize(across(where(is.numeric), sum), .by = hgnc_symbol) %>% 
  column_to_rownames("hgnc_symbol")

spp1_seq_meta = spp1_seq_meta_raw %>% 
  rename(sample_id = Your.Sample.Name, comment = Comment, cell_line = Cell.line, spp1_profile = SPP1.profile, cisplatine = Cisplatine) %>% 
  select(id, sample_id, cell_line, spp1_profile, cisplatine, comment)

#log2 transform expression data
spp1_salmon_tpm_aggregated = as.matrix(spp1_salmon_tpm_aggregated)
spp1_salmon_tpm_log2 = log2(spp1_salmon_tpm_aggregated + 1)

#export data
expr_data = as.data.frame(spp1_salmon_tpm_log2)
metadata = spp1_seq_meta
save(expr_data, file = "data/expr_data.Rdata")
save(metadata, file = "data/metadata.Rdata")

