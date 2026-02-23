#load data
load("C:/Users/matts/Desktop/GIT_REPOS/spp1_analysis/data/expr_data_getmm.Rdata")
load("C:/Users/matts/Desktop/GIT_REPOS/spp1_analysis/data/metadata.Rdata")

#source function
source('C:/Users/matts/Desktop/GIT_REPOS/spp1_analysis/funcitons/check_spp1_correlation_volcano.R')

#load libraries
library(ggplot2)
library(dplyr)
library(ggrepel)
library(patchwork)

#subset sample base info
metadata$sample_base <- sub("_[0-9]+$", "", metadata$sample_id)

#get unique sample abses
sample_bases <- unique(metadata$sample_base)
cor_results_list <- setNames(vector("list", length(sample_bases)), sample_bases)

#run on all sample base
for (sb in sample_bases) {
  cor_results_list[[sb]] <- check_spp1_correlation_volcano(expr_data_getmm, metadata, sb)
}

#combine all plots into one grid (2 rows x 5 columns for 10 plots, adjust as needed)
combined_plot <- 
  cor_results_list[["SCABER_CTRL"]]$plot +
  cor_results_list[["VMCUB1_CTRL"]]$plot +
  cor_results_list[["SCABER_OESPP1"]]$plot +
  cor_results_list[["VMCUB1_SISPP1"]]$plot +
  cor_results_list[["SCABER_RSPP1"]]$plot +
  cor_results_list[["SCABER_CISPL"]]$plot +
  cor_results_list[["VMCUB1_CISPL"]]$plot +
  cor_results_list[["SCABER_OESPP1_CISPL"]]$plot +
  cor_results_list[["SCABER_RSPP1_CISPL"]]$plot +
  cor_results_list[["VMCUB1_SISPP1_CISPL"]]$plot +
  plot_layout(ncol = 5, guides = "collect") & 
  theme(legend.position = "bottom")

# Save to PDF (or PNG)
ggsave("all_volcano_plots.pdf", combined_plot, width = 20, height = 8)

top_bottom_df <- do.call(
  rbind,
  lapply(names(cor_results_list), function(sb) {
    df <- cor_results_list[[sb]]$results
    # Top 5 positively correlated
    top5 <- df %>% arrange(-correlation) %>% head(5)
    # Top 5 negatively correlated
    bottom5 <- df %>% arrange(correlation) %>% head(5)
    # Combine and annotate
    rbind(top5, bottom5) %>% mutate(sample_base = sb)
  })
)
