#load data
load("C:/Users/matts/Desktop/GIT_REPOS/spp1_analysis/data/expr_data_getmm.Rdata")
load("C:/Users/matts/Desktop/GIT_REPOS/spp1_analysis/data/metadata.Rdata")

#source function
source('C:/Users/matts/Desktop/GIT_REPOS/spp1_analysis/funcitons/deg_analysis.R')

#load packages
library(limma)
library(ComplexHeatmap)
library(dplyr)
library(tibble)
library(circlize)

#subset sample base info
metadata$sample_base <- sub("_[0-9]+$", "", metadata$sample_id)

#run DEG analysis
hm1 <- deg_analysis_by_sample_base(sample_base1 = "VMCUB1_CISPL", sample_base2 =  "VMCUB1_SISPP1_CISPL")
hm2 <- deg_analysis_by_sample_base(sample_base1 = "VMCUB1_CTRL", sample_base2 =  "VMCUB1_CISPL")
hm3 <- deg_analysis_by_sample_base(sample_base1 = "VMCUB1_CTRL", sample_base2 =  "VMCUB1_SISPP1")
hm4 <- deg_analysis_by_sample_base(sample_base1 = "SCABER_CISPL", sample_base2 =  "SCABER_OESPP1_CISPL")
hm5 <- deg_analysis_by_sample_base(sample_base1 = "SCABER_CTRL", sample_base2 =  "SCABER_CISPL")
hm6 <- deg_analysis_by_sample_base(sample_base1 = "SCABER_CTRL", sample_base2 =  "SCABER_OESPP1")
hm7 <- deg_analysis_by_sample_base(sample_base1 = "SCABER_CTRL", sample_base2 =  "SCABER_RSPP1")
hm8 <- deg_analysis_by_sample_base(sample_base1 = "SCABER_CISPL", sample_base2 =  "SCABER_RSPP1_CISPL")

draw(hm1$heatmap)
draw(hm2$heatmap)
draw(hm3$heatmap)
draw(hm4$heatmap)
draw(hm5$heatmap)
draw(hm6$heatmap)
draw(hm7$heatmap)
draw(hm8$heatmap)

for (i in 1:8) {
  png(sprintf("heatmap_%d.png", i), width = 8, height = 8, units = "in", res = 300)
  draw(get(paste0("hm", i))$heatmap)
  dev.off()
}

deg_list <- list(
  hm1 = hm1$deg_results,
  hm2 = hm2$deg_results,
  hm3 = hm3$deg_results,
  hm4 = hm4$deg_results,
  hm5 = hm5$deg_results,
  hm6 = hm6$deg_results,
  hm7 = hm7$deg_results,
  hm8 = hm8$deg_results
)
