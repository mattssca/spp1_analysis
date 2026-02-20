#load data
load("C:/Users/matts/Desktop/projects/spp1/new_analysis/data/expr_data.Rdata")
load("C:/Users/matts/Desktop/projects/spp1/new_analysis/data/metadata.Rdata")

#libraries
library(LundTaxR)

#classify samples
predicted = LundTaxR::classify_samples(this_data = expr_data, 
                                       log_transform = FALSE, 
                                       include_data = TRUE)
#export data
save(predicted, file = "data/predicted.Rdata")
