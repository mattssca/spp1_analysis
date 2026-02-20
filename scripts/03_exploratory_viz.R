# Load required libraries
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

####################################################################################################
## Load Data
####################################################################################################
load("C:/Users/matts/Desktop/projects/spp1/new_analysis/data/expr_data.Rdata")
load("C:/Users/matts/Desktop/projects/spp1/new_analysis/data/metadata.Rdata")
load("C:/Users/matts/Desktop/projects/spp1/new_analysis/data/predicted.Rdata")

####################################################################################################
## Prepare Metadata and Annotations
####################################################################################################
# Add replicate information and create sorting variables
metadata <- metadata %>%
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

# Prepare annotation data frame
annotation_col <- metadata %>%
  select(sample_id, cell_line, spp1_profile, cisplatine, comment, replicate) %>%
  column_to_rownames("sample_id")

# Define colors for annotations
annotation_colors <- list(
  cell_line = c("SCaBER" = "#E41A1C", "VMCUB1" = "#377EB8"),
  spp1_profile = c(
    "No" = "#999999",
    "Recombinant protein" = "#FF7F00",
    "Stable overexpression" = "#984EA3",
    "Express by the cell line" = "#4DAF4A",
    "Inhibtion with siRNA" = "#E6AB02"
  ),
  cisplatine = c("Yes" = "#D62728", "No" = "#BCBCBC"),
  comment = c(
    "Untreated cells" = "#1F77B4",
    "treated cells" = "#FF7F0E"
  ),
  replicate = c("1" = "#F26076", "2" = "#FF9760", "3" = "#458B73")
)

# Create shared annotation object
ha_col <- HeatmapAnnotation(
  replicate = annotation_col$replicate,
  comment = annotation_col$comment,
  cisplatine = annotation_col$cisplatine,
  spp1_profile = annotation_col$spp1_profile,
  cell_line = annotation_col$cell_line,
  col = annotation_colors,
  annotation_name_side = "left",
  gap = unit(2, "mm")
)

####################################################################################################
## Heatmap 1: Lund Taxonomy Prediction Scores
####################################################################################################
# Get prediction scores and reorder to match metadata
pred_scores <- predicted$subtype_scores[metadata$sample_id, c("Uro", "GU", "BaSq", "Mes", "ScNE", "prediction_delta_collapsed")]


# Define color function
col_fun_pred <- colorRamp2(c(0, 0.5, 1), c("white", "steelblue", "navy"))

# Create prediction scores heatmap
ht1 <- Heatmap(
  t(pred_scores),
  name = "Prediction\nScore",
  col = col_fun_pred,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_split = annotation_col$replicate,
  column_order = rownames(annotation_col),
  show_column_names = TRUE,
  show_row_names = TRUE,
  top_annotation = ha_col,
  column_title = "Lund Taxonomy Prediction Scores",
  row_names_side = "left"
)

####################################################################################################
## Heatmap 2: SPP1 Isoform Expression
####################################################################################################
# Load transcript-level data
transcript_tpm <- read_tsv("C:/Users/matts/Desktop/DATA/SPP1/RNA-seq/salmon.merged.transcript_tpm.tsv")
tx2gene <- read_tsv("C:/Users/matts/Desktop/DATA/SPP1/RNA-seq/tx2gene.tsv", 
                    col_names = c("transcript_id", "gene_id", "gene_name"))

# Filter for SPP1 transcripts
spp1_transcripts <- tx2gene %>%
  filter(gene_name == "SPP1" | gene_id == "ENSG00000118785")

spp1_tpm <- transcript_tpm %>%
  filter(tx %in% spp1_transcripts$transcript_id) %>%
  select(-gene_id) %>%
  column_to_rownames("tx")

# Add transcript labels and order by metadata
spp1_tpm <- spp1_tpm %>%
  rownames_to_column("tx") %>%
  left_join(tx2gene %>% select(transcript_id, gene_name), 
            by = c("tx" = "transcript_id")) %>%
  mutate(label = paste0(tx, " (", gene_name, ")")) %>%
  select(-tx, -gene_name) %>%
  column_to_rownames("label") %>%
  select(all_of(metadata$sample_id))

# Log2 transform
spp1_log_tpm <- log2(spp1_tpm + 1)

# Define color function
col_fun_tpm <- colorRamp2(c(0, 5, 10), c("white", "orange", "red"))

# Create SPP1 isoform heatmap
ht2 <- Heatmap(
  as.matrix(spp1_log_tpm),
  name = "log2(TPM+1)",
  col = col_fun_tpm,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_split = annotation_col$replicate,
  column_order = rownames(annotation_col),
  show_column_names = TRUE,
  show_row_names = TRUE,
  column_title = "SPP1 Isoform Expression",
  row_names_side = "left",
  column_names_gp = gpar(fontsize = 8),
  row_names_gp = gpar(fontsize = 9),
  heatmap_legend_param = list(
    title = "Expression",
    direction = "vertical"
  )
)

####################################################################################################
## Draw Combined Heatmap
####################################################################################################
draw(ht1 %v% ht2, column_title = "SPP1 Analysis: Lund Taxonomy & Isoform Expression")


####################################################################################################
## Replicate Correlation Visualization - Dot Plot
####################################################################################################

# Prepare data for plotting
replicate_plot_data <- list()

for(iso in rownames(spp1_tpm_expressed)) {
  # Extract TPM values as a vector
  tpm_values <- as.numeric(spp1_tpm_expressed[iso, metadata$sample_id])
  
  iso_data <- metadata %>%
    as.data.frame() %>%
    mutate(TPM = tpm_values) %>%
    select(sample_id, cell_line, spp1_profile, cisplatine, comment, replicate, TPM) %>%
    mutate(condition = paste(cell_line, spp1_profile, cisplatine, comment, sep = "_"))
  
  replicate_plot_data[[iso]] <- iso_data
}

# Combine all isoforms
all_replicate_data <- bind_rows(replicate_plot_data, .id = "isoform")

# Create pairwise comparisons for dot plots
replicate_pairs <- all_replicate_data %>%
  select(isoform, condition, replicate, TPM) %>%
  pivot_wider(names_from = replicate, values_from = TPM, names_prefix = "Rep_") %>%
  as.data.frame()

# Calculate correlations for each isoform
cor_labels <- replicate_pairs %>%
  group_by(isoform) %>%
  summarise(
    cor_1_2 = cor(Rep_1, Rep_2, use = "complete.obs"),
    cor_1_3 = cor(Rep_1, Rep_3, use = "complete.obs"),
    cor_2_3 = cor(Rep_2, Rep_3, use = "complete.obs")
  )

# Plot Rep1 vs Rep2
p1 <- ggplot(replicate_pairs, aes(x = Rep_1, y = Rep_2)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  geom_smooth(method = "lm", se = TRUE, color = "blue", alpha = 0.2) +
  facet_wrap(~isoform, scales = "free", ncol = 3) +
  theme_bw() +
  labs(
    title = "Replicate 1 vs Replicate 2",
    x = "Replicate 1 (TPM)",
    y = "Replicate 2 (TPM)"
  ) +
  theme(strip.text = element_text(size = 8))

# Add correlation labels
p1 <- p1 + geom_text(
  data = cor_labels,
  aes(x = -Inf, y = Inf, label = paste0("r = ", round(cor_1_2, 3))),
  hjust = -0.1, vjust = 1.5, size = 3, inherit.aes = FALSE
)

print(p1)

# Plot Rep1 vs Rep3
p2 <- ggplot(replicate_pairs, aes(x = Rep_1, y = Rep_3)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  geom_smooth(method = "lm", se = TRUE, color = "blue", alpha = 0.2) +
  facet_wrap(~isoform, scales = "free", ncol = 3) +
  theme_bw() +
  labs(
    title = "Replicate 1 vs Replicate 3",
    x = "Replicate 1 (TPM)",
    y = "Replicate 3 (TPM)"
  ) +
  theme(strip.text = element_text(size = 8))

p2 <- p2 + geom_text(
  data = cor_labels,
  aes(x = -Inf, y = Inf, label = paste0("r = ", round(cor_1_3, 3))),
  hjust = -0.1, vjust = 1.5, size = 3, inherit.aes = FALSE
)

print(p2)

# Plot Rep2 vs Rep3
p3 <- ggplot(replicate_pairs, aes(x = Rep_2, y = Rep_3)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  geom_smooth(method = "lm", se = TRUE, color = "blue", alpha = 0.2) +
  facet_wrap(~isoform, scales = "free", ncol = 3) +
  theme_bw() +
  labs(
    title = "Replicate 2 vs Replicate 3",
    x = "Replicate 2 (TPM)",
    y = "Replicate 3 (TPM)"
  ) +
  theme(strip.text = element_text(size = 8))

p3 <- p3 + geom_text(
  data = cor_labels,
  aes(x = -Inf, y = Inf, label = paste0("r = ", round(cor_2_3, 3))),
  hjust = -0.1, vjust = 1.5, size = 3, inherit.aes = FALSE
)

print(p3)
