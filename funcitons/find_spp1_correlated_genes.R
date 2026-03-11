find_spp1_correlated_genes <- function(
    expr_data = expr_data_getmm,
    this_metadata = metadata,
    cell_line = NULL,
    spp1_profile = NULL,
    cisplatine = NULL,
    comment = NULL,
    replicate = NULL,
    cor_method = "pearson",     # "pearson", "spearman", or "kendall"
    min_cor = 0.5,              # Minimum absolute correlation threshold
    p_value_threshold = 0.05,   # P-value threshold for significance
    top_n = NULL,               # Optional: return only top N genes for each direction
    label_top_n = 10            # Number of top genes to label on volcano plot
) {
  
  # Start with all metadata
  meta_sub <- this_metadata
  
  # Apply filters if specified
  if (!is.null(cell_line)) {
    meta_sub <- meta_sub %>% filter(cell_line %in% !!cell_line)
  }
  if (!is.null(spp1_profile)) {
    meta_sub <- meta_sub %>% filter(spp1_profile %in% !!spp1_profile)
  }
  if (!is.null(cisplatine)) {
    meta_sub <- meta_sub %>% filter(cisplatine %in% !!cisplatine)
  }
  if (!is.null(comment)) {
    meta_sub <- meta_sub %>% filter(comment %in% !!comment)
  }
  if (!is.null(replicate)) {
    meta_sub <- meta_sub %>% filter(replicate %in% !!replicate)
  }
  
  # Check if any samples remain after filtering
  if (nrow(meta_sub) == 0) {
    stop("No samples match the specified filters")
  }
  
  # Get sample IDs
  samples <- meta_sub$sample_id
  expr_sub <- expr_data[, samples, drop = FALSE]
  
  # Check if SPP1 exists in the data
  if (!"SPP1" %in% rownames(expr_sub)) {
    stop("SPP1 not found in expression data")
  }
  
  # Get SPP1 expression
  spp1_expr <- as.numeric(expr_sub["SPP1", ])
  
  # Remove genes with zero variance
  gene_vars <- apply(expr_sub, 1, var, na.rm = TRUE)
  expr_sub <- expr_sub[gene_vars > 0, , drop = FALSE]
  
  # Calculate correlations with SPP1 for all genes
  cor_results <- data.frame(
    gene = rownames(expr_sub),
    correlation = NA,
    p_value = NA,
    stringsAsFactors = FALSE
  )
  
  for (i in 1:nrow(expr_sub)) {
    gene_expr <- as.numeric(expr_sub[i, ])
    
    # Skip if it's SPP1 itself
    if (rownames(expr_sub)[i] == "SPP1") {
      cor_results$correlation[i] <- 1
      cor_results$p_value[i] <- 0
      next
    }
    
    # Calculate correlation
    cor_test <- cor.test(spp1_expr, gene_expr, method = cor_method)
    cor_results$correlation[i] <- cor_test$estimate
    cor_results$p_value[i] <- cor_test$p.value
  }
  
  # Add adjusted p-values and -log10(p-value)
  cor_results$adj_p_value <- p.adjust(cor_results$p_value, method = "BH")
  cor_results$neg_log10_p <- -log10(cor_results$p_value)
  cor_results$neg_log10_adj_p <- -log10(cor_results$adj_p_value)
  
  # Classify genes for volcano plot
  cor_results <- cor_results %>%
    mutate(
      significance = case_when(
        adj_p_value < p_value_threshold & correlation > min_cor ~ "Positive",
        adj_p_value < p_value_threshold & correlation < -min_cor ~ "Negative",
        TRUE ~ "Not Significant"
      )
    )
  
  # Filter for significant correlations
  significant <- cor_results %>%
    filter(adj_p_value < p_value_threshold,
           abs(correlation) >= min_cor,
           gene != "SPP1") %>%
    arrange(desc(abs(correlation)))
  
  # Split into positive and negative correlations
  positive_cor <- significant %>%
    filter(correlation > 0) %>%
    arrange(desc(correlation))
  
  negative_cor <- significant %>%
    filter(correlation < 0) %>%
    arrange(correlation)
  
  # Identify genes to label on volcano plot
  top_positive <- positive_cor %>% head(label_top_n) %>% pull(gene)
  top_negative <- negative_cor %>% head(label_top_n) %>% pull(gene)
  genes_to_label <- c(top_positive, top_negative)
  
  cor_results <- cor_results %>%
    mutate(label = ifelse(gene %in% genes_to_label, gene, ""))
  
  # Create volcano plot
  volcano_plot <- ggplot(cor_results, aes(x = correlation, y = neg_log10_adj_p)) +
    geom_point(aes(color = significance), alpha = 0.6, size = 2) +
    scale_color_manual(
      values = c(
        "Positive" = "#D62728",
        "Negative" = "#1F77B4",
        "Not Significant" = "#999999"
      )
    ) +
    geom_vline(xintercept = c(-min_cor, min_cor), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(p_value_threshold), linetype = "dashed", color = "black") +
    geom_text_repel(
      aes(label = label),
      size = 3,
      max.overlaps = 20,
      box.padding = 0.5,
      point.padding = 0.2
    ) +
    labs(
      title = "Volcano Plot: Correlation with SPP1",
      subtitle = paste0("Method: ", cor_method, " | Samples: ", length(samples)),
      x = "Correlation with SPP1",
      y = expression(-log[10](adjusted~p-value)),
      color = "Significance"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10)
    )
  
  # Apply top_n if specified
  if (!is.null(top_n)) {
    positive_cor <- positive_cor %>% head(top_n)
    negative_cor <- negative_cor %>% head(top_n)
  }
  
  # Create summary statistics
  summary_stats <- list(
    n_samples = length(samples),
    n_genes_tested = nrow(cor_results) - 1,  # Exclude SPP1 itself
    n_positive = nrow(positive_cor),
    n_negative = nrow(negative_cor),
    mean_spp1_expression = mean(spp1_expr, na.rm = TRUE),
    sd_spp1_expression = sd(spp1_expr, na.rm = TRUE),
    cor_method = cor_method,
    min_cor_threshold = min_cor,
    p_value_threshold = p_value_threshold
  )
  
  return(list(
    positive_correlation = positive_cor,
    negative_correlation = negative_cor,
    all_correlations = cor_results,
    volcano_plot = volcano_plot,
    summary = summary_stats,
    metadata_used = meta_sub,
    spp1_expression = spp1_expr
  ))
}

