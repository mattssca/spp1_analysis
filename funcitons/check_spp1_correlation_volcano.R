check_spp1_correlation_volcano <- function(expr_data = expr_data_getmm, 
                                           this_metadata = metadata, 
                                           sample_base_value, 
                                           method = "pearson",
                                           cor_cutoff = 0.7, 
                                           pval_cutoff = 0.05, 
                                           top_n_labels = 5){
  
  # Get sample IDs for the specified sample base
  samples <- this_metadata %>%
    filter(sample_base == sample_base_value) %>%
    pull(sample_id)
  
  print(samples)
  
  # Subset expression data for these samples
  expr_sub <- expr_data[, samples, drop = FALSE]
  
  # Check if SPP1 is present
  if (!"SPP1" %in% rownames(expr_sub)) {
    stop("SPP1 not found in expression data.")
  }
  
  spp1_expr <- as.numeric(expr_sub["SPP1", ])
  
  # Compute correlation and p-value for each gene with SPP1
  cor_list <- apply(expr_sub, 1, function(x) {
    res <- tryCatch(
      cor.test(as.numeric(x), spp1_expr, method = method),
      error = function(e) return(list(estimate = NA, p.value = NA))
    )
    c(cor = as.numeric(res$estimate), pval = as.numeric(res$p.value))
  })
  
  cor_df <- data.frame(
    gene = rownames(expr_sub),
    correlation = as.numeric(cor_list["cor", ]),
    pval = as.numeric(cor_list["pval", ])
  )
  
  # Exclude SPP1 from results
  cor_df <- cor_df %>% filter(gene != "SPP1")
  
  # Add significance column
  cor_df <- cor_df %>%
    mutate(
      significance = case_when(
        pval < pval_cutoff & correlation > cor_cutoff ~ "Positively Correlated",
        pval < pval_cutoff & correlation < -cor_cutoff ~ "Negatively Correlated",
        TRUE ~ "Not Significant"
      )
    )
  
  # Select top N positive and negative significant genes for labeling
  label_df <- bind_rows(
    cor_df %>% filter(significance == "Positively Correlated") %>% arrange(-correlation) %>% head(top_n_labels),
    cor_df %>% filter(significance == "Negatively Correlated") %>% arrange(correlation) %>% head(top_n_labels)
  )
  
  # Volcano plot
  p <- ggplot(cor_df, aes(x = correlation, y = -log10(pval), color = significance)) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = c(
      "Positively Correlated" = "red",
      "Negatively Correlated" = "blue",
      "Not Significant" = "grey"
    )) +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(label = gene),
      size = 3,
      max.overlaps = Inf
    ) +
    theme_minimal() +
    labs(
      title = paste0(sample_base_value),
      x = "Correlation with SPP1",
      y = "-log10(p-value)"
    ) +
    theme(legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5))
  
  print(p)
  return(list(results = cor_df, plot = p))
}