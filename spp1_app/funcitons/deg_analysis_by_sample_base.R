deg_analysis_by_sample_base <- function(
    expr_data = expr_data_getmm, 
    this_metadata = metadata, 
    sample_base1, 
    sample_base2, 
    top_n = 20,
    heatmap_genes = "significant",
    max_sig_genes = 20
) {

  # Color palettes for sample base and replicate
  sample_base_colors <- c("#09637E", "#FF5B5B")
  replicate_colors <- c("#3D45AA", "#F8843F", "#427A43")
  
  # Filter metadata for the selected sample bases (all replicates)
  meta_sub <- this_metadata %>%
    filter(sample_base %in% c(sample_base1, sample_base2))
  
  # Get sample IDs for each group
  samples1 <- meta_sub %>% filter(sample_base == sample_base1) %>% pull(sample_id)
  samples2 <- meta_sub %>% filter(sample_base == sample_base2) %>% pull(sample_id)
  samples <- c(samples1, samples2)
  
  # Check that both groups have at least two samples
  if(length(samples1) < 2 | length(samples2) < 2) {
    stop("Both groups must have at least two samples for DEG analysis.")
  }
  
  # Subset expression data
  expr_sub <- expr_data[, samples, drop = FALSE]
  
  # Remove genes with zero variance across selected samples
  expr_sub <- expr_sub[apply(expr_sub, 1, function(x) var(as.numeric(x)) > 0), , drop = FALSE]
  
  # Create design matrix
  group <- factor(c(rep("group1", length(samples1)), rep("group2", length(samples2))))
  design <- model.matrix(~ group)
  
  # Run limma
  fit <- lmFit(expr_sub, design)
  fit <- eBayes(fit)
  res <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
  res <- rownames_to_column(res, "gene")
  
  # Top and bottom DEGs
  top_deg <- res %>% arrange(P.Value) %>% head(top_n)
  bottom_deg <- res %>% arrange(desc(P.Value)) %>% head(top_n)
  
  # Decide which genes to plot in the heatmap
  genes_to_plot <- NULL
  if (is.character(heatmap_genes) && length(heatmap_genes) == 1) {
    if (heatmap_genes == "significant") {
      sig_genes <- res %>% filter(adj.P.Val < 0.05) %>% arrange(P.Value) %>% pull(gene)
      genes_to_plot <- head(sig_genes, max_sig_genes)
    } else if (heatmap_genes == "top_bottom") {
      genes_to_plot <- unique(c(top_deg$gene, bottom_deg$gene))
    }
  } else if (is.character(heatmap_genes)) {
    genes_to_plot <- intersect(heatmap_genes, rownames(expr_sub))
  }
  
  heatmap_obj <- NULL
  if (!is.null(genes_to_plot) && length(genes_to_plot) > 1) {
    ann_df <- res %>%
      filter(gene %in% genes_to_plot) %>%
      select(gene, logFC, adj.P.Val) %>%
      column_to_rownames("gene")
    ann_df$adj.P.Val <- signif(ann_df$adj.P.Val, 2)
    
    mat <- as.matrix(expr_sub[genes_to_plot, , drop = FALSE])
    
    # Dynamically set color ramps based on data
    expr_min <- min(mat, na.rm = TRUE)
    expr_med <- median(mat, na.rm = TRUE)
    expr_max <- max(mat, na.rm = TRUE)
    expr_col_fun <- circlize::colorRamp2(c(expr_min, expr_med, expr_max), c("#4DF76F", "black", "#F74D4D"))
    
    logfc_min <- min(ann_df$logFC, na.rm = TRUE)
    logfc_max <- max(ann_df$logFC, na.rm = TRUE)
    logfc_col_fun <- circlize::colorRamp2(c(logfc_min, 0, logfc_max), c("#B2182B", "white", "#2166AC"))
    
    pval_min <- min(ann_df$adj.P.Val, na.rm = TRUE)
    pval_max <- max(ann_df$adj.P.Val, na.rm = TRUE)
    pval_col_fun <- circlize::colorRamp2(c(pval_min, 0.05, pval_max), c("#B2182B", "white", "#2166AC"))
    
    # Row annotation (right side, gene names on left)
    row_ha <- rowAnnotation(
      log2FC = anno_simple(ann_df$logFC, col = logfc_col_fun),
      adj.P.Val = anno_simple(ann_df$adj.P.Val, col = pval_col_fun),
      annotation_name_side = "bottom"
    )
    
    # Column annotation (sample base and replicate)
    col_sample_bases <- this_metadata %>%
      filter(sample_id %in% colnames(expr_sub)) %>%
      arrange(match(sample_id, colnames(expr_sub)))
    sample_base_colormap <- setNames(sample_base_colors, c(sample_base1, sample_base2))
    replicate_levels <- sort(unique(col_sample_bases$replicate))
    replicate_colormap <- setNames(replicate_colors[seq_along(replicate_levels)], replicate_levels)
    
    col_ha <- HeatmapAnnotation(
      SampleBase = col_sample_bases$sample_base,
      Replicate = as.factor(col_sample_bases$replicate),
      col = list(
        SampleBase = sample_base_colormap,
        Replicate = replicate_colormap
      ),
      annotation_name_side = "left"
    )
    
    heatmap_obj <- Heatmap(
      mat,
      name = "Expression",
      col = expr_col_fun,
      show_heatmap_legend = FALSE,
      show_row_names = TRUE,
      row_names_side = "left",
      show_column_names = FALSE,
      column_title = paste(sample_base1, "vs", sample_base2),
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      right_annotation = row_ha,
      top_annotation = col_ha
    )
  }
  
  return(list(
    deg_results = res,
    top_deg = top_deg,
    bottom_deg = bottom_deg,
    heatmap = heatmap_obj
  ))
}