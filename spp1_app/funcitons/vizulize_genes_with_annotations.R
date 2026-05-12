visualize_genes_with_annotations <- function(
    expr_data = expr_data_getmm, 
    this_metadata = metadata, 
    genes_of_interest,
    cell_line = NULL,           # Filter: c("SCaBER", "VMCUB1")
    spp1_profile = NULL,        # Filter: c("No", "Recombinant protein", etc.)
    cisplatine = NULL,          # Filter: c("Yes", "No")
    comment = NULL,             # Filter: c("Untreated cells", "treated cells")
    replicate = NULL,           # Filter: c("1", "2", "3")
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    split_by = "replicate",     # Split columns by: "replicate", "cell_line", "cisplatine", "comment", "spp1_profile", or NULL
    show_annotations = c("cell_line", "spp1_profile", "cisplatine", "comment", "replicate")
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
  
  # Filter for genes of interest
  genes_present <- intersect(genes_of_interest, rownames(expr_sub))
  
  if (length(genes_present) == 0) {
    stop("None of the genes of interest found in expression data")
  }
  
  # Filter expression matrix and remove rows with all zeros
  mat <- as.matrix(expr_sub[genes_present, , drop = FALSE])
  mat <- mat[rowSums(mat, na.rm = TRUE) > 0, , drop = FALSE]
  
  # Clean up any NA/Inf values
  mat[is.na(mat)] <- 0
  mat[is.infinite(mat)] <- 0
  
  # Reorder matrix to match metadata
  mat <- mat[, meta_sub$sample_id, drop = FALSE]
  
  # Set up color function for expression
  expr_min <- min(mat, na.rm = TRUE)
  expr_med <- median(mat, na.rm = TRUE)
  expr_max <- max(mat, na.rm = TRUE)
  expr_col_fun <- circlize::colorRamp2(
    c(expr_min, expr_med, expr_max), 
    c("#4DF76F", "black", "#F74D4D")
  )
  
  # Prepare annotation data frame
  annotation_col <- meta_sub %>%
    select(sample_id, cell_line, spp1_profile, cisplatine, comment, replicate) %>%
    column_to_rownames("sample_id")
  
  # Define colors for annotations
  annotation_colors <- list(
    cell_line = c("SCaBER" = "#6F8F72", "VMCUB1" = "#FF7444"),
    spp1_profile = c(
      "No" = "#BFC6C4",
      "Recombinant protein" = "#5A9CB5",
      "Stable overexpression" = "#FACE68",
      "Express by the cell line" = "#FAAC68",
      "Inhibtion with siRNA" = "#FA6868"
    ),
    cisplatine = c("Yes" = "#9E3B3B", "No" = "#BFC6C4"),
    comment = c(
      "Untreated cells" = "#7F55B1",
      "treated cells" = "#F49BAB"
    ),
    replicate = c("1" = "#7C444F", "2" = "#E16A54", "3" = "#F39E60")
  )
  
  # Build annotation object based on show_annotations parameter
  anno_list <- list()
  anno_colors <- list()
  
  if ("replicate" %in% show_annotations) {
    anno_list$replicate <- annotation_col$replicate
    anno_colors$replicate <- annotation_colors$replicate
  }
  if ("comment" %in% show_annotations) {
    anno_list$comment <- annotation_col$comment
    anno_colors$comment <- annotation_colors$comment
  }
  if ("cisplatine" %in% show_annotations) {
    anno_list$cisplatine <- annotation_col$cisplatine
    anno_colors$cisplatine <- annotation_colors$cisplatine
  }
  if ("spp1_profile" %in% show_annotations) {
    anno_list$spp1_profile <- annotation_col$spp1_profile
    anno_colors$spp1_profile <- annotation_colors$spp1_profile
  }
  if ("cell_line" %in% show_annotations) {
    anno_list$cell_line <- annotation_col$cell_line
    anno_colors$cell_line <- annotation_colors$cell_line
  }
  
  # Create annotation object
  ha_col <- do.call(HeatmapAnnotation, c(
    anno_list,
    list(
      col = anno_colors,
      annotation_name_side = "left",
      gap = unit(2, "mm")
    )
  ))
  
  # Determine column split and order based on split_by parameter
  col_split <- NULL
  col_order <- NULL
  
  if (!cluster_columns) {
    if (!is.null(split_by) && split_by %in% colnames(annotation_col)) {
      col_split <- annotation_col[[split_by]]
    }
    col_order <- rownames(annotation_col)
  }
  
  # Create heatmap
  heatmap_obj <- Heatmap(
    mat,
    name = "Expression",
    col = expr_col_fun,
    show_row_names = TRUE,
    row_names_side = "left",
    show_column_names = FALSE,
    column_names_gp = gpar(fontsize = 8),
    row_names_gp = gpar(fontsize = 10),
    column_title = "CD44 and Integrin Expression",
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    column_split = col_split,
    column_order = col_order,
    top_annotation = ha_col
  )
  
  return(list(
    heatmap = heatmap_obj,
    genes_plotted = rownames(mat),
    expression_matrix = mat,
    metadata = meta_sub
  ))
}
