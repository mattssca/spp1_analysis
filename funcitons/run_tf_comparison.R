####################################################################################################
## run_tf_comparison.R
## Modular TF activity comparison across user-defined contrasts
## Returns: list(results_df, heatmap, volcano)
####################################################################################################

run_tf_comparison <- function(
    contrasts,          # list of lists: each must have sb1, sb2, label, context
    acts_mat,           # TF x sample matrix of ULM activity scores
    this_metadata,      # metadata with sample_id and sample_base columns
    tf_set    = NULL,   # optional character vector to restrict to specific TFs
    top_n_hm  = 30,     # number of top variable TFs to show in heatmap
    label_n   = 10,     # number of TFs to label in volcano
    scale_hm  = TRUE    # z-score rows in heatmap to remove baseline differences
) {
  
  library(limma)
  library(tibble)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
  
  ##──────────────────────────────────────────────────────────────────────────
  ## 1. Run differential TF activity for each contrast
  ##──────────────────────────────────────────────────────────────────────────
  
  .diff_activity <- function(sb1, sb2, label, context) {
    samples1 <- this_metadata %>% filter(sample_base == sb1) %>% pull(sample_id)
    samples2 <- this_metadata %>% filter(sample_base == sb2) %>% pull(sample_id)
    
    if (length(samples1) < 2 | length(samples2) < 2)
      stop(paste("Need >= 2 samples per group:", sb1, "vs", sb2))
    
    mat_sub <- acts_mat[, c(samples1, samples2), drop = FALSE]
    if (!is.null(tf_set)) mat_sub <- mat_sub[rownames(mat_sub) %in% tf_set, , drop = FALSE]
    mat_sub <- mat_sub[complete.cases(mat_sub), , drop = FALSE]
    mat_sub <- mat_sub[apply(mat_sub, 1, var) > 0, , drop = FALSE]
    
    group  <- factor(c(rep("grp1", length(samples1)), rep("grp2", length(samples2))))
    design <- model.matrix(~ group)
    fit    <- eBayes(lmFit(mat_sub, design))
    
    topTable(fit, coef = 2, number = Inf, sort.by = "P") %>%
      rownames_to_column("TF") %>%
      mutate(contrast = label, context = context,
             sb1 = sb1, sb2 = sb2)
  }
  
  results_df <- bind_rows(
    lapply(contrasts, function(c) {
      .diff_activity(c$sb1, c$sb2, c$label, c$context)
    })
  )
  
  ##──────────────────────────────────────────────────────────────────────────
  ## 2. Volcano plot (one facet per contrast)
  ##──────────────────────────────────────────────────────────────────────────
  
  plot_df <- results_df %>%
    mutate(
      sig = case_when(
        adj.P.Val < 0.05 & logFC > 0 ~ "Up",
        adj.P.Val < 0.05 & logFC < 0 ~ "Down",
        TRUE ~ "NS"
      )
    ) %>%
    group_by(contrast) %>%
    mutate(
      rank_p = rank(P.Value),
      label  = ifelse(sig != "NS" & rank_p <= label_n, TF, "")
    ) %>%
    ungroup()
  
  volcano <- ggplot(plot_df,
                    aes(x = logFC, y = -log10(P.Value),
                        color = sig, label = label)) +
    geom_point(alpha = 0.7, size = 1.8) +
    geom_text_repel(size = 3, max.overlaps = 20) +
    scale_color_manual(values = c("Up" = "#D73027",
                                  "Down" = "#4575B4",
                                  "NS"   = "grey70")) +
    geom_vline(xintercept = 0,         linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    facet_wrap(~ contrast, scales = "free") +
    labs(x = "log2 fold change (TF activity)",
         y = "-log10(p-value)",
         color = "Direction") +
    theme_bw() +
    theme(strip.text = element_text(size = 8))
  
  ##──────────────────────────────────────────────────────────────────────────
  ## 3. Heatmap: top variable TFs across samples in selected contrasts
  ##──────────────────────────────────────────────────────────────────────────
  
  all_samples <- unique(unlist(
    lapply(contrasts, function(c) {
      this_metadata %>%
        filter(sample_base %in% c(c$sb1, c$sb2)) %>%
        pull(sample_id)
    })
  ))
  
  # Top N most variable TFs within this sample set
  candidate_tfs <- if (!is.null(tf_set)) {
    rownames(acts_mat)[rownames(acts_mat) %in% tf_set]
  } else {
    rownames(acts_mat)
  }
  
  top_tfs <- apply(acts_mat[candidate_tfs, all_samples, drop = FALSE], 1, var) %>%
    sort(decreasing = TRUE) %>%
    head(top_n_hm) %>%
    names()
  
  acts_sub <- acts_mat[top_tfs, all_samples, drop = FALSE]
  if (scale_hm) acts_sub <- t(scale(t(acts_sub)))
  
  meta_sub <- this_metadata %>%
    filter(sample_id %in% all_samples) %>%
    column_to_rownames("sample_id")
  
  annotation_colors <- list(
    cell_line = c("SCaBER" = "#E41A1C", "VMCUB1" = "#377EB8"),
    spp1_profile = c(
      "No"                       = "#999999",
      "Recombinant protein"      = "#FF7F00",
      "Stable overexpression"    = "#984EA3",
      "Express by the cell line" = "#4DAF4A",
      "Inhibtion with siRNA"     = "#E6AB02"
    ),
    cisplatine = c("Yes" = "#D62728", "No" = "#BCBCBC"),
    replicate  = c("1" = "#F26076", "2" = "#FF9760", "3" = "#458B73")
  )
  
  ha_col <- HeatmapAnnotation(
    replicate    = meta_sub[colnames(acts_sub), "replicate"],
    cisplatine   = meta_sub[colnames(acts_sub), "cisplatine"],
    spp1_profile = meta_sub[colnames(acts_sub), "spp1_profile"],
    cell_line    = meta_sub[colnames(acts_sub), "cell_line"],
    col          = annotation_colors,
    annotation_name_side = "left",
    gap          = unit(2, "mm")
  )
  
  score_lim     <- max(abs(acts_sub), na.rm = TRUE)
  score_col_fun <- colorRamp2(c(-score_lim, 0, score_lim),
                              c("#4575B4", "white", "#D73027"))
  
  heatmap <- Heatmap(
    acts_sub,
    name              = ifelse(scale_hm, "TF activity\n(z-score)", "TF activity\n(ULM score)"),
    col               = score_col_fun,
    show_column_names = FALSE,
    show_row_names    = TRUE,
    row_names_side    = "left",
    row_names_gp      = gpar(fontsize = 7),
    cluster_rows      = TRUE,
    cluster_columns   = FALSE,
    column_split      = meta_sub[colnames(acts_sub), "sample_base"],
    column_title_gp   = gpar(fontsize = 7),
    column_title_rot  = 90,
    top_annotation    = ha_col
  )
  
  return(list(
    results_df = results_df,
    heatmap    = heatmap,
    volcano    = volcano
  ))
}
