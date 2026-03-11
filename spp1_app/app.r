# Load required libraries
library(shiny)
library(shinydashboard)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)

# ============================================================================
# BUNDLED DATA
# ============================================================================

# Load data from the data folder
load("data/expr_data_getmm.Rdata")  # Should contain expr_data_getmm
load("data/metadata.Rdata")   # Should contain metadata

# Verify data is loaded
if (!exists("expr_data_getmm")) {
  stop("expr_data_getmm not found in data/expr_data.Rdata")
}
if (!exists("metadata")) {
  stop("metadata not found in data/metadata.Rdata")
}

# ============================================================================
# ANALYSIS FUNCTION
# ============================================================================

visualize_genes_with_annotations <- function(
    expr_data = expr_data_getmm, 
    this_metadata = metadata, 
    genes_of_interest,
    cell_line = NULL,
    spp1_profile = NULL,
    cisplatine = NULL,
    comment = NULL,
    replicate = NULL,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    split_by = "replicate",
    show_annotations = c("cell_line", "spp1_profile", "cisplatine", "comment", "replicate")
) {
  meta_sub <- this_metadata
  
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
  
  if (nrow(meta_sub) == 0) {
    stop("No samples match the specified filters")
  }
  
  samples <- meta_sub$sample_id
  expr_sub <- expr_data[, samples, drop = FALSE]
  
  genes_present <- intersect(genes_of_interest, rownames(expr_sub))
  
  if (length(genes_present) == 0) {
    stop("None of the genes of interest found in expression data")
  }
  
  mat <- as.matrix(expr_sub[genes_present, , drop = FALSE])
  mat <- mat[rowSums(mat, na.rm = TRUE) > 0, , drop = FALSE]
  
  mat[is.na(mat)] <- 0
  mat[is.infinite(mat)] <- 0
  
  mat <- mat[, meta_sub$sample_id, drop = FALSE]
  
  expr_min <- min(mat, na.rm = TRUE)
  expr_med <- median(mat, na.rm = TRUE)
  expr_max <- max(mat, na.rm = TRUE)
  expr_col_fun <- circlize::colorRamp2(
    c(expr_min, expr_med, expr_max), 
    c("#4DF76F", "black", "#F74D4D")
  )
  
  annotation_col <- meta_sub %>%
    select(sample_id, cell_line, spp1_profile, cisplatine, comment, replicate) %>%
    column_to_rownames("sample_id")
  
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
  
  ha_col <- do.call(HeatmapAnnotation, c(
    anno_list,
    list(
      col = anno_colors,
      annotation_name_side = "left",
      gap = unit(2, "mm")
    )
  ))
  
  col_split <- NULL
  col_order <- NULL
  
  if (!cluster_columns) {
    if (!is.null(split_by) && split_by %in% colnames(annotation_col)) {
      col_split <- annotation_col[[split_by]]
    }
    col_order <- rownames(annotation_col)
  }
  
  heatmap_obj <- Heatmap(
    mat,
    name = "Expression",
    col = expr_col_fun,
    show_row_names = TRUE,
    row_names_side = "left",
    show_column_names = FALSE,
    column_names_gp = gpar(fontsize = 8),
    row_names_gp = gpar(fontsize = 10),
    column_title = "Gene Expression",
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

# ============================================================================
# SHINY UI
# ============================================================================

ui <- dashboardPage(
  skin = "blue",
  
  dashboardHeader(title = "Gene Expression Viewer"),
  
  dashboardSidebar(disable = TRUE),
  
  dashboardBody(
        fluidRow(
          box(
            title = "Visualization Parameters",
            width = 3,
            solidHeader = TRUE,
            status = "primary",
            
            h4("Genes Input"),
            textAreaInput("genes_input", "Genes of Interest (one per line, default is CD44 and ITGa/b):",
                          value = "CD44\nITGA1\nITGA10\nITGA11\nITGA2\nITGA2B\nITGA3\nITGA4\nITGA5\nITGA6\nITGA7\nITGA8\nITGA9\nITGAD\nITGAE\nITGAL\nITGAM\nITGAV\nITGAX\nITGB1\nITGB1BP1\nITGB1BP2\nITGB2\nITGB3\nITGB3BP\nITGB4\nITGB5\nITGB6\nITGB7\nITGB8\nITGBL1",
                          rows = 8),
            
            h4("Filters:"),
            checkboxGroupInput("filter_cell_line", "Cell Line:",
                               choices = NULL),
            checkboxGroupInput("filter_spp1_profile", "SPP1 Profile:",
                               choices = NULL),
            checkboxGroupInput("filter_cisplatine", "Cisplatine:",
                               choices = NULL),
            checkboxGroupInput("filter_comment", "Comment:",
                               choices = NULL),
            checkboxGroupInput("filter_replicate", "Replicate:",
                               choices = NULL),
            
            h4("Display Options:"),
            selectInput("split_by_var", "Split by:",
                        choices = c("None" = "none", 
                                    "Replicate" = "replicate", 
                                    "Cell Line" = "cell_line", 
                                    "SPP1 Profile" = "spp1_profile",
                                    "Cisplatine" = "cisplatine",
                                    "Comment" = "comment"),
                        selected = "replicate"),
            checkboxInput("cluster_rows_viz", "Cluster rows", value = TRUE),
            checkboxInput("cluster_cols_viz", "Cluster columns", value = FALSE),
            
            actionButton("run_viz", "Visualize Genes", class = "btn-primary")
          ),
          
          box(
            title = "Gene Expression Heatmap",
            width = 9,
            solidHeader = TRUE,
            status = "info",
            plotOutput("gene_viz_heatmap", height = "700px")
          )
        )
      )
    )

# ============================================================================
# SHINY SERVER
# ============================================================================

server <- function(input, output, session) {
  
  # Initialize choices based on metadata
  observe({
    req(exists("metadata"))
    
    # Update filter choices
    if ("cell_line" %in% colnames(metadata)) {
      cell_lines <- unique(metadata$cell_line)
      updateCheckboxGroupInput(session, "filter_cell_line", choices = cell_lines, selected = cell_lines)
    }
    
    if ("spp1_profile" %in% colnames(metadata)) {
      spp1_profiles <- unique(metadata$spp1_profile)
      updateCheckboxGroupInput(session, "filter_spp1_profile", choices = spp1_profiles, selected = spp1_profiles)
    }
    
    if ("cisplatine" %in% colnames(metadata)) {
      cisplatine_vals <- unique(metadata$cisplatine)
      updateCheckboxGroupInput(session, "filter_cisplatine", choices = cisplatine_vals, selected = cisplatine_vals)
    }
    
    if ("comment" %in% colnames(metadata)) {
      comments <- unique(metadata$comment)
      updateCheckboxGroupInput(session, "filter_comment", choices = comments, selected = comments)
    }
    
    if ("replicate" %in% colnames(metadata)) {
      replicates <- unique(metadata$replicate)
      updateCheckboxGroupInput(session, "filter_replicate", choices = replicates, selected = replicates)
    }
  })
  
  # Gene Visualization
  viz_results <- eventReactive(input$run_viz, {
    genes <- trimws(strsplit(input$genes_input, "\n")[[1]])
    genes <- genes[genes != ""]
    
    req(length(genes) > 0)
    
    withProgress(message = 'Creating visualization...', value = 0, {
      tryCatch({
        visualize_genes_with_annotations(
          expr_data = expr_data_getmm,
          this_metadata = metadata,
          genes_of_interest = genes,
          cell_line = if(length(input$filter_cell_line) > 0) input$filter_cell_line else NULL,
          spp1_profile = if(length(input$filter_spp1_profile) > 0) input$filter_spp1_profile else NULL,
          cisplatine = if(length(input$filter_cisplatine) > 0) input$filter_cisplatine else NULL,
          comment = if(length(input$filter_comment) > 0) input$filter_comment else NULL,
          replicate = if(length(input$filter_replicate) > 0) input$filter_replicate else NULL,
          cluster_rows = input$cluster_rows_viz,
          cluster_columns = input$cluster_cols_viz,
          split_by = if(input$split_by_var == "none") NULL else input$split_by_var
        )
      }, error = function(e) {
        showNotification(paste("Error:", e$message), type = "error")
        NULL
      })
    })
  })
  
  output$gene_viz_heatmap <- renderPlot({
    req(viz_results())
    draw(viz_results()$heatmap)
  })
}

# ============================================================================
# RUN APP
# ============================================================================

shinyApp(ui = ui, server = server)