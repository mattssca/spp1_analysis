# Load required libraries
library(shiny)
library(shinydashboard)
library(tidyverse)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
library(circlize)
library(grid)
library(plotly)
library(DT)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# ============================================================================
# BUNDLED DATA
# ============================================================================

# Load data from the data folder
load("data/expr_data_getmm.Rdata")  
load("data/metadata.Rdata")
load("data/deg_results.Rdata")

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
    sort_by = NULL,
    auto_scale = TRUE,
    expr_min = NULL,
    expr_max = NULL,
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
  
  # Set expression color scale
  if (auto_scale) {
    expr_min <- min(mat, na.rm = TRUE)
    expr_max <- max(mat, na.rm = TRUE)
  } else {
    if (is.null(expr_min)) expr_min <- 0
    if (is.null(expr_max)) expr_max <- 5
  }
  expr_med <- (expr_min + expr_max) / 2
  
  expr_col_fun <- circlize::colorRamp2(
    c(expr_min, expr_med, expr_max), 
    c("#4DF76F", "black", "#F74D4D")
  )
  
  annotation_col <- meta_sub %>%
    dplyr::select(sample_id, cell_line, spp1_profile, cisplatine, comment, replicate) %>%
    column_to_rownames("sample_id")
  
  # Sort samples if sort_by is specified (before creating annotations)
  if (!cluster_columns && !is.null(sort_by)) {
    # Remove "none" values and keep only valid column names
    sort_vars <- sort_by[sort_by != "none" & sort_by %in% colnames(annotation_col)]
    
    if (length(sort_vars) > 0) {
      # Create a multi-level sort using do.call and order
      sort_order <- do.call(order, lapply(sort_vars, function(var) annotation_col[[var]]))
      annotation_col <- annotation_col[sort_order, , drop = FALSE]
      mat <- mat[, rownames(annotation_col), drop = FALSE]
    }
  }
  
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
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Gene Expression Heatmap", tabName = "heatmap", icon = icon("table")),
      menuItem("DEG Comparison", tabName = "deg_comparison", icon = icon("chart-line")),
      menuItem("Volcano Plot", tabName = "volcano_plot", icon = icon("volcano")),
      menuItem("GO Enrichment Analysis", tabName = "go_enrichment", icon = icon("sitemap")),
      menuItem("SPP1 Correlation", tabName = "spp1_correlation", icon = icon("project-diagram")),
      menuItem("Gene Expression Profile", tabName = "gene_profile", icon = icon("chart-area")),
      menuItem("PCA Analysis", tabName = "pca_analysis", icon = icon("project-diagram"))
    )
  ),
  
  dashboardBody(
    tabItems(
      # First tab - Gene Expression Heatmap
      tabItem(tabName = "heatmap",
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
            
            h4("Sort Samples (hierarchical):"),
            selectInput("sort_by_1", "Primary sort:",
                        choices = c("None" = "none", 
                                    "Replicate" = "replicate", 
                                    "Cell Line" = "cell_line", 
                                    "SPP1 Profile" = "spp1_profile",
                                    "Cisplatine" = "cisplatine",
                                    "Comment" = "comment"),
                        selected = "none"),
            selectInput("sort_by_2", "Secondary sort:",
                        choices = c("None" = "none", 
                                    "Replicate" = "replicate", 
                                    "Cell Line" = "cell_line", 
                                    "SPP1 Profile" = "spp1_profile",
                                    "Cisplatine" = "cisplatine",
                                    "Comment" = "comment"),
                        selected = "none"),
            selectInput("sort_by_3", "Tertiary sort:",
                        choices = c("None" = "none", 
                                    "Replicate" = "replicate", 
                                    "Cell Line" = "cell_line", 
                                    "SPP1 Profile" = "spp1_profile",
                                    "Cisplatine" = "cisplatine",
                                    "Comment" = "comment"),
                        selected = "none"),
            checkboxInput("cluster_rows_viz", "Cluster rows", value = TRUE),
            checkboxInput("cluster_cols_viz", "Cluster columns", value = FALSE),
            
            h4("Expression Scale:"),
            checkboxInput("auto_scale_expr", "Auto-scale to data", value = TRUE),
            conditionalPanel(
              condition = "!input.auto_scale_expr",
              sliderInput("expr_range", "Manual expression range:",
                         min = 0, max = 20, value = c(0, 5), step = 0.1)
            ),
            
            actionButton("run_viz", "Visualize Genes", class = "btn-primary")
          ),
          
          box(
            title = "Gene Expression Heatmap",
            width = 9,
            solidHeader = TRUE,
            status = "info",
            downloadButton("download_pdf", "Download PDF", class = "btn-success", 
                          style = "margin-bottom: 10px;"),
            InteractiveComplexHeatmapOutput("gene_viz_heatmap")
          )
        )
      ),
      
      # Second tab - DEG Comparison
      tabItem(tabName = "deg_comparison",
        fluidRow(
          box(
            title = "DEG Comparison Scatter Plot",
            width = 12,
            solidHeader = TRUE,
            status = "primary",
            
            fluidRow(
              column(3,
                h4("Select Comparisons"),
                selectInput("x_comparison", "X-axis (logFC):",
                            choices = NULL),
                selectInput("y_comparison", "Y-axis (logFC):",
                            choices = NULL),
                hr(),
                h4("Filtering Options"),
                numericInput("pvalue_threshold", "P-value threshold (genes with p < threshold):",
                            min = 0, max = 1, value = 0.05, step = 0.001),
                numericInput("adjp_threshold", "Adjusted P-value threshold (genes with adj.p < threshold):",
                            min = 0, max = 1, value = 0.05, step = 0.01),
                checkboxInput("show_sig_only", "Show only significant genes", value = FALSE),
                helpText("Lower thresholds = more stringent filtering (e.g., 0.001 for very significant genes)")
              ),
              column(9,
                downloadButton("download_deg_data", "Download DEG Comparison Data", 
                              class = "btn-success", style = "margin-bottom: 10px;"),
                plotlyOutput("deg_scatter", height = "600px"),
                hr(),
                h4("Top 10 Genes by Absolute logFC"),
                DTOutput("top_genes_table")
              )
            )
          )
        )
      ),
      
      # Volcano Plot tab
      tabItem(tabName = "volcano_plot",
        fluidRow(
          box(
            title = "Volcano Plot Visualization",
            width = 12,
            solidHeader = TRUE,
            status = "primary",
            
            fluidRow(
              column(3,
                h4("Select Comparison"),
                selectInput("volcano_comparison", "DEG Comparison:",
                            choices = NULL),
                hr(),
                h4("Significance Thresholds"),
                numericInput("volcano_pval_threshold", "-log10(adj.p) threshold:",
                            min = 0, max = 10, value = 1.3, step = 0.1),
                helpText("Default: 1.3 equals adj.p = 0.05"),
                numericInput("volcano_fc_threshold", "|logFC| threshold:",
                            min = 0, max = 10, value = 1, step = 0.1),
                hr(),
                h4("Labeling Options"),
                numericInput("volcano_top_n", "Label top N genes:",
                            min = 0, max = 50, value = 10, step = 1),
                helpText("Labels most significant genes"),
                hr(),
                h4("Color Scheme"),
                selectInput("volcano_colors", "Color by:",
                            choices = c("Significance" = "sig",
                                        "Expression change" = "fc"),
                            selected = "sig")
              ),
              column(9,
                downloadButton("download_volcano_data", "Download Volcano Data", 
                              class = "btn-success", style = "margin-bottom: 10px;"),
                plotlyOutput("volcano_plot", height = "600px"),
                hr(),
                h4("Significant Genes Summary"),
                DTOutput("volcano_summary_table")
              )
            )
          )
        )
      ),
      
      # Third tab - GO Enrichment Analysis
      tabItem(tabName = "go_enrichment",
        fluidRow(
          box(
            title = "Pathway Enrichment Parameters",
            width = 3,
            solidHeader = TRUE,
            status = "primary",
            
            h4("Select Comparison"),
            selectInput("go_comparison", "DEG Comparison:",
                        choices = NULL),
            hr(),
            h4("Database Selection"),
            radioButtons("pathway_database", "Enrichment Database:",
                        choices = c("Gene Ontology (GO)" = "GO",
                                    "KEGG Pathways" = "KEGG"),
                        selected = "GO"),
            helpText(HTML("<b>GO:</b> Gene Ontology terms (BP, MF, CC)<br><b>KEGG:</b> KEGG pathway database")),
            hr(),
            h4("Analysis Method"),
            radioButtons("go_method", "Enrichment Method:",
                        choices = c("Over-Representation Analysis (ORA)" = "ORA",
                                    "Gene Set Enrichment Analysis (GSEA)" = "GSEA"),
                        selected = "ORA"),
            helpText(HTML("<b>ORA:</b> Tests if terms are over-represented in genes passing thresholds. 
                          <br><b>GSEA:</b> Tests if term genes are enriched at top/bottom of ranked gene list. 
                          Uses all genes, no cutoffs needed.")),
            hr(),
            conditionalPanel(
              condition = "input.go_method == 'ORA'",
              h4("Gene Selection Thresholds"),
              numericInput("go_logfc_threshold", "Log Fold Change threshold (|logFC| > threshold):",
                          min = 0, max = 10, value = 1, step = 0.1),
              numericInput("go_adjp_threshold", "Adjusted P-value threshold (adj.p < threshold):",
                          min = 0, max = 1, value = 0.05, step = 0.01),
              helpText("Genes exceeding these thresholds will be used for enrichment"),
              hr()
            ),
            conditionalPanel(
              condition = "input.go_method == 'GSEA'",
              h4("Gene Ranking"),
              radioButtons("gsea_rank_by", "Rank genes by:",
                          choices = c("Log Fold Change" = "logFC",
                                      "t-statistic" = "t",
                                      "-log10(p) × sign(logFC)" = "signed_pval"),
                          selected = "logFC"),
              helpText("Determines how genes are ranked for GSEA"),
              hr()
            ),
            conditionalPanel(
              condition = "input.pathway_database == 'GO'",
              h4("GO Ontology"),
              selectInput("go_ontology", "Ontology:",
                          choices = c("Biological Process" = "BP",
                                      "Molecular Function" = "MF",
                                      "Cellular Component" = "CC"),
                          selected = "BP")
            ),
            h4("Enrichment Options"),
            numericInput("go_pval_cutoff", "P-value cutoff:",
                        min = 0, max = 1, value = 0.05, step = 0.01),
            numericInput("go_qval_cutoff", "Q-value cutoff:",
                        min = 0, max = 1, value = 0.2, step = 0.05),
            helpText("P-value and Q-value cutoffs for term significance"),
            hr(),
            actionButton("run_go_enrichment", "Run Enrichment", class = "btn-primary"),
            hr(),
            h4("Summary"),
            textOutput("go_gene_count"),
            textOutput("go_term_count")
          ),
          
          box(
            title = "Pathway Enrichment Results",
            width = 9,
            solidHeader = TRUE,
            status = "info",
            
            tabsetPanel(
              tabPanel("Bubble Plot",
                plotlyOutput("go_bubble_plot", height = "600px")
              ),
              tabPanel("Summary Table",
                DTOutput("go_summary_table")
              ),
              tabPanel("Full Results",
                downloadButton("download_go_results", "Download Full Results", 
                              class = "btn-success", style = "margin-bottom: 10px;"),
                DTOutput("go_full_table")
              )
            )
          )
        )
      ),
      
      # Fourth tab - SPP1 Correlation Volcano Plot
      tabItem(tabName = "spp1_correlation",
        fluidRow(
          box(
            title = "Sample Selection",
            width = 3,
            solidHeader = TRUE,
            status = "primary",
            
            h4("Filters:"),
            checkboxGroupInput("filter_cell_line_spp1", "Cell Line:",
                               choices = NULL),
            checkboxGroupInput("filter_spp1_profile_spp1", "SPP1 Profile:",
                               choices = NULL),
            checkboxGroupInput("filter_cisplatine_spp1", "Cisplatine:",
                               choices = NULL),
            checkboxGroupInput("filter_comment_spp1", "Comment:",
                               choices = NULL),
            checkboxGroupInput("filter_replicate_spp1", "Replicate:",
                               choices = NULL),
            hr(),
            h4("Correlation Thresholds:"),
            numericInput("cor_threshold", "Correlation threshold (|r| > threshold):",
                        min = 0, max = 1, value = 0.5, step = 0.05),
            numericInput("pval_threshold_cor", "P-value threshold (p < threshold):",
                        min = 0, max = 1, value = 0.05, step = 0.001),
            helpText("Higher correlation threshold = more stringent filtering"),
            hr(),
            h4("Query Specific Gene"),
            selectizeInput("query_gene_spp1", "Search for a gene:",
                          choices = NULL,
                          options = list(
                            placeholder = 'Type to search...',
                            maxOptions = 50
                          )),
            conditionalPanel(
              condition = "input.query_gene_spp1 != null && input.query_gene_spp1 != ''",
              actionButton("clear_gene_query", "Clear Selection", class = "btn-warning btn-sm")
            ),
            hr(),
            actionButton("run_spp1_cor", "Generate Correlation Plot", class = "btn-primary")
          ),
          
          box(
            title = "SPP1 Correlation Volcano Plot",
            width = 9,
            solidHeader = TRUE,
            status = "info",
            downloadButton("download_spp1_cor", "Download Correlation Data", 
                          class = "btn-success", style = "margin-bottom: 10px;"),
            plotlyOutput("spp1_volcano", height = "500px"),
            conditionalPanel(
              condition = "input.query_gene_spp1 != null && input.query_gene_spp1 != ''",
              hr(),
              h4("Correlation Details for Selected Gene"),
              fluidRow(
                column(6,
                  DTOutput("gene_correlation_details")
                ),
                column(6,
                  plotlyOutput("gene_spp1_scatter", height = "300px")
                )
              )
            ),
            hr(),
            fluidRow(
              column(6,
                h4("Top 10 Positively Correlated Genes"),
                DTOutput("top_positive_genes")
              ),
              column(6,
                h4("Top 10 Negatively Correlated Genes"),
                DTOutput("top_negative_genes")
              )
            )
          )
        )
      ),
      
      # Fifth tab - Gene Expression Profile
      tabItem(tabName = "gene_profile",
        fluidRow(
          box(
            title = "Gene Selection",
            width = 3,
            solidHeader = TRUE,
            status = "primary",
            
            h4("Select Gene"),
            selectizeInput("selected_gene", "Gene:",
                          choices = NULL,
                          options = list(
                            placeholder = 'Type to search for a gene...',
                            maxOptions = 50
                          )),
            hr(),
            h4("Grouping Variable"),
            selectInput("group_by_var", "Group samples by:",
                        choices = c("Cell Line" = "cell_line",
                                    "SPP1 Profile" = "spp1_profile",
                                    "Cisplatine" = "cisplatine",
                                    "Comment" = "comment",
                                    "Replicate" = "replicate"),
                        selected = "spp1_profile"),
            selectInput("color_by_var", "Color by:",
                        choices = c("Cell Line" = "cell_line",
                                    "SPP1 Profile" = "spp1_profile",
                                    "Cisplatine" = "cisplatine",
                                    "Comment" = "comment",
                                    "Replicate" = "replicate",
                                    "Same as grouping" = "same"),
                        selected = "cell_line"),
            hr(),
            h4("Filters:"),
            checkboxGroupInput("filter_cell_line_profile", "Cell Line:",
                               choices = NULL),
            checkboxGroupInput("filter_spp1_profile_profile", "SPP1 Profile:",
                               choices = NULL),
            checkboxGroupInput("filter_cisplatine_profile", "Cisplatine:",
                               choices = NULL),
            checkboxGroupInput("filter_comment_profile", "Comment:",
                               choices = NULL),
            checkboxGroupInput("filter_replicate_profile", "Replicate:",
                               choices = NULL),
            hr(),
            checkboxInput("show_points", "Show individual points", value = TRUE),
            actionButton("run_gene_profile", "Generate Profile", class = "btn-primary")
          ),
          
          box(
            title = "Gene Expression Profile",
            width = 9,
            solidHeader = TRUE,
            status = "info",
            downloadButton("download_gene_profile", "Download Profile Data", 
                          class = "btn-success", style = "margin-bottom: 10px;"),
            plotlyOutput("gene_profile_plot", height = "600px"),
            hr(),
            h4("Summary Statistics"),
            DTOutput("gene_profile_stats")
          )
        )
      ),
      
      # PCA Analysis tab  
      tabItem(tabName = "pca_analysis",
        fluidRow(
          box(
            title = "PCA Parameters",
            width = 3,
            solidHeader = TRUE,
            status = "primary",
            
            h4("Sample Filters:"),
            checkboxGroupInput("filter_cell_line_pca", "Cell Line:",
                               choices = NULL),
            checkboxGroupInput("filter_spp1_profile_pca", "SPP1 Profile:",
                               choices = NULL),
            checkboxGroupInput("filter_cisplatine_pca", "Cisplatine:",
                               choices = NULL),
            checkboxGroupInput("filter_comment_pca", "Comment:",
                               choices = NULL),
            checkboxGroupInput("filter_replicate_pca", "Replicate:",
                               choices = NULL),
            hr(),
            h4("PCA Options"),
            selectInput("pca_color_by", "Color samples by:",
                        choices = c("Cell Line" = "cell_line",
                                    "SPP1 Profile" = "spp1_profile",
                                    "Cisplatine" = "cisplatine",
                                    "Comment" = "comment",
                                    "Replicate" = "replicate"),
                        selected = "spp1_profile"),
            selectInput("pca_shape_by", "Shape samples by:",
                        choices = c("None" = "none",
                                    "Cell Line" = "cell_line",
                                    "SPP1 Profile" = "spp1_profile",
                                    "Cisplatine" = "cisplatine",
                                    "Comment" = "comment",
                                    "Replicate" = "replicate"),
                        selected = "cell_line"),
            numericInput("pca_top_genes", "Number of most variable genes:",
                        min = 100, max = 10000, value = 500, step = 100),
            helpText("Uses top N genes by variance for PCA"),
            hr(),
            checkboxInput("pca_scale", "Scale genes", value = TRUE),
            actionButton("run_pca", "Run PCA", class = "btn-primary")
          ),
          
          box(
            title = "PCA Visualization",
            width = 9,
            solidHeader = TRUE,
            status = "info",
            
            tabsetPanel(
              tabPanel("PCA Plot",
                downloadButton("download_pca_plot", "Download PCA Data", 
                              class = "btn-success", style = "margin: 10px;"),
                plotlyOutput("pca_plot", height = "600px")
              ),
              tabPanel("Scree Plot",
                plotlyOutput("scree_plot", height = "500px"),
                hr(),
                h4("Variance Explained"),
                DTOutput("variance_table")
              ),
              tabPanel("Loadings",
                h4("Top Gene Loadings for Selected PC"),
                fluidRow(
                  column(3,
                    selectInput("pc_for_loadings", "Select PC:",
                                choices = c("PC1" = 1, "PC2" = 2, "PC3" = 3, "PC4" = 4, "PC5" = 5),
                                selected = 1)
                  )
                ),
                DTOutput("loadings_table")
              )
            )
          )
        )
      )
    )
  )
)

# ============================================================================
# SHINY SERVER
# ============================================================================

server <- function(input, output, session) {
  
  # ============================================================================
  # CACHE FOR PERFORMANCE OPTIMIZATION
  # ============================================================================
  
  # Cache for gene ID conversions to speed up repeated GO analyses
  gene_id_cache <- reactiveValues(
    symbol_to_entrez = NULL
  )
  
  # Pre-build gene ID mapping on first use (lazy loading)
  observe({
    # This will be triggered when data is available
    req(exists("expr_data_getmm"))
    
    # Only build once
    if (is.null(gene_id_cache$symbol_to_entrez)) {
      # Build mapping for all genes in the dataset to speed up future lookups
      all_genes <- rownames(expr_data_getmm)
      
      tryCatch({
        gene_id_cache$symbol_to_entrez <- bitr(
          all_genes,
          fromType = "SYMBOL",
          toType = "ENTREZID",
          OrgDb = org.Hs.eg.db
        )
      }, error = function(e) {
        # If it fails, set to empty data frame to avoid repeated attempts
        gene_id_cache$symbol_to_entrez <- data.frame(SYMBOL = character(), ENTREZID = character())
      })
    }
  })
  
  # Initialize choices based on metadata
  observe({
    req(exists("metadata"))
    
    # Update filter choices for all tabs
    if ("cell_line" %in% colnames(metadata)) {
      cell_lines <- unique(metadata$cell_line)
      updateCheckboxGroupInput(session, "filter_cell_line", choices = cell_lines, selected = cell_lines)
      updateCheckboxGroupInput(session, "filter_cell_line_spp1", choices = cell_lines, selected = cell_lines)
      updateCheckboxGroupInput(session, "filter_cell_line_profile", choices = cell_lines, selected = cell_lines)
    }
    
    if ("spp1_profile" %in% colnames(metadata)) {
      spp1_profiles <- unique(metadata$spp1_profile)
      updateCheckboxGroupInput(session, "filter_spp1_profile", choices = spp1_profiles, selected = spp1_profiles)
      updateCheckboxGroupInput(session, "filter_spp1_profile_spp1", choices = spp1_profiles, selected = spp1_profiles)
      updateCheckboxGroupInput(session, "filter_spp1_profile_profile", choices = spp1_profiles, selected = spp1_profiles)
    }
    
    if ("cisplatine" %in% colnames(metadata)) {
      cisplatine_vals <- unique(metadata$cisplatine)
      updateCheckboxGroupInput(session, "filter_cisplatine", choices = cisplatine_vals, selected = cisplatine_vals)
      updateCheckboxGroupInput(session, "filter_cisplatine_spp1", choices = cisplatine_vals, selected = cisplatine_vals)
      updateCheckboxGroupInput(session, "filter_cisplatine_profile", choices = cisplatine_vals, selected = cisplatine_vals)
    }
    
    if ("comment" %in% colnames(metadata)) {
      comments <- unique(metadata$comment)
      updateCheckboxGroupInput(session, "filter_comment", choices = comments, selected = comments)
      updateCheckboxGroupInput(session, "filter_comment_spp1", choices = comments, selected = comments)
      updateCheckboxGroupInput(session, "filter_comment_profile", choices = comments, selected = comments)
    }
    
    if ("replicate" %in% colnames(metadata)) {
      replicates <- unique(metadata$replicate)
      updateCheckboxGroupInput(session, "filter_replicate", choices = replicates, selected = replicates)
      updateCheckboxGroupInput(session, "filter_replicate_spp1", choices = replicates, selected = replicates)
      updateCheckboxGroupInput(session, "filter_replicate_profile", choices = replicates, selected = replicates)
    }
  })
  
  # Initialize gene choices for gene profile tab
  observe({
    req(exists("expr_data_getmm"))
    gene_list <- sort(rownames(expr_data_getmm))
    updateSelectizeInput(session, "selected_gene", choices = gene_list, server = TRUE)
    updateSelectizeInput(session, "query_gene_spp1", choices = gene_list, server = TRUE)
  })
  
  # Clear gene query when button is clicked
  observeEvent(input$clear_gene_query, {
    updateSelectizeInput(session, "query_gene_spp1", selected = "")
  })
  
  # Gene Visualization
  viz_results <- eventReactive(input$run_viz, {
    genes <- trimws(strsplit(input$genes_input, "\n")[[1]])
    genes <- genes[genes != ""]
    
    req(length(genes) > 0)
    
    withProgress(message = 'Creating visualization...', value = 0, {
      tryCatch({
        # Get expression range values, use defaults if auto_scale is TRUE or slider not available
        expr_min_val <- if(!is.null(input$expr_range)) input$expr_range[1] else 0
        expr_max_val <- if(!is.null(input$expr_range)) input$expr_range[2] else 5
        
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
          split_by = if(input$split_by_var == "none") NULL else input$split_by_var,
          sort_by = c(input$sort_by_1, input$sort_by_2, input$sort_by_3),
          auto_scale = input$auto_scale_expr,
          expr_min = expr_min_val,
          expr_max = expr_max_val
        )
      }, error = function(e) {
        showNotification(paste("Error:", e$message), type = "error")
        NULL
      })
    })
  })
  
  # Gene Visualization - Interactive Heatmap
  observeEvent(viz_results(), {
    req(viz_results())
    makeInteractiveComplexHeatmap(input, output, session, viz_results()$heatmap, "gene_viz_heatmap")
  })
  
  # Download PDF handler
  output$download_pdf <- downloadHandler(
    filename = function() {
      paste0("gene_expression_heatmap_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      req(viz_results())
      pdf(file, width = 14, height = 10)
      draw(viz_results()$heatmap)
      dev.off()
    }
  )
  
  # ============================================================================
  # DEG COMPARISON TAB
  # ============================================================================
  
  # Initialize comparison choices based on deg_results
  observe({
    req(exists("deg_results"))
    
    comparison_names <- names(deg_results)
    updateSelectInput(session, "x_comparison", 
                      choices = comparison_names,
                      selected = comparison_names[1])
    updateSelectInput(session, "y_comparison", 
                      choices = comparison_names,
                      selected = comparison_names[min(2, length(comparison_names))])
  })
  
  # Prepare merged data for scatter plot
  scatter_data <- reactive({
    req(input$x_comparison, input$y_comparison)
    
    x_data <- deg_results[[input$x_comparison]] %>%
      dplyr::select(gene, logFC_x = logFC, P.Value_x = P.Value, adj.P.Val_x = adj.P.Val)
    
    y_data <- deg_results[[input$y_comparison]] %>%
      dplyr::select(gene, logFC_y = logFC, P.Value_y = P.Value, adj.P.Val_y = adj.P.Val)
    
    # Merge the two datasets
    merged <- inner_join(x_data, y_data, by = "gene") %>%
      mutate(
        sig_category = case_when(
          P.Value_x < input$pvalue_threshold & P.Value_y < input$pvalue_threshold & 
            adj.P.Val_x < input$adjp_threshold & adj.P.Val_y < input$adjp_threshold ~ "Sig in both (p & adj.p)",
          P.Value_x < input$pvalue_threshold & P.Value_y < input$pvalue_threshold ~ "Sig in both (p-value only)",
          adj.P.Val_x < input$adjp_threshold & adj.P.Val_y < input$adjp_threshold ~ "Sig in both (adj.p only)",
          P.Value_x < input$pvalue_threshold | P.Value_y < input$pvalue_threshold ~ "Sig in one (p-value)",
          adj.P.Val_x < input$adjp_threshold | adj.P.Val_y < input$adjp_threshold ~ "Sig in one (adj.p)",
          TRUE ~ "Not significant"
        ),
        min_pvalue = pmin(P.Value_x, P.Value_y),
        min_adjp = pmin(adj.P.Val_x, adj.P.Val_y)
      )
    
    if (input$show_sig_only) {
      merged <- merged %>%
        filter(sig_category != "Not significant")
    }
    
    merged
  })
  
  # Create scatter plot
  output$deg_scatter <- renderPlotly({
    req(scatter_data())
    
    data <- scatter_data()
    
    # Define colors for significance categories
    color_mapping <- c(
      "Sig in both (p & adj.p)" = "#e41a1c",
      "Sig in both (p-value only)" = "#ff7f00",
      "Sig in both (adj.p only)" = "#984ea3",
      "Sig in one (p-value)" = "#4daf4a",
      "Sig in one (adj.p)" = "#377eb8",
      "Not significant" = "#999999"
    )
    
    p <- plot_ly(data, 
                 x = ~logFC_x, 
                 y = ~logFC_y,
                 type = 'scatter',
                 mode = 'markers',
                 color = ~sig_category,
                 colors = color_mapping,
                 marker = list(
                   size = 6,
                   opacity = 0.6,
                   line = list(width = 0)
                 ),
                 text = ~paste0(
                   "Gene: ", gene,
                   "<br>", input$x_comparison, " logFC: ", round(logFC_x, 3),
                   "<br>", input$y_comparison, " logFC: ", round(logFC_y, 3),
                   "<br>P-value (x): ", format(P.Value_x, scientific = TRUE, digits = 3),
                   "<br>P-value (y): ", format(P.Value_y, scientific = TRUE, digits = 3),
                   "<br>Adj.P (x): ", format(adj.P.Val_x, scientific = TRUE, digits = 3),
                   "<br>Adj.P (y): ", format(adj.P.Val_y, scientific = TRUE, digits = 3)
                 ),
                 hoverinfo = 'text'
    ) %>%
      layout(
        title = "DEG Comparison",
        xaxis = list(title = paste(input$x_comparison, "(logFC)")),
        yaxis = list(title = paste(input$y_comparison, "(logFC)")),
        hovermode = 'closest',
        showlegend = TRUE,
        legend = list(title = list(text = "Significance"))
      ) %>%
      add_segments(x = 0, xend = 0, y = min(data$logFC_y), yend = max(data$logFC_y),
                   line = list(color = "black", dash = "dash", width = 1),
                   showlegend = FALSE, inherit = FALSE) %>%
      add_segments(x = min(data$logFC_x), xend = max(data$logFC_x), y = 0, yend = 0,
                   line = list(color = "black", dash = "dash", width = 1),
                   showlegend = FALSE, inherit = FALSE)
    
    p
  })
  
  # Top genes table
  output$top_genes_table <- renderDT({
    req(input$x_comparison)
    
    # Get data from the x-axis comparison
    top_genes <- deg_results[[input$x_comparison]] %>%
      arrange(desc(abs(logFC))) %>%
      head(10) %>%
      mutate(
        logFC = round(logFC, 3),
        AveExpr = round(AveExpr, 3),
        t = round(t, 3),
        P.Value = format(P.Value, scientific = TRUE, digits = 3),
        adj.P.Val = format(adj.P.Val, scientific = TRUE, digits = 3),
        B = round(B, 3)
      )
    
    datatable(
      top_genes,
      options = list(
        pageLength = 10,
        dom = 't',
        ordering = FALSE
      ),
      rownames = FALSE,
      caption = paste("Top 10 genes by absolute logFC in", input$x_comparison)
    )
  })
  
  # Download handler for DEG comparison data
  output$download_deg_data <- downloadHandler(
    filename = function() {
      paste0("DEG_comparison_", input$x_comparison, "_vs_", input$y_comparison, "_", 
             Sys.Date(), ".csv")
    },
    content = function(file) {
      req(scatter_data())
      write.csv(scatter_data(), file, row.names = FALSE)
    }
  )
  
  # ============================================================================
  # VOLCANO PLOT TAB
  # ============================================================================
  
  # Initialize comparison choices for volcano plot
  observe({
    req(exists("deg_results"))
    
    comparison_names <- names(deg_results)
    updateSelectInput(session, "volcano_comparison", 
                      choices = comparison_names,
                      selected = comparison_names[1])
  })
  
  # Prepare volcano plot data
  volcano_data <- reactive({
    req(input$volcano_comparison)
    
    deg_data <- deg_results[[input$volcano_comparison]]
    
    # Calculate -log10(adj.P.Val) and classify genes
    volcano_df <- deg_data %>%
      mutate(
        neg_log10_pval = -log10(adj.P.Val),
        abs_logFC = abs(logFC),
        # Classify genes by significance
        significance = case_when(
          abs(logFC) >= input$volcano_fc_threshold & neg_log10_pval >= input$volcano_pval_threshold ~ "Significant",
          abs(logFC) >= input$volcano_fc_threshold ~ "FC only",
          neg_log10_pval >= input$volcano_pval_threshold ~ "P-val only",
          TRUE ~ "Not significant"
        ),
        # Direction for coloring
        direction = case_when(
          logFC > input$volcano_fc_threshold & neg_log10_pval >= input$volcano_pval_threshold ~ "Up-regulated",
          logFC < -input$volcano_fc_threshold & neg_log10_pval >= input$volcano_pval_threshold ~ "Down-regulated",
          TRUE ~ "Not significant"
        )
      )
    
    volcano_df
  })
  
  # Render volcano plot
  output$volcano_plot <- renderPlotly({
    req(volcano_data())
    
    df <- volcano_data()
    
    # Select color scheme
    if (input$volcano_colors == "sig") {
      color_col <- "significance"
      colors <- c("Significant" = "#E74C3C", "FC only" = "#F39C12", 
                  "P-val only" = "#3498DB", "Not significant" = "grey70")
    } else {
      color_col <- "direction"
      colors <- c("Up-regulated" = "#E74C3C", "Down-regulated" = "#3498DB", 
                  "Not significant" = "grey70")
    }
    
    # Get top N genes to label
    top_genes <- df %>%
      filter(significance == "Significant") %>%
      arrange(desc(abs_logFC * neg_log10_pval)) %>%
      head(input$volcano_top_n)
    
    p <- plot_ly(df, x = ~logFC, y = ~neg_log10_pval, 
                 type = 'scatter', mode = 'markers',
                 color = ~get(color_col), colors = colors,
                 marker = list(size = 6, opacity = 0.6),
                 text = ~paste("Gene:", gene,
                              "<br>logFC:", round(logFC, 2),
                              "<br>-log10(adj.p):", round(neg_log10_pval, 2),
                              "<br>adj.P.Val:", format(adj.P.Val, digits = 3)),
                 hoverinfo = 'text') %>%
      add_trace(data = top_genes, 
                mode = 'text',
                text = ~gene,
                textposition = "top center",
                textfont = list(size = 10, color = "black"),
                showlegend = FALSE,
                hoverinfo = 'skip') %>%
      layout(
        title = paste("Volcano Plot:", input$volcano_comparison),
        xaxis = list(title = "Log2 Fold Change", zeroline = TRUE),
        yaxis = list(title = "-log10(Adjusted P-value)"),
        hovermode = 'closest',
        shapes = list(
          # Vertical lines for FC threshold
          list(type = "line", x0 = input$volcano_fc_threshold, x1 = input$volcano_fc_threshold,
               y0 = 0, y1 = 1, yref = "paper",
               line = list(color = "grey", dash = "dash", width = 1)),
          list(type = "line", x0 = -input$volcano_fc_threshold, x1 = -input$volcano_fc_threshold,
               y0 = 0, y1 = 1, yref = "paper",
               line = list(color = "grey", dash = "dash", width = 1)),
          # Horizontal line for p-value threshold
          list(type = "line", x0 = 0, x1 = 1, xref = "paper",
               y0 = input$volcano_pval_threshold, y1 = input$volcano_pval_threshold,
               line = list(color = "grey", dash = "dash", width = 1))
        )
      )
    
    p
  })
  
  # Summary table for volcano plot
  output$volcano_summary_table <- renderDT({
    req(volcano_data())
    
    df <- volcano_data()
    
    summary_stats <- df %>%
      summarise(
        `Total genes` = n(),
        `Significant (both)` = sum(significance == "Significant"),
        `Up-regulated` = sum(direction == "Up-regulated"),
        `Down-regulated` = sum(direction == "Down-regulated"),
        `FC only` = sum(significance == "FC only"),
        `P-value only` = sum(significance == "P-val only")
      ) %>%
      tidyr::pivot_longer(everything(), names_to = "Category", values_to = "Count")
    
    datatable(summary_stats, 
              options = list(dom = 't', ordering = FALSE),
              rownames = FALSE)
  })
  
  # Download handler for volcano data
  output$download_volcano_data <- downloadHandler(
    filename = function() {
      paste0("volcano_plot_", input$volcano_comparison, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(volcano_data())
      write.csv(volcano_data(), file, row.names = FALSE)
    }
  )
  
  # ============================================================================
  # GO ENRICHMENT ANALYSIS TAB
  # ============================================================================
  
  # Initialize comparison choices for GO analysis
  observe({
    req(exists("deg_results"))
    
    comparison_names <- names(deg_results)
    updateSelectInput(session, "go_comparison", 
                      choices = comparison_names,
                      selected = comparison_names[1])
  })
  
  # Perform GO/KEGG enrichment when button is clicked
  go_enrichment_results <- eventReactive(input$run_go_enrichment, {
    withProgress(message = 'Running pathway enrichment...', value = 0, {
      tryCatch({
        req(input$go_comparison, input$go_method, input$pathway_database)
        
        # Get DEG data for the selected comparison
        deg_data <- deg_results[[input$go_comparison]]
        
        if (input$go_method == "ORA") {
          # ============ Over-Representation Analysis ============
          incProgress(0.1, detail = "Filtering genes...")
          
          # Filter genes based on thresholds
          selected_genes <- deg_data %>%
            filter(
              abs(logFC) > input$go_logfc_threshold,
              adj.P.Val < input$go_adjp_threshold
            )
          
          if (nrow(selected_genes) == 0) {
            return(list(
              error = TRUE,
              message = "No genes meet the specified thresholds. Please adjust the thresholds.",
              method = "ORA",
              database = input$pathway_database
            ))
          }
          
          # Get gene symbols
          gene_list <- selected_genes$gene
          
          incProgress(0.2, detail = "Converting gene IDs...")
          
          # Use cached gene ID mapping for faster conversion
          if (!is.null(gene_id_cache$symbol_to_entrez) && nrow(gene_id_cache$symbol_to_entrez) > 0) {
            # Use pre-built cache
            gene_entrez <- gene_id_cache$symbol_to_entrez %>%
              filter(SYMBOL %in% gene_list)
          } else {
            # Fallback: convert on-the-fly if cache not available
            gene_entrez <- bitr(gene_list, 
                               fromType = "SYMBOL",
                               toType = "ENTREZID",
                               OrgDb = org.Hs.eg.db)
          }
          
          if (nrow(gene_entrez) == 0) {
            return(list(
              error = TRUE,
              message = "No genes could be mapped to Entrez IDs. Check gene symbols.",
              method = "ORA",
              database = input$pathway_database
            ))
          }
          
          incProgress(0.4, detail = "Running ORA enrichment...")
          
          # Perform enrichment based on database selection
          if (input$pathway_database == "GO") {
            # GO enrichment
            enrich_result <- enrichGO(gene = gene_entrez$ENTREZID,
                           OrgDb = org.Hs.eg.db,
                           ont = input$go_ontology,
                           pAdjustMethod = "BH",
                           pvalueCutoff = input$go_pval_cutoff,
                           qvalueCutoff = input$go_qval_cutoff,
                           readable = TRUE,
                           pool = FALSE)
          } else {
            # KEGG pathway enrichment
            enrich_result <- enrichKEGG(gene = gene_entrez$ENTREZID,
                           organism = "hsa",  # Homo sapiens
                           pAdjustMethod = "BH",
                           pvalueCutoff = input$go_pval_cutoff,
                           qvalueCutoff = input$go_qval_cutoff)
          }
          
          incProgress(0.3, detail = "Processing results...")
          
          if (is.null(enrich_result) || nrow(enrich_result@result) == 0) {
            return(list(
              error = TRUE,
              message = paste("No significant", input$pathway_database, "terms found. Try relaxing the cutoff values."),
              method = "ORA",
              database = input$pathway_database
            ))
          }
          
          # Return results
          list(
            error = FALSE,
            results = enrich_result,
            method = "ORA",
            database = input$pathway_database,
            n_genes = length(gene_list),
            n_mapped = nrow(gene_entrez),
            n_terms = nrow(enrich_result@result)
          )
          
        } else {
          # ============ Gene Set Enrichment Analysis (GSEA) ============
          incProgress(0.1, detail = "Preparing gene list...")
          
          # Create ranked gene list based on selected metric
          if (input$gsea_rank_by == "logFC") {
            ranked_values <- deg_data$logFC
          } else if (input$gsea_rank_by == "t") {
            ranked_values <- deg_data$t
          } else {  # signed_pval
            ranked_values <- -log10(deg_data$P.Value) * sign(deg_data$logFC)
          }
          
          # Create named vector and remove any NAs or Inf values
          gene_list_ranked <- setNames(ranked_values, deg_data$gene)
          gene_list_ranked <- gene_list_ranked[!is.na(gene_list_ranked) & 
                                                !is.infinite(gene_list_ranked)]
          
          # Sort in decreasing order
          gene_list_ranked <- sort(gene_list_ranked, decreasing = TRUE)
          
          incProgress(0.2, detail = "Converting gene IDs...")
          
          # Convert gene symbols to Entrez IDs
          if (!is.null(gene_id_cache$symbol_to_entrez) && nrow(gene_id_cache$symbol_to_entrez) > 0) {
            # Use pre-built cache
            gene_mapping <- gene_id_cache$symbol_to_entrez
          } else {
            # Fallback: convert on-the-fly
            gene_mapping <- bitr(names(gene_list_ranked), 
                                fromType = "SYMBOL",
                                toType = "ENTREZID",
                                OrgDb = org.Hs.eg.db)
          }
          
          if (nrow(gene_mapping) == 0) {
            return(list(
              error = TRUE,
              message = "No genes could be mapped to Entrez IDs. Check gene symbols.",
              method = "GSEA",
              database = input$pathway_database
            ))
          }
          
          # Create mapping preserving order
          # Match genes to get their ranking values
          gene_df <- data.frame(
            SYMBOL = names(gene_list_ranked),
            rank_value = as.numeric(gene_list_ranked),
            stringsAsFactors = FALSE
          )
          
          # Merge with Entrez IDs
          gene_df <- merge(gene_df, gene_mapping, by = "SYMBOL", all.x = FALSE)
          
          # Remove duplicates (keep first/highest ranked)
          gene_df <- gene_df[!duplicated(gene_df$ENTREZID), ]
          
          # Sort by rank_value in decreasing order
          gene_df <- gene_df[order(gene_df$rank_value, decreasing = TRUE), ]
          
          # Create final gene list
          gene_list_entrez <- setNames(gene_df$rank_value, gene_df$ENTREZID)
          
          # Final check: ensure it's properly sorted and has no issues
          if (!all(diff(gene_list_entrez) <= 0)) {
            # Re-sort if somehow not sorted
            gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)
          }
          
          incProgress(0.4, detail = "Running GSEA...")
          
          # Perform GSEA based on database selection
          if (input$pathway_database == "GO") {
            # GO GSEA
            gsea_result <- gseGO(
              geneList = gene_list_entrez,
              OrgDb = org.Hs.eg.db,
              ont = input$go_ontology,
              pvalueCutoff = input$go_pval_cutoff,
              pAdjustMethod = "BH",
              nPermSimple = 1000,
              verbose = FALSE
            )
          } else {
            # KEGG GSEA
            gsea_result <- gseKEGG(
              geneList = gene_list_entrez,
              organism = "hsa",
              pvalueCutoff = input$go_pval_cutoff,
              pAdjustMethod = "BH",
              nPermSimple = 1000,
              verbose = FALSE
            )
          }
          
          incProgress(0.3, detail = "Processing results...")
          
          if (is.null(gsea_result) || nrow(gsea_result@result) == 0) {
            return(list(
              error = TRUE,
              message = paste("No significant", input$pathway_database, "terms found. Try relaxing the p-value cutoff."),
              method = "GSEA",
              database = input$pathway_database
            ))
          }
          
          # Return results
          list(
            error = FALSE,
            results = gsea_result,
            method = "GSEA",
            database = input$pathway_database,
            n_genes = length(gene_list_ranked),
            n_mapped = length(gene_list_entrez),
            n_terms = nrow(gsea_result@result)
          )
        }
        
      }, error = function(e) {
        return(list(
          error = TRUE,
          message = paste("Error:", e$message),
          method = input$go_method,
          database = input$pathway_database
        ))
      })
    })
  })
  
  # Display gene count summary
  output$go_gene_count <- renderText({
    req(go_enrichment_results())
    results <- go_enrichment_results()
    
    if (results$error) {
      return("")
    }
    
    if (results$method == "ORA") {
      paste0("Genes selected: ", results$n_genes, " (", results$n_mapped, " mapped to Entrez ID)")
    } else {
      paste0("Total genes ranked: ", results$n_genes, " (", results$n_mapped, " mapped to Entrez ID)")
    }
  })
  
  # Display term count summary
  output$go_term_count <- renderText({
    req(go_enrichment_results())
    results <- go_enrichment_results()
    
    if (results$error) {
      return("")
    }
    
    term_type <- if (results$database == "GO") "GO terms" else "KEGG pathways"
    paste0("Significant ", term_type, ": ", results$n_terms)
  })
  
  # Create bubble plot
  output$go_bubble_plot <- renderPlotly({
    req(go_enrichment_results())
    results <- go_enrichment_results()
    
    if (results$error) {
      # Show error message
      plot_ly() %>%
        layout(
          title = list(text = results$message, x = 0.5, xanchor = "center"),
          xaxis = list(visible = FALSE),
          yaxis = list(visible = FALSE)
        )
    } else if (results$method == "ORA") {
      # ============ ORA Bubble Plot ============
      go_data <- results$results@result %>%
        head(20) %>%  # Show top 20 terms
        mutate(
          GeneRatio_num = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2])),
          log_pval = -log10(p.adjust),
          Description = forcats::fct_reorder(Description, GeneRatio_num)
        )
      
      plot_ly(go_data,
              x = ~GeneRatio_num,
              y = ~Description,
              type = 'scatter',
              mode = 'markers',
              marker = list(
                size = ~Count,
                sizemode = 'diameter',
                sizeref = max(go_data$Count) / 50,
                color = ~log_pval,
                colorscale = 'Viridis',
                colorbar = list(title = "-log10(adj.p)"),
                line = list(width = 1, color = 'rgba(0,0,0,0.3)')
              ),
              text = ~paste0(
                "GO Term: ", Description,
                "<br>Gene Ratio: ", GeneRatio,
                "<br>Count: ", Count,
                "<br>P-value: ", format(pvalue, scientific = TRUE, digits = 3),
                "<br>Adj. P-value: ", format(p.adjust, scientific = TRUE, digits = 3)
              ),
              hoverinfo = 'text'
      ) %>%
        layout(
          title = paste("GO Over-Representation -", input$go_ontology),
          xaxis = list(title = "Gene Ratio"),
          yaxis = list(title = ""),
          margin = list(l = 300),
          showlegend = FALSE
        )
    } else {
      # ============ GSEA Bubble Plot ============
      gsea_data <- results$results@result %>%
        head(20) %>%  # Show top 20 terms
        mutate(
          log_pval = -log10(p.adjust),
          Description = forcats::fct_reorder(Description, NES)
        )
      
      plot_ly(gsea_data,
              x = ~NES,
              y = ~Description,
              type = 'scatter',
              mode = 'markers',
              marker = list(
                size = ~setSize,
                sizemode = 'diameter',
                sizeref = max(gsea_data$setSize) / 50,
                color = ~log_pval,
                colorscale = 'Viridis',
                colorbar = list(title = "-log10(adj.p)"),
                line = list(width = 1, color = 'rgba(0,0,0,0.3)')
              ),
              text = ~paste0(
                "GO Term: ", Description,
                "<br>NES: ", round(NES, 3),
                "<br>Set Size: ", setSize,
                "<br>Enrichment Score: ", round(enrichmentScore, 3),
                "<br>P-value: ", format(pvalue, scientific = TRUE, digits = 3),
                "<br>Adj. P-value: ", format(p.adjust, scientific = TRUE, digits = 3)
              ),
              hoverinfo = 'text'
      ) %>%
        layout(
          title = paste("GSEA Results -", input$go_ontology),
          xaxis = list(title = "Normalized Enrichment Score (NES)", zeroline = TRUE),
          yaxis = list(title = ""),
          margin = list(l = 300),
          showlegend = FALSE
        ) %>%
        add_segments(x = 0, xend = 0, 
                    y = 0, yend = nrow(gsea_data) + 1,
                    line = list(color = "red", dash = "dash", width = 1),
                    showlegend = FALSE, inherit = FALSE)
    }
  })
  
  # Create summary table
  output$go_summary_table <- renderDT({
    req(go_enrichment_results())
    results <- go_enrichment_results()
    
    if (results$error) {
      # Return empty table with error message
      datatable(data.frame(Message = results$message),
                options = list(dom = 't'),
                rownames = FALSE)
    } else if (results$method == "ORA") {
      # ORA summary table
      summary_data <- results$results@result %>%
        head(20) %>%
        dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, Count) %>%
        mutate(
          pvalue = format(pvalue, scientific = TRUE, digits = 3),
          p.adjust = format(p.adjust, scientific = TRUE, digits = 3)
        )
      
      datatable(
        summary_data,
        options = list(
          pageLength = 20,
          scrollX = TRUE,
          ordering = TRUE
        ),
        rownames = FALSE,
        caption = "Top 20 GO terms (ORA)"
      )
    } else {
      # GSEA summary table
      summary_data <- results$results@result %>%
        head(20) %>%
        dplyr::select(ID, Description, setSize, enrichmentScore, NES, pvalue, p.adjust) %>%
        mutate(
          enrichmentScore = round(enrichmentScore, 3),
          NES = round(NES, 3),
          pvalue = format(pvalue, scientific = TRUE, digits = 3),
          p.adjust = format(p.adjust, scientific = TRUE, digits = 3)
        )
      
      datatable(
        summary_data,
        options = list(
          pageLength = 20,
          scrollX = TRUE,
          ordering = TRUE
        ),
        rownames = FALSE,
        colnames = c("ID", "Description", "Set Size", "ES", "NES", "P-value", "Adj. P-value"),
        caption = "Top 20 GO terms (GSEA)"
      )
    }
  })
  
  # Create full results table
  output$go_full_table <- renderDT({
    req(go_enrichment_results())
    results <- go_enrichment_results()
    
    if (results$error) {
      # Return empty table with error message
      datatable(data.frame(Message = results$message),
                options = list(dom = 't'),
                rownames = FALSE)
    } else {
      # Show all results
      full_data <- results$results@result %>%
        mutate(
          pvalue = format(pvalue, scientific = TRUE, digits = 3),
          p.adjust = format(p.adjust, scientific = TRUE, digits = 3)
        )
      
      # Format qvalue for ORA (GSEA doesn't always have qvalue)
      if ("qvalue" %in% colnames(full_data)) {
        full_data <- full_data %>%
          mutate(qvalue = format(qvalue, scientific = TRUE, digits = 3))
      }
      
      datatable(
        full_data,
        options = list(
          pageLength = 25,
          scrollX = TRUE,
          ordering = TRUE
        ),
        rownames = FALSE,
        caption = paste("Full", results$method, "Results")
      )
    }
  })
  
  # Download handler for GO results
  output$download_go_results <- downloadHandler(
    filename = function() {
      req(go_enrichment_results())
      results <- go_enrichment_results()
      paste0("GO_", results$method, "_", input$go_comparison, "_", 
             input$go_ontology, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(go_enrichment_results())
      results <- go_enrichment_results()
      
      if (!results$error) {
        write.csv(results$results@result, file, row.names = FALSE)
      }
    }
  )
  
  # ============================================================================
  # SPP1 CORRELATION TAB
  # ============================================================================
  
  # SPP1 Correlation Analysis
  spp1_cor_results <- eventReactive(input$run_spp1_cor, {
    withProgress(message = 'Calculating correlations...', value = 0, {
      tryCatch({
        # Filter metadata
        meta_sub <- metadata
        
        if (length(input$filter_cell_line_spp1) > 0) {
          meta_sub <- meta_sub %>% filter(cell_line %in% input$filter_cell_line_spp1)
        }
        if (length(input$filter_spp1_profile_spp1) > 0) {
          meta_sub <- meta_sub %>% filter(spp1_profile %in% input$filter_spp1_profile_spp1)
        }
        if (length(input$filter_cisplatine_spp1) > 0) {
          meta_sub <- meta_sub %>% filter(cisplatine %in% input$filter_cisplatine_spp1)
        }
        if (length(input$filter_comment_spp1) > 0) {
          meta_sub <- meta_sub %>% filter(comment %in% input$filter_comment_spp1)
        }
        if (length(input$filter_replicate_spp1) > 0) {
          meta_sub <- meta_sub %>% filter(replicate %in% input$filter_replicate_spp1)
        }
        
        if (nrow(meta_sub) == 0) {
          showNotification("No samples match the specified filters", type = "error")
          return(NULL)
        }
        
        # Get expression data for selected samples
        samples <- meta_sub$sample_id
        expr_sub <- expr_data_getmm[, samples, drop = FALSE]
        
        # Check if SPP1 exists
        if (!"SPP1" %in% rownames(expr_sub)) {
          showNotification("SPP1 not found in expression data", type = "error")
          return(NULL)
        }
        
        # Get SPP1 expression
        spp1_expr <- as.numeric(expr_sub["SPP1", ])
        
        # Calculate correlation for each gene with SPP1
        cor_results <- data.frame(
          gene = rownames(expr_sub),
          correlation = NA,
          pvalue = NA,
          stringsAsFactors = FALSE
        )
        
        for (i in 1:nrow(expr_sub)) {
          gene_expr <- as.numeric(expr_sub[i, ])
          
          # Skip if no variance
          if (sd(gene_expr, na.rm = TRUE) == 0 || sd(spp1_expr, na.rm = TRUE) == 0) {
            next
          }
          
          # Calculate Pearson correlation
          cor_test <- cor.test(gene_expr, spp1_expr, method = "pearson")
          cor_results$correlation[i] <- cor_test$estimate
          cor_results$pvalue[i] <- cor_test$p.value
        }
        
        # Remove NA values and SPP1 itself
        cor_results <- cor_results %>%
          filter(!is.na(correlation), !is.na(pvalue), gene != "SPP1") %>%
          mutate(
            neg_log10_pval = -log10(pvalue),
            significance = case_when(
              pvalue < input$pval_threshold_cor & correlation > input$cor_threshold ~ "Positive correlation",
              pvalue < input$pval_threshold_cor & correlation < -input$cor_threshold ~ "Negative correlation",
              pvalue < input$pval_threshold_cor ~ "Significant (weak correlation)",
              TRUE ~ "Not significant"
            )
          )
        
        list(
          data = cor_results,
          n_samples = ncol(expr_sub)
        )
      }, error = function(e) {
        showNotification(paste("Error:", e$message), type = "error")
        NULL
      })
    })
  })
  
  # Create volcano plot
  output$spp1_volcano <- renderPlotly({
    req(spp1_cor_results())
    
    cor_results <- spp1_cor_results()$data
    
    # Check if a gene is queried
    query_gene <- input$query_gene_spp1
    if (!is.null(query_gene) && query_gene != "") {
      # Highlight the queried gene
      cor_results <- cor_results %>%
        mutate(highlight = ifelse(gene == query_gene, "Queried Gene", significance))
    } else {
      cor_results <- cor_results %>%
        mutate(highlight = significance)
    }
    
    # Define colors
    color_mapping <- c(
      "Queried Gene" = "#000000",
      "Positive correlation" = "#e41a1c",
      "Negative correlation" = "#377eb8",
      "Significant (weak correlation)" = "#ff7f00",
      "Not significant" = "#999999"
    )
    
    # Define sizes
    cor_results <- cor_results %>%
      mutate(point_size = ifelse(highlight == "Queried Gene", 12, 5))
    
    p <- plot_ly(
      cor_results,
      x = ~correlation,
      y = ~neg_log10_pval,
      type = 'scatter',
      mode = 'markers',
      color = ~highlight,
      colors = color_mapping,
      marker = list(
        size = ~point_size,
        opacity = ~ifelse(highlight == "Queried Gene", 1, 0.6),
        line = list(width = ~ifelse(highlight == "Queried Gene", 2, 0),
                    color = "white")
      ),
      text = ~paste0(
        "Gene: ", gene,
        "<br>Correlation: ", round(correlation, 3),
        "<br>P-value: ", format(pvalue, scientific = TRUE, digits = 3),
        "<br>-log10(p): ", round(neg_log10_pval, 2)
      ),
      hoverinfo = 'text'
    ) %>%
      layout(
        title = paste("Gene Correlation with SPP1 (", spp1_cor_results()$n_samples, " samples)"),
        xaxis = list(
          title = "Pearson Correlation with SPP1",
          zeroline = TRUE,
          zerolinewidth = 2,
          zerolinecolor = 'black'
        ),
        yaxis = list(
          title = "-log10(P-value)"
        ),
        hovermode = 'closest',
        showlegend = TRUE,
        legend = list(title = list(text = "Significance"))
      ) %>%
      # Add threshold lines
      add_segments(
        x = input$cor_threshold, xend = input$cor_threshold,
        y = 0, yend = max(cor_results$neg_log10_pval),
        line = list(color = "black", dash = "dash", width = 1),
        showlegend = FALSE, inherit = FALSE
      ) %>%
      add_segments(
        x = -input$cor_threshold, xend = -input$cor_threshold,
        y = 0, yend = max(cor_results$neg_log10_pval),
        line = list(color = "black", dash = "dash", width = 1),
        showlegend = FALSE, inherit = FALSE
      ) %>%
      add_segments(
        x = min(cor_results$correlation), xend = max(cor_results$correlation),
        y = -log10(input$pval_threshold_cor), yend = -log10(input$pval_threshold_cor),
        line = list(color = "black", dash = "dash", width = 1),
        showlegend = FALSE, inherit = FALSE
      )
    
    p
  })
  
  # Top positively correlated genes table
  output$top_positive_genes <- renderDT({
    req(spp1_cor_results())
    
    top_pos <- spp1_cor_results()$data %>%
      filter(significance == "Positive correlation") %>%
      arrange(desc(correlation)) %>%
      head(10) %>%
      mutate(
        correlation = round(correlation, 3),
        pvalue = format(pvalue, scientific = TRUE, digits = 3),
        neg_log10_pval = round(neg_log10_pval, 2)
      ) %>%
      dplyr::select(gene, correlation, pvalue, neg_log10_pval)
    
    datatable(
      top_pos,
      options = list(
        pageLength = 10,
        dom = 't',
        ordering = FALSE
      ),
      rownames = FALSE,
      colnames = c("Gene", "Correlation", "P-value", "-log10(p)")
    )
  })
  
  # Top negatively correlated genes table
  output$top_negative_genes <- renderDT({
    req(spp1_cor_results())
    
    top_neg <- spp1_cor_results()$data %>%
      filter(significance == "Negative correlation") %>%
      arrange(correlation) %>%
      head(10) %>%
      mutate(
        correlation = round(correlation, 3),
        pvalue = format(pvalue, scientific = TRUE, digits = 3),
        neg_log10_pval = round(neg_log10_pval, 2)
      ) %>%
      dplyr::select(gene, correlation, pvalue, neg_log10_pval)
    
    datatable(
      top_neg,
      options = list(
        pageLength = 10,
        dom = 't',
        ordering = FALSE
      ),
      rownames = FALSE,
      colnames = c("Gene", "Correlation", "P-value", "-log10(p)")
    )
  })
  
  # Download handler for SPP1 correlation data
  output$download_spp1_cor <- downloadHandler(
    filename = function() {
      paste0("SPP1_correlation_results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(spp1_cor_results())
      
      # Get full correlation data
      cor_data <- spp1_cor_results()$data %>%
        mutate(
          correlation = round(correlation, 4),
          pvalue = format(pvalue, scientific = TRUE, digits = 4),
          neg_log10_pval = round(neg_log10_pval, 3)
        )
      
      write.csv(cor_data, file, row.names = FALSE)
    }
  )
  
  # Gene correlation details table
  output$gene_correlation_details <- renderDT({
    req(spp1_cor_results(), input$query_gene_spp1)
    req(input$query_gene_spp1 != "")
    
    gene <- input$query_gene_spp1
    gene_data <- spp1_cor_results()$data %>%
      filter(gene == !!gene)
    
    if (nrow(gene_data) == 0) {
      return(NULL)
    }
    
    # Create a detailed summary
    details <- data.frame(
      Metric = c("Gene", "Correlation", "P-value", "-log10(p)", "Significance"),
      Value = c(
        gene_data$gene,
        round(gene_data$correlation, 4),
        format(gene_data$pvalue, scientific = TRUE, digits = 4),
        round(gene_data$neg_log10_pval, 3),
        gene_data$significance
      )
    )
    
    datatable(
      details,
      options = list(
        dom = 't',
        ordering = FALSE,
        pageLength = 5
      ),
      rownames = FALSE
    )
  })
  
  # Gene vs SPP1 scatter plot
  output$gene_spp1_scatter <- renderPlotly({
    req(spp1_cor_results(), input$query_gene_spp1)
    req(input$query_gene_spp1 != "")
    
    gene <- input$query_gene_spp1
    
    # Get the filtered metadata from correlation analysis
    meta_sub <- metadata
    
    if (length(input$filter_cell_line_spp1) > 0) {
      meta_sub <- meta_sub %>% filter(cell_line %in% input$filter_cell_line_spp1)
    }
    if (length(input$filter_spp1_profile_spp1) > 0) {
      meta_sub <- meta_sub %>% filter(spp1_profile %in% input$filter_spp1_profile_spp1)
    }
    if (length(input$filter_cisplatine_spp1) > 0) {
      meta_sub <- meta_sub %>% filter(cisplatine %in% input$filter_cisplatine_spp1)
    }
    if (length(input$filter_comment_spp1) > 0) {
      meta_sub <- meta_sub %>% filter(comment %in% input$filter_comment_spp1)
    }
    if (length(input$filter_replicate_spp1) > 0) {
      meta_sub <- meta_sub %>% filter(replicate %in% input$filter_replicate_spp1)
    }
    
    samples <- meta_sub$sample_id
    
    # Get expression data
    if (!gene %in% rownames(expr_data_getmm)) {
      return(NULL)
    }
    
    spp1_expr <- as.numeric(expr_data_getmm["SPP1", samples])
    gene_expr <- as.numeric(expr_data_getmm[gene, samples])
    
    # Get correlation info
    gene_cor <- spp1_cor_results()$data %>%
      filter(gene == !!gene) %>%
      pull(correlation)
    
    scatter_data <- data.frame(
      spp1 = spp1_expr,
      gene = gene_expr,
      sample = samples
    )
    
    plot_ly(scatter_data,
            x = ~spp1,
            y = ~gene,
            type = 'scatter',
            mode = 'markers',
            marker = list(size = 8, opacity = 0.7),
            text = ~paste("Sample:", sample),
            hoverinfo = 'text'
    ) %>%
      add_lines(
        x = ~spp1,
        y = fitted(lm(gene ~ spp1, data = scatter_data)),
        line = list(color = 'red', dash = 'dash'),
        showlegend = FALSE,
        hoverinfo = 'skip'
      ) %>%
      layout(
        title = paste(gene, "vs SPP1 (r =", round(gene_cor, 3), ")"),
        xaxis = list(title = "SPP1 Expression"),
        yaxis = list(title = paste(gene, "Expression")),
        hovermode = 'closest'
      )
  })
  
  # ============================================================================
  # GENE EXPRESSION PROFILE TAB
  # ============================================================================
  
  # Gene profile data
  gene_profile_data <- eventReactive(input$run_gene_profile, {
    req(input$selected_gene)
    
    withProgress(message = 'Generating profile...', value = 0, {
      tryCatch({
        # Filter metadata
        meta_sub <- metadata
        
        if (length(input$filter_cell_line_profile) > 0) {
          meta_sub <- meta_sub %>% filter(cell_line %in% input$filter_cell_line_profile)
        }
        if (length(input$filter_spp1_profile_profile) > 0) {
          meta_sub <- meta_sub %>% filter(spp1_profile %in% input$filter_spp1_profile_profile)
        }
        if (length(input$filter_cisplatine_profile) > 0) {
          meta_sub <- meta_sub %>% filter(cisplatine %in% input$filter_cisplatine_profile)
        }
        if (length(input$filter_comment_profile) > 0) {
          meta_sub <- meta_sub %>% filter(comment %in% input$filter_comment_profile)
        }
        if (length(input$filter_replicate_profile) > 0) {
          meta_sub <- meta_sub %>% filter(replicate %in% input$filter_replicate_profile)
        }
        
        if (nrow(meta_sub) == 0) {
          showNotification("No samples match the specified filters", type = "error")
          return(NULL)
        }
        
        # Check if gene exists
        if (!input$selected_gene %in% rownames(expr_data_getmm)) {
          showNotification("Gene not found in expression data", type = "error")
          return(NULL)
        }
        
        # Get expression data
        samples <- meta_sub$sample_id
        gene_expr <- as.numeric(expr_data_getmm[input$selected_gene, samples])
        
        # Prepare plot data
        plot_data <- meta_sub %>%
          mutate(expression = gene_expr)
        
        list(
          data = plot_data,
          gene = input$selected_gene
        )
      }, error = function(e) {
        showNotification(paste("Error:", e$message), type = "error")
        NULL
      })
    })
  })
  
  # Create gene profile plot
  output$gene_profile_plot <- renderPlotly({
    req(gene_profile_data())
    
    plot_data <- gene_profile_data()$data
    gene_name <- gene_profile_data()$gene
    
    # Determine color variable
    color_var <- if(input$color_by_var == "same") input$group_by_var else input$color_by_var
    
    # Create plot
    p <- plot_ly(plot_data, 
                 x = as.formula(paste0("~", input$group_by_var)),
                 y = ~expression,
                 color = as.formula(paste0("~", color_var)),
                 type = 'violin',
                 box = list(visible = TRUE),
                 meanline = list(visible = TRUE),
                 text = ~paste0(
                   "Sample: ", sample_id,
                   "<br>Expression: ", round(expression, 3),
                   "<br>Cell Line: ", cell_line,
                   "<br>SPP1 Profile: ", spp1_profile,
                   "<br>Cisplatine: ", cisplatine,
                   "<br>Comment: ", comment,
                   "<br>Replicate: ", replicate
                 ),
                 hoverinfo = 'text'
    )
    
    # Add individual points if requested
    if (input$show_points) {
      p <- p %>%
        add_trace(
          type = 'box',
          boxpoints = 'all',
          jitter = 0.3,
          pointpos = 0,
          marker = list(size = 4, opacity = 0.6),
          showlegend = FALSE
        )
    }
    
    p %>%
      layout(
        title = paste("Expression Profile:", gene_name),
        xaxis = list(title = input$group_by_var),
        yaxis = list(title = "Expression (log2 TPM)"),
        hovermode = 'closest',
        showlegend = TRUE
      )
  })
  
  # Summary statistics table
  output$gene_profile_stats <- renderDT({
    req(gene_profile_data())
    
    plot_data <- gene_profile_data()$data
    
    # Calculate summary stats by group
    stats <- plot_data %>%
      group_by(!!sym(input$group_by_var)) %>%
      summarise(
        n_samples = n(),
        mean = mean(expression, na.rm = TRUE),
        median = median(expression, na.rm = TRUE),
        sd = sd(expression, na.rm = TRUE),
        min = min(expression, na.rm = TRUE),
        max = max(expression, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      mutate(
        mean = round(mean, 3),
        median = round(median, 3),
        sd = round(sd, 3),
        min = round(min, 3),
        max = round(max, 3)
      )
    
    datatable(
      stats,
      options = list(
        pageLength = 20,
        dom = 't',
        ordering = FALSE
      ),
      rownames = FALSE,
      colnames = c("Group", "N", "Mean", "Median", "SD", "Min", "Max")
    )
  })
  
  # Download handler for gene profile data
  output$download_gene_profile <- downloadHandler(
    filename = function() {
      req(input$selected_gene)
      paste0("Gene_profile_", input$selected_gene, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(gene_profile_data())
      
      # Get the full profile data with all sample information
      profile_data <- gene_profile_data()$data
      
      write.csv(profile_data, file, row.names = FALSE)
    }
  )
  
  # ============================================================================
  # PCA ANALYSIS TAB
  # ============================================================================
  
  # Initialize filter choices for PCA
  observe({
    req(exists("metadata"))
    
    updateCheckboxGroupInput(session, "filter_cell_line_pca",
                             choices = unique(metadata$cell_line),
                             selected = unique(metadata$cell_line))
    updateCheckboxGroupInput(session, "filter_spp1_profile_pca",
                             choices = unique(metadata$spp1_profile),
                             selected = unique(metadata$spp1_profile))
    updateCheckboxGroupInput(session, "filter_cisplatine_pca",
                             choices = unique(metadata$cisplatine),
                             selected = unique(metadata$cisplatine))
    updateCheckboxGroupInput(session, "filter_comment_pca",
                             choices = unique(metadata$comment),
                             selected = unique(metadata$comment))
    updateCheckboxGroupInput(session, "filter_replicate_pca",
                             choices = as.character(unique(metadata$replicate)),
                             selected = as.character(unique(metadata$replicate)))
  })
  
  # Perform PCA
  pca_results <- eventReactive(input$run_pca, {
    withProgress(message = 'Running PCA...', value = 0, {
      tryCatch({
        incProgress(0.1, detail = "Filtering samples...")
        
        # Filter metadata
        meta_sub <- metadata
        
        if (!is.null(input$filter_cell_line_pca) && length(input$filter_cell_line_pca) > 0) {
          meta_sub <- meta_sub %>% filter(cell_line %in% input$filter_cell_line_pca)
        }
        if (!is.null(input$filter_spp1_profile_pca) && length(input$filter_spp1_profile_pca) > 0) {
          meta_sub <- meta_sub %>% filter(spp1_profile %in% input$filter_spp1_profile_pca)
        }
        if (!is.null(input$filter_cisplatine_pca) && length(input$filter_cisplatine_pca) > 0) {
          meta_sub <- meta_sub %>% filter(cisplatine %in% input$filter_cisplatine_pca)
        }
        if (!is.null(input$filter_comment_pca) && length(input$filter_comment_pca) > 0) {
          meta_sub <- meta_sub %>% filter(comment %in% input$filter_comment_pca)
        }
        if (!is.null(input$filter_replicate_pca) && length(input$filter_replicate_pca) > 0) {
          meta_sub <- meta_sub %>% filter(as.character(replicate) %in% input$filter_replicate_pca)
        }
        
        if (nrow(meta_sub) < 3) {
          return(list(error = TRUE, message = "Need at least 3 samples for PCA"))
        }
        
        incProgress(0.2, detail = "Selecting genes...")
        
        # Get expression data for selected samples
        expr_sub <- expr_data_getmm[, meta_sub$sample_id, drop = FALSE]
        
        # Select most variable genes
        gene_vars <- apply(expr_sub, 1, var, na.rm = TRUE)
        top_genes <- names(sort(gene_vars, decreasing = TRUE)[1:min(input$pca_top_genes, nrow(expr_sub))])
        expr_pca <- expr_sub[top_genes, ]
        
        # Remove genes with any NA or zero variance
        expr_pca <- expr_pca[complete.cases(expr_pca), ]
        gene_vars_filtered <- apply(expr_pca, 1, var)
        expr_pca <- expr_pca[gene_vars_filtered > 0, ]
        
        if (nrow(expr_pca) < 10) {
          return(list(error = TRUE, message = "Not enough variable genes for PCA"))
        }
        
        incProgress(0.3, detail = "Running PCA...")
        
        # Transpose for PCA (samples as rows, genes as columns)
        expr_t <- t(expr_pca)
        
        # Scale if requested
        if (input$pca_scale) {
          pca_result <- prcomp(expr_t, center = TRUE, scale. = TRUE)
        } else {
          pca_result <- prcomp(expr_t, center = TRUE, scale. = FALSE)
        }
        
        incProgress(0.3, detail = "Processing results...")
        
        # Calculate variance explained
        var_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100
        
        # Get PC scores
        pca_scores <- as.data.frame(pca_result$x)
        pca_scores$sample_id <- rownames(pca_scores)
        
        # Merge with metadata
        pca_data <- merge(pca_scores, meta_sub, by = "sample_id")
        
        list(
          error = FALSE,
          pca_result = pca_result,
          pca_data = pca_data,
          var_explained = var_explained,
          n_samples = nrow(meta_sub),
          n_genes = nrow(expr_pca)
        )
        
      }, error = function(e) {
        return(list(error = TRUE, message = paste("Error:", e$message)))
      })
    })
  })
  
  # PCA plot
  output$pca_plot <- renderPlotly({
    req(pca_results())
    results <- pca_results()
    
    if (results$error) {
      plot_ly() %>%
        layout(
          title = list(text = results$message, x = 0.5, xanchor = "center"),
          xaxis = list(visible = FALSE),
          yaxis = list(visible = FALSE)
        )
    } else {
      pca_data <- results$pca_data
      var_exp <- results$var_explained
      
      # Select color and shape variables
      color_var <- input$pca_color_by
      shape_var <- if (input$pca_shape_by == "none") NULL else input$pca_shape_by
      
      # Create plot
      p <- plot_ly(pca_data, 
                   x = ~PC1, y = ~PC2,
                   type = 'scatter', mode = 'markers',
                   color = ~get(color_var),
                   symbol = if (!is.null(shape_var)) ~get(shape_var) else NULL,
                   marker = list(size = 10),
                   text = ~paste("Sample:", sample_id,
                                "<br>Cell Line:", cell_line,
                                "<br>SPP1:", spp1_profile,
                                "<br>Cisplatine:", cisplatine,
                                "<br>Comment:", comment,
                                "<br>Replicate:", replicate),
                   hoverinfo = 'text') %>%
        layout(
          title = paste0("PCA Analysis (", results$n_samples, " samples, ", results$n_genes, " genes)"),
          xaxis = list(title = paste0("PC1 (", round(var_exp[1], 1), "%)")),
          yaxis = list(title = paste0("PC2 (", round(var_exp[2], 1), "%)")),
          hovermode = 'closest'
        )
      
      p
    }
  })
  
  # Scree plot
  output$scree_plot <- renderPlotly({
    req(pca_results())
    results <- pca_results()
    
    if (results$error) {
      return(NULL)
    }
    
    var_exp <- results$var_explained
    n_pcs <- min(20, length(var_exp))
    
    scree_data <- data.frame(
      PC = paste0("PC", 1:n_pcs),
      Variance = var_exp[1:n_pcs],
      PC_num = 1:n_pcs
    )
    
    plot_ly(scree_data, x = ~PC_num, y = ~Variance,
            type = 'scatter', mode = 'lines+markers',
            marker = list(size = 8),
            line = list(width = 2)) %>%
      layout(
        title = "Scree Plot",
        xaxis = list(title = "Principal Component", dtick = 1),
        yaxis = list(title = "Variance Explained (%)"),
        hovermode = 'closest'
      )
  })
  
  # Variance table
  output$variance_table <- renderDT({
    req(pca_results())
    results <- pca_results()
    
    if (results$error) {
      return(NULL)
    }
    
    var_exp <- results$var_explained
    n_pcs <- min(10, length(var_exp))
    
    var_table <- data.frame(
      PC = paste0("PC", 1:n_pcs),
      `Variance Explained (%)` = round(var_exp[1:n_pcs], 2),
      `Cumulative (%)` = round(cumsum(var_exp)[1:n_pcs], 2),
      check.names = FALSE
    )
    
    datatable(var_table, 
              options = list(dom = 't', ordering = FALSE),
              rownames = FALSE)
  })
  
  # Loadings table
  output$loadings_table <- renderDT({
    req(pca_results())
    results <- pca_results()
    
    if (results$error) {
      return(NULL)
    }
    
    pc_num <- as.numeric(input$pc_for_loadings)
    loadings <- results$pca_result$rotation[, pc_num]
    
    # Get top 50 genes by absolute loading
    top_loadings <- sort(abs(loadings), decreasing = TRUE)[1:min(50, length(loadings))]
    top_genes <- names(top_loadings)
    
    loadings_table <- data.frame(
      Gene = top_genes,
      Loading = round(loadings[top_genes], 4),
      `Abs Loading` = round(abs(loadings[top_genes]), 4),
      check.names = FALSE
    ) %>%
      arrange(desc(`Abs Loading`))
    
    datatable(loadings_table,
              options = list(pageLength = 20, dom = 'tp'),
              rownames = FALSE)
  })
  
  # Download PCA data
  output$download_pca_plot <- downloadHandler(
    filename = function() {
      paste0("PCA_analysis_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(pca_results())
      results <- pca_results()
      if (!results$error) {
        write.csv(results$pca_data, file, row.names = FALSE)
      }
    }
  )
}

# ============================================================================
# RUN APP
# ============================================================================

shinyApp(ui = ui, server = server)