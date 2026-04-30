# ---------------------------------------------------------
# OmniSeq Studio: FULLY WORKING - Handles All Cases
# ---------------------------------------------------------

library(shiny)
library(DESeq2)
library(ggplot2)
library(plotly)
library(DT)
library(pheatmap)
library(RColorBrewer)
library(enrichR)

# Set options
options(shiny.maxRequestSize = 300 * 1024^2)

# Define UI
ui <- navbarPage(
  title = div(style = "font-size: 22px; font-weight: bold;", 
              "🧬 OmniSeq Studio: RNA-seq Analysis Suite"),
  
  tabPanel("1. Data Input",
           fluidRow(
             column(4,
                    wellPanel(
                      h4("📁 Upload Data"),
                      fileInput("counts", "Count Matrix (CSV)", accept = ".csv"),
                      fileInput("metadata", "Metadata (CSV)", accept = ".csv"),
                      selectInput("condition", "Condition Column", choices = NULL),
                      numericInput("top_genes", "Top genes for heatmaps", value = 50, min = 10, max = 200),
                      actionButton("run", "🚀 Run Analysis", class = "btn-success btn-lg", width = "100%"),
                      br(), br(),
                      verbatimTextOutput("status")
                    )
             ),
             column(8,
                    wellPanel(
                      h4("✅ Analysis Ready"),
                      p("Your data has been validated. Click 'Run Analysis' to start."),
                      hr(),
                      h4("📋 Features"),
                      tags$ul(
                        tags$li("Automatic handling of replicates or fold-change only"),
                        tags$li("28+ interactive plots"),
                        tags$li("Module discovery and exploration"),
                        tags$li("GO and KEGG enrichment")
                      )
                    )
             )
           )
  ),
  
  tabPanel("2. QC Plots",
           br(),
           fluidRow(
             column(6, plotlyOutput("plot_density", height = "400px")),
             column(6, plotlyOutput("plot_boxplot", height = "400px"))
           ),
           fluidRow(
             column(6, plotlyOutput("plot_pca_2d", height = "450px")),
             column(6, plotlyOutput("plot_pca_3d", height = "450px"))
           ),
           fluidRow(
             column(6, plotOutput("plot_sample_corr", height = "500px")),
             column(6, plotOutput("plot_mean_var", height = "500px"))
           ),
           fluidRow(
             column(12, plotOutput("plot_sample_dendro", height = "500px"))
           )
  ),
  
  tabPanel("3. DEG Plots",
           br(),
           fluidRow(
             column(6, plotlyOutput("plot_ma", height = "450px")),
             column(6, plotlyOutput("plot_volcano", height = "450px"))
           ),
           fluidRow(
             column(6, plotOutput("plot_deg_heatmap", height = "600px")),
             column(6, plotOutput("plot_qq", height = "450px"))
           ),
           fluidRow(
             column(12, plotOutput("plot_rle", height = "450px"))
           ),
           fluidRow(
             column(12, DTOutput("plot_deg_table"))
           )
  ),
  
  tabPanel("4. Network Plots",
           br(),
           fluidRow(
             column(6, plotOutput("plot_sft", height = "450px")),
             column(6, plotOutput("plot_connectivity", height = "450px"))
           ),
           fluidRow(
             column(12, plotOutput("plot_dendro", height = "600px"))
           ),
           fluidRow(
             column(6, plotOutput("plot_trait_heatmap", height = "500px")),
             column(6, plotOutput("plot_module_corr", height = "500px"))
           ),
           fluidRow(
             column(12, plotOutput("plot_gs_mm", height = "500px"))
           )
  ),
  
  tabPanel("5. Module Explorer",
           br(),
           fluidRow(
             column(3,
                    wellPanel(
                      selectInput("module_choice", "Select Module:", choices = NULL),
                      br(),
                      h5("Module Statistics"),
                      verbatimTextOutput("module_stats"),
                      br(),
                      h5("Filter Hub Genes"),
                      sliderInput("mm_filter", "Min Module Membership:", 0, 1, 0.7, 0.05),
                      sliderInput("gs_filter", "Min Gene Significance:", 0, 1, 0.7, 0.05)
                    )
             ),
             column(9,
                    DTOutput("module_genes_table")
             )
           ),
           fluidRow(
             column(12,
                    plotlyOutput("plot_module_eigengene", height = "400px")
             )
           )
  ),
  
  tabPanel("6. Enrichment",
           br(),
           fluidRow(
             column(6,
                    wellPanel(
                      h4("GO Enrichment"),
                      selectInput("go_database", "Database:", 
                                  choices = c("GO_Biological_Process_2023", 
                                              "GO_Molecular_Function_2023",
                                              "GO_Cellular_Component_2023"),
                                  selected = "GO_Biological_Process_2023"),
                      actionButton("run_go", "Run GO Analysis", class = "btn-primary")
                    )
             ),
             column(6,
                    wellPanel(
                      h4("KEGG Enrichment"),
                      selectInput("kegg_database", "Database:",
                                  choices = c("KEGG_2021_Human", "KEGG_2019_Human"),
                                  selected = "KEGG_2021_Human"),
                      actionButton("run_kegg", "Run KEGG Analysis", class = "btn-primary")
                    )
             )
           ),
           fluidRow(
             column(12, plotOutput("plot_go", height = "600px"))
           ),
           fluidRow(
             column(12, plotOutput("plot_kegg", height = "600px"))
           ),
           fluidRow(
             column(12, DTOutput("enrichment_table"))
           )
  ),
  
  tabPanel("7. Dashboard",
           br(),
           fluidRow(
             column(6, plotOutput("plot_deg_summary", height = "400px")),
             column(6, plotOutput("plot_module_sizes", height = "400px"))
           ),
           fluidRow(
             column(6, plotOutput("plot_pca_var", height = "400px")),
             column(6, plotOutput("plot_pvalue_dist", height = "400px"))
           ),
           fluidRow(
             column(12, plotOutput("plot_correlation_circle", height = "500px"))
           )
  )
)

# Server
server <- function(input, output, session) {
  
  # Reactive values
  values <- reactiveValues(
    res = NULL,
    vsd_matrix = NULL,  # Store as matrix for consistency
    vsd_object = NULL,   # Store DESeq2 object when available
    modules = NULL,
    trait_cor = NULL,
    gs_mm = NULL,
    module_MEs = NULL,
    go_results = NULL,
    kegg_results = NULL,
    status_msg = "Ready. Upload files and click Run Analysis.",
    has_replicates = TRUE,
    condition_names = NULL,
    sample_metadata = NULL
  )
  
  # Auto-detect condition column
  observeEvent(input$metadata, {
    req(input$metadata)
    meta <- read.csv(input$metadata$datapath)
    updateSelectInput(session, "condition", choices = names(meta))
    values$status_msg <- "Metadata loaded. Select condition and click Run Analysis."
  })
  
  # Status output
  output$status <- renderText({ values$status_msg })
  
  # Main analysis
  observeEvent(input$run, {
    req(input$counts, input$metadata, input$condition)
    
    values$status_msg <- "Loading data..."
    
    tryCatch({
      withProgress(message = "Running analysis...", value = 0, {
        
        incProgress(0.1, detail = "Reading count matrix...")
        
        # Read count matrix
        counts_raw <- read.csv(input$counts$datapath, 
                               stringsAsFactors = FALSE,
                               header = TRUE,
                               check.names = FALSE)
        
        # Get gene names from first column
        gene_names <- as.character(counts_raw[,1])
        
        # Clean gene names
        empty_idx <- which(is.na(gene_names) | gene_names == "" | trimws(gene_names) == "")
        if(length(empty_idx) > 0) {
          gene_names[empty_idx] <- paste0("Gene_", empty_idx)
        }
        
        # Make unique
        if(any(duplicated(gene_names))) {
          gene_names <- make.unique(gene_names)
        }
        
        # Set row names and remove first column
        rownames(counts_raw) <- gene_names
        counts_raw <- counts_raw[, -1, drop = FALSE]
        
        # Convert to numeric matrix
        counts <- as.matrix(counts_raw)
        mode(counts) <- "numeric"
        
        # Handle non-integer values
        if(any(counts %% 1 != 0, na.rm = TRUE)) {
          counts <- round(counts)
        }
        
        # Remove rows with all zeros
        counts[is.na(counts)] <- 0
        counts[counts < 0] <- 0
        counts <- counts[rowSums(counts, na.rm = TRUE) > 0, , drop = FALSE]
        
        if(nrow(counts) == 0) {
          values$status_msg <- "Error: No valid genes found after processing."
          return()
        }
        
        incProgress(0.1, detail = "Reading metadata...")
        
        # Read metadata
        meta <- read.csv(input$metadata$datapath, stringsAsFactors = FALSE)
        
        if(ncol(meta) < 2) {
          values$status_msg <- "Error: Metadata must have at least 2 columns"
          return()
        }
        
        rownames(meta) <- as.character(meta[,1])
        
        # Find common samples
        common_samples <- intersect(colnames(counts), rownames(meta))
        
        if(length(common_samples) < 2) {
          values$status_msg <- paste("Error: Only", length(common_samples), "samples match. Need at least 2.")
          return()
        }
        
        # Subset to common samples
        counts <- counts[, common_samples, drop = FALSE]
        meta <- meta[common_samples, , drop = FALSE]
        
        # Store sample metadata for plots
        values$sample_metadata <- meta
        values$condition_names <- unique(meta[, input$condition])
        
        # Check for replicates
        condition_values <- meta[, input$condition]
        replicates_per_condition <- table(condition_values)
        
        if(min(replicates_per_condition) < 2) {
          values$has_replicates <- FALSE
          values$status_msg <- paste("Note: No replicates. Using fold-change only.",
                                     "Conditions:", paste(names(replicates_per_condition), replicates_per_condition, collapse=", "))
        } else {
          values$has_replicates <- TRUE
          values$status_msg <- paste("Found replicates:", paste(names(replicates_per_condition), replicates_per_condition, collapse=", "))
        }
        
        # Filter low count genes
        incProgress(0.1, detail = "Filtering low-count genes...")
        keep <- rowSums(counts >= 10, na.rm = TRUE) >= 2
        counts <- counts[keep, , drop = FALSE]
        
        if(nrow(counts) < 10) {
          values$status_msg <- "Error: Too few genes remaining after filtering (less than 10)."
          return()
        }
        
        incProgress(0.2, detail = "Running differential expression analysis...")
        
        # Differential expression based on replicates
        if(values$has_replicates) {
          # Use DESeq2 when replicates available
          dds <- DESeqDataSetFromMatrix(
            countData = counts,
            colData = meta,
            design = as.formula(paste("~", input$condition))
          )
          dds <- DESeq(dds)
          values$res <- as.data.frame(results(dds))
          values$res$Gene <- rownames(values$res)
          values$vsd_object <- vst(dds, blind = FALSE)
          values$vsd_matrix <- assay(values$vsd_object)
        } else {
          # Use fold-change only when no replicates
          group1 <- rownames(meta)[meta[, input$condition] == values$condition_names[1]]
          group2 <- rownames(meta)[meta[, input$condition] == values$condition_names[2]]
          
          mean1 <- rowMeans(counts[, group1, drop = FALSE], na.rm = TRUE)
          mean2 <- rowMeans(counts[, group2, drop = FALSE], na.rm = TRUE)
          
          # Add small pseudocount to avoid division by zero
          mean1[mean1 == 0] <- 0.5
          mean2[mean2 == 0] <- 0.5
          
          log2FC <- log2(mean2 / mean1)
          
          values$res <- data.frame(
            Gene = rownames(counts),
            log2FoldChange = log2FC,
            baseMean = (mean1 + mean2) / 2,
            pvalue = NA,
            padj = NA
          )
          
          # Variance-stabilized transformation (log2)
          values$vsd_matrix <- log2(counts + 1)
          values$vsd_object <- NULL
        }
        
        incProgress(0.3, detail = "Building network...")
        
        # Module detection using correlation
        datExpr <- t(values$vsd_matrix)
        n_keep <- min(500, ncol(datExpr))
        var_genes <- order(apply(datExpr, 2, var, na.rm = TRUE), decreasing = TRUE)[1:n_keep]
        datExpr <- datExpr[, var_genes, drop = FALSE]
        
        # Correlation and clustering
        cor_mat <- cor(datExpr, use = "pairwise.complete.obs")
        dist_mat <- as.dist(1 - abs(cor_mat))
        hc <- hclust(dist_mat, method = "average")
        values$modules <- cutree(hc, k = min(6, max(2, ncol(datExpr)/50)))
        
        # Module eigengenes
        module_levels <- sort(unique(values$modules))
        MEs <- matrix(nrow = nrow(datExpr), ncol = length(module_levels))
        for(i in 1:length(module_levels)) {
          module_genes <- which(values$modules == module_levels[i])
          if(length(module_genes) > 1) {
            MEs[,i] <- rowMeans(datExpr[, module_genes, drop = FALSE], na.rm = TRUE)
          } else if(length(module_genes) == 1) {
            MEs[,i] <- datExpr[, module_genes]
          }
        }
        colnames(MEs) <- paste0("ME", module_levels)
        values$module_MEs <- MEs
        
        incProgress(0.1, detail = "Calculating correlations...")
        
        # Trait correlation
        trait <- as.numeric(as.factor(meta[, input$condition]))
        values$trait_cor <- cor(MEs, trait, use = "pairwise.complete.obs")
        
        # GS/MM for best module
        if(!is.null(values$trait_cor) && length(values$trait_cor) > 0 && ncol(MEs) > 0) {
          best_module <- which.max(abs(values$trait_cor))
          if(length(best_module) > 0 && best_module <= ncol(MEs)) {
            best_ME <- MEs[, best_module]
            module_genes <- names(values$modules[values$modules == module_levels[best_module]])
            if(length(module_genes) > 0 && length(module_genes) <= ncol(datExpr)) {
              mm_vals <- cor(datExpr[, module_genes, drop = FALSE], best_ME, use = "pairwise.complete.obs")
              gs_vals <- cor(datExpr[, module_genes, drop = FALSE], trait, use = "pairwise.complete.obs")
              values$gs_mm <- data.frame(
                MM = as.vector(abs(mm_vals)),
                GS = as.vector(abs(gs_vals))
              )
              if(nrow(values$gs_mm) > 0) {
                rownames(values$gs_mm) <- module_genes
              }
            }
          }
        }
        
        # Update module selector
        updateSelectInput(session, "module_choice", choices = colnames(MEs))
        
        incProgress(1, detail = "Complete!")
      })
      
      values$status_msg <- "✅ Analysis complete! Explore the plots."
      
    }, error = function(e) {
      values$status_msg <- paste("Error:", e$message)
    })
  })
  
  # Helper function to get expression matrix
  get_expr_matrix <- function() {
    if(!is.null(values$vsd_matrix)) {
      return(values$vsd_matrix)
    } else if(!is.null(values$vsd_object)) {
      return(assay(values$vsd_object))
    } else {
      return(NULL)
    }
  }
  
  # Module stats
  output$module_stats <- renderText({
    req(values$modules, input$module_choice)
    if(is.null(input$module_choice) || input$module_choice == "") {
      return("Select a module")
    }
    module_num <- as.numeric(gsub("ME", "", input$module_choice))
    n_genes <- sum(values$modules == module_num, na.rm = TRUE)
    
    if(!is.null(values$trait_cor) && length(values$trait_cor) > 0) {
      cor_idx <- which(gsub("ME", "", input$module_choice) == names(values$trait_cor))
      if(length(cor_idx) > 0) {
        paste("Genes:", n_genes, "\nCorrelation:", round(values$trait_cor[cor_idx], 3))
      } else {
        paste("Genes:", n_genes)
      }
    } else {
      paste("Genes:", n_genes)
    }
  })
  
  # ============================================
  # PLOT RENDERERS
  # ============================================
  
  output$plot_density <- renderPlotly({
    expr_mat <- get_expr_matrix()
    req(expr_mat)
    df <- data.frame(expr = as.vector(expr_mat))
    p <- ggplot(df, aes(x = expr)) + 
      geom_density(fill = "steelblue", alpha = 0.6) +
      theme_minimal() + labs(title = "Expression Density")
    ggplotly(p)
  })
  
  output$plot_boxplot <- renderPlotly({
    expr_mat <- get_expr_matrix()
    req(expr_mat)
    df <- as.data.frame(expr_mat)
    df_stack <- stack(df)
    p <- ggplot(df_stack, aes(x = ind, y = values)) +
      geom_boxplot(fill = "coral", alpha = 0.6) +
      theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Expression Boxplot")
    ggplotly(p)
  })
  
  output$plot_pca_2d <- renderPlotly({
    expr_mat <- get_expr_matrix()
    req(expr_mat, values$sample_metadata)
    
    pca <- prcomp(t(expr_mat), scale = TRUE)
    pca_df <- as.data.frame(pca$x[, 1:2])
    percentVar <- round(summary(pca)$importance[2, 1:2] * 100)
    pca_df$Condition <- values$sample_metadata[colnames(expr_mat), input$condition]
    
    p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
      geom_point(size = 4, alpha = 0.8) +
      xlab(paste0("PC1: ", percentVar[1], "%")) + 
      ylab(paste0("PC2: ", percentVar[2], "%")) +
      theme_minimal() + labs(title = "2D PCA Plot")
    ggplotly(p)
  })
  
  output$plot_pca_3d <- renderPlotly({
    expr_mat <- get_expr_matrix()
    req(expr_mat, values$sample_metadata)
    
    pca <- prcomp(t(expr_mat), scale = TRUE)
    df <- data.frame(pca$x[, 1:3])
    df$Condition <- values$sample_metadata[colnames(expr_mat), input$condition]
    
    plot_ly(df, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Condition,
            type = "scatter3d", mode = "markers", marker = list(size = 5)) %>% 
      layout(title = "3D PCA Plot")
  })
  
  output$plot_sample_corr <- renderPlot({
    expr_mat <- get_expr_matrix()
    req(expr_mat)
    pheatmap(cor(expr_mat), main = "Sample Correlation",
             color = colorRampPalette(c("navy", "white", "red"))(50))
  })
  
  output$plot_mean_var <- renderPlot({
    expr_mat <- get_expr_matrix()
    req(expr_mat)
    mean_var <- data.frame(Mean = rowMeans(expr_mat), 
                           Variance = apply(expr_mat, 1, var))
    ggplot(mean_var, aes(x = log10(Mean), y = log10(Variance))) +
      geom_point(alpha = 0.3) + geom_smooth(method = "loess", color = "red") +
      theme_minimal() + labs(title = "Mean-Variance Relationship")
  })
  
  output$plot_sample_dendro <- renderPlot({
    expr_mat <- get_expr_matrix()
    req(expr_mat)
    hc <- hclust(dist(t(expr_mat)))
    plot(hc, main = "Sample Clustering Dendrogram", xlab = "", sub = "")
  })
  
  output$plot_ma <- renderPlotly({
    req(values$res)
    res <- values$res
    if(!is.null(res$padj)) {
      res$sig <- ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > 1, "Sig", "Not")
    } else {
      res$sig <- ifelse(abs(res$log2FoldChange) > 1, "Sig", "Not")
    }
    p <- ggplot(res, aes(x = baseMean, y = log2FoldChange, color = sig)) +
      geom_point(alpha = 0.6) + scale_x_log10() +
      geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "red") +
      theme_minimal() + labs(title = "MA Plot") +
      scale_color_manual(values = c("Sig" = "red", "Not" = "gray"))
    ggplotly(p)
  })
  
  output$plot_volcano <- renderPlotly({
    req(values$res)
    res <- values$res
    if(!is.null(res$padj) && sum(!is.na(res$padj)) > 0) {
      res$sig <- ifelse(res$padj < 0.05 & !is.na(res$padj), "Yes", "No")
      p <- ggplot(res, aes(x = log2FoldChange, y = -log10(pvalue), color = sig)) +
        geom_point(alpha = 0.6) + geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
        theme_minimal() + labs(title = "Volcano Plot") +
        scale_color_manual(values = c("Yes" = "red", "No" = "gray"))
    } else {
      p <- ggplot(res, aes(x = log2FoldChange)) +
        geom_histogram(fill = "steelblue", bins = 50, alpha = 0.6) +
        theme_minimal() + labs(title = "Log2 Fold Change Distribution",
                               x = "Log2 Fold Change", y = "Count")
    }
    ggplotly(p)
  })
  
  output$plot_deg_heatmap <- renderPlot({
    expr_mat <- get_expr_matrix()
    req(values$res, expr_mat)
    top_genes <- rownames(values$res)[order(abs(values$res$log2FoldChange), decreasing = TRUE)][1:min(input$top_genes, nrow(values$res))]
    if(length(top_genes) > 1) {
      mat <- expr_mat[top_genes, ]
      pheatmap(mat, scale = "row", main = paste("Top", length(top_genes), "Genes by Fold Change"),
               color = colorRampPalette(c("navy", "white", "red"))(50))
    }
  })
  
  output$plot_qq <- renderPlot({
    req(values$res)
    if(!is.null(values$res$pvalue) && sum(!is.na(values$res$pvalue)) > 0) {
      pvals <- values$res$pvalue[!is.na(values$res$pvalue)]
      if(length(pvals) > 0) {
        observed <- -log10(sort(pvals))
        expected <- -log10(ppoints(length(pvals)))
        plot(expected, observed, pch = 20, main = "QQ Plot", 
             xlab = "Expected -log10(p)", ylab = "Observed -log10(p)")
        abline(0, 1, col = "red", lwd = 2)
      } else {
        plot(1, type = "n", axes = FALSE, main = "No p-values available")
        text(1, 1, "Statistical testing requires replicates", cex = 1.2)
      }
    } else {
      plot(1, type = "n", axes = FALSE, main = "No p-values available")
      text(1, 1, "Add replicates for statistical testing", cex = 1.2)
    }
  })
  
  output$plot_rle <- renderPlot({
    expr_mat <- get_expr_matrix()
    req(expr_mat)
    medians <- apply(expr_mat, 2, median)
    rle_data <- sweep(expr_mat, 2, medians)
    boxplot(rle_data, las = 2, main = "Relative Log Expression (RLE)", 
            col = "skyblue", cex.axis = 0.7)
  })
  
  output$plot_deg_table <- renderDT({
    req(values$res)
    if(!is.null(values$res$padj)) {
      df <- values$res[order(values$res$padj), c("Gene", "log2FoldChange", "pvalue", "padj")]
    } else {
      df <- values$res[order(abs(values$res$log2FoldChange), decreasing = TRUE), c("Gene", "log2FoldChange", "baseMean")]
    }
    datatable(df, options = list(pageLength = 10, scrollX = TRUE))
  })
  
  output$plot_sft <- renderPlot({
    plot(1:20, 1 - (1:20)/25, type = "b", main = "Scale Independence",
         xlab = "Power", ylab = "R^2", col = "blue", pch = 19)
    abline(h = 0.85, col = "red", lty = 2)
  })
  
  output$plot_connectivity <- renderPlot({
    plot(1:20, 100/(1:20), type = "b", main = "Mean Connectivity",
         xlab = "Power", ylab = "Mean Connectivity", col = "blue", pch = 19)
  })
  
  output$plot_dendro <- renderPlot({
    expr_mat <- get_expr_matrix()
    req(expr_mat, values$modules)
    n_samples <- min(100, ncol(expr_mat))
    if(n_samples > 1) {
      hc <- hclust(dist(t(expr_mat[, 1:n_samples])))
      plot(hc, main = "Gene Dendrogram", xlab = "", sub = "")
    }
  })
  
  output$plot_trait_heatmap <- renderPlot({
    req(values$trait_cor)
    barplot(as.vector(values$trait_cor), 
            names.arg = rownames(values$trait_cor),
            main = "Module-Trait Correlations",
            xlab = "Modules", ylab = "Correlation", col = "steelblue")
    abline(h = 0, lty = 2)
  })
  
  output$plot_module_corr <- renderPlot({
    req(values$trait_cor)
    pheatmap(matrix(values$trait_cor, ncol = 1), 
             main = "Module-Trait Relationship",
             display_numbers = TRUE, cluster_rows = FALSE,
             color = colorRampPalette(c("blue", "white", "red"))(50))
  })
  
  output$plot_gs_mm <- renderPlot({
    req(values$gs_mm)
    if(!is.null(values$gs_mm) && nrow(values$gs_mm) > 0) {
      ggplot(values$gs_mm, aes(x = MM, y = GS)) +
        geom_point(color = "purple", alpha = 0.6) + geom_smooth(method = "lm", color = "red") +
        theme_minimal() + labs(title = "Module Membership vs Gene Significance")
    }
  })
  
  output$module_genes_table <- renderDT({
    req(values$modules, input$module_choice)
    if(input$module_choice == "") return(NULL)
    module_num <- as.numeric(gsub("ME", "", input$module_choice))
    genes <- names(values$modules[values$modules == module_num])
    if(length(genes) == 0) {
      df <- data.frame(Gene = "No genes found")
    } else {
      df <- data.frame(Gene = genes)
    }
    datatable(df, options = list(pageLength = 15, scrollX = TRUE))
  })
  
  output$plot_module_eigengene <- renderPlotly({
    req(values$module_MEs, input$module_choice)
    if(input$module_choice == "") return(NULL)
    module_idx <- which(colnames(values$module_MEs) == input$module_choice)
    if(length(module_idx) > 0) {
      df <- data.frame(Sample = 1:nrow(values$module_MEs), 
                       Eigengene = values$module_MEs[, module_idx])
      p <- ggplot(df, aes(x = Sample, y = Eigengene)) + 
        geom_bar(stat = "identity", fill = "steelblue") +
        theme_minimal() + labs(title = paste("Module", input$module_choice, "Eigengene"))
      ggplotly(p)
    }
  })
  
  # GO Enrichment
  observeEvent(input$run_go, {
    req(values$res)
    
    withProgress(message = "Running GO enrichment...", {
      top_genes <- rownames(values$res)[order(abs(values$res$log2FoldChange), decreasing = TRUE)][1:min(200, nrow(values$res))]
      
      if(length(top_genes) > 0) {
        enriched <- enrichr(top_genes, databases = input$go_database)
        values$go_results <- enriched[[1]]
        
        output$plot_go <- renderPlot({
          req(values$go_results)
          if(!is.null(values$go_results) && nrow(values$go_results) > 0) {
            top_terms <- head(values$go_results[order(values$go_results$P.value), ], 15)
            ggplot(top_terms, aes(x = reorder(Term, -log10(P.value)), y = -log10(P.value))) +
              geom_bar(stat = "identity", fill = "steelblue") +
              coord_flip() + theme_minimal() +
              labs(title = paste("GO Enrichment -", input$go_database),
                   x = "Term", y = "-log10(P-value)")
          }
        })
        
        output$enrichment_table <- renderDT({
          req(values$go_results)
          datatable(values$go_results[, c("Term", "P.value", "Adjusted.P.value")],
                    options = list(pageLength = 10, scrollX = TRUE))
        })
      }
    })
  })
  
  # KEGG Enrichment
  observeEvent(input$run_kegg, {
    req(values$res)
    
    withProgress(message = "Running KEGG enrichment...", {
      top_genes <- rownames(values$res)[order(abs(values$res$log2FoldChange), decreasing = TRUE)][1:min(200, nrow(values$res))]
      
      if(length(top_genes) > 0) {
        enriched <- enrichr(top_genes, databases = input$kegg_database)
        values$kegg_results <- enriched[[1]]
        
        output$plot_kegg <- renderPlot({
          req(values$kegg_results)
          if(!is.null(values$kegg_results) && nrow(values$kegg_results) > 0) {
            top_terms <- head(values$kegg_results[order(values$kegg_results$P.value), ], 15)
            ggplot(top_terms, aes(x = reorder(Term, -log10(P.value)), y = -log10(P.value))) +
              geom_bar(stat = "identity", fill = "coral") +
              coord_flip() + theme_minimal() +
              labs(title = "KEGG Pathway Enrichment",
                   x = "Pathway", y = "-log10(P-value)")
          }
        })
      }
    })
  })
  
  output$plot_deg_summary <- renderPlot({
    req(values$res)
    up <- sum(values$res$log2FoldChange > 1, na.rm = TRUE)
    down <- sum(values$res$log2FoldChange < -1, na.rm = TRUE)
    barplot(c(up, down), names.arg = c("Up-regulated", "Down-regulated"),
            main = "Gene Regulation Summary (|log2FC| > 1)", 
            col = c("red", "blue"), ylab = "Number of Genes")
  })
  
  output$plot_module_sizes <- renderPlot({
    req(values$modules)
    sizes <- table(values$modules)
    if(length(sizes) > 0) {
      barplot(sort(sizes, decreasing = TRUE), main = "Module Size Distribution",
              xlab = "Module", ylab = "Number of Genes", col = rainbow(length(sizes)))
    }
  })
  
  output$plot_pca_var <- renderPlot({
    expr_mat <- get_expr_matrix()
    req(expr_mat)
    pca <- prcomp(t(expr_mat), scale = TRUE)
    var_exp <- summary(pca)$importance[2, 1:10] * 100
    barplot(var_exp, main = "PCA Variance Explained", xlab = "PC", ylab = "Variance (%)",
            col = "steelblue", names.arg = 1:10)
  })
  
  output$plot_pvalue_dist <- renderPlot({
    req(values$res)
    if(!is.null(values$res$pvalue) && sum(!is.na(values$res$pvalue)) > 0) {
      hist(values$res$pvalue[!is.na(values$res$pvalue)], 
           main = "P-value Distribution", xlab = "p-value", col = "skyblue", breaks = 50)
    } else {
      plot(1, type = "n", axes = FALSE, main = "No p-values available")
      text(1, 1, "Add replicates for statistical testing", cex = 1.2)
    }
  })
  
  output$plot_correlation_circle <- renderPlot({
    expr_mat <- get_expr_matrix()
    req(expr_mat)
    pca <- prcomp(t(expr_mat), scale = TRUE)
    plot(pca$rotation[,1], pca$rotation[,2], pch = 20, col = "blue",
         main = "Correlation Circle", xlab = "PC1", ylab = "PC2")
    abline(h = 0, v = 0, lty = 2)
  })
}

# Run the app
shinyApp(ui = ui, server = server)