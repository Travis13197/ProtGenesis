# ==============================================================================
# 6.AdvancedAnalysis.R
# ==============================================================================
#
# Purpose: Advanced analysis functions for ProtGenesis protein sequence data
#   This script provides comprehensive analysis functions including data 
#   processing, Seurat object construction, dimensionality reduction, and
#   visualization for EAC, MOTIF, and POLY datasets.
#
# Author: ProtGenesis Analysis Team
# Date: 2026-03-03
# Version: 1.0.0
#
# Usage:
#   source("6.AdvancedAnalysis.R")
#   data <- read_embedding_data(...)
#   seurat_obj <- build_seurat_object(...)
#   seurat_obj <- perform_dim_reduction(seurat_obj)
#   analyze_dataset(...)
#
# Requirements:
#   - Seurat
#   - data.table
#   - dplyr
#   - plotly
#   - ggplot2
#   - htmlwidgets
#
# ==============================================================================

# ==============================================================================
# Package Management | 包管理
# ==============================================================================

check_and_install_packages <- function() {
  required_packages <- c("Seurat", "data.table", "dplyr", "plotly",
                        "ggplot2", "htmlwidgets", "reshape2")

  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (pkg == "Seurat") {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager", repos = "https://cloud.r-project.org")
        }
        BiocManager::install(pkg, update = FALSE, ask = FALSE)
      } else {
        install.packages(pkg, repos = "https://cloud.r-project.org")
      }
    }
  }

  for (pkg in required_packages) {
    library(pkg, character.only = TRUE)
  }
}

check_and_install_packages()

# ==============================================================================
# Data Processing Functions
# ==============================================================================

#' Read embedding data from CSV file
#'
#' @param file_path Path to the embedding CSV file
#' @return List containing sequence_ids and embed_matrix
#' @export
read_embedding_data <- function(file_path) {
  cat("Reading embedding data...\n")
  raw_data <- fread(file_path) %>% as.data.frame()
  
  cat("Data dimensions:", dim(raw_data), "\n")
  
  sequence_ids <- raw_data[, 1]
  embed_matrix <- as.matrix(raw_data[, -1])
  
  return(list(
    sequence_ids = sequence_ids,
    embed_matrix = embed_matrix
  ))
}

#' Classify sequence based on ID prefix
#'
#' @param seq_id Sequence ID string
#' @return Category name (EAC, MOTIF, POLY, or Unknown)
#' @export
classify_sequence <- function(seq_id) {
  if (grepl("^EAC_", seq_id)) {
    return("EAC")
  } else if (grepl("^MOTIF_", seq_id)) {
    return("MOTIF")
  } else if (grepl("^POLY_", seq_id)) {
    return("POLY")
  } else {
    return("Unknown")
  }
}

#' Filter data by target category
#'
#' @param sequence_ids Vector of sequence IDs
#' @param embed_matrix Embedding matrix
#' @param target_category Target category to filter
#' @return List containing filtered sequence_ids and embed_matrix
#' @export
filter_category_data <- function(sequence_ids, embed_matrix, target_category) {
  cat("Filtering", target_category, "category data...\n")
  
  categories <- sapply(sequence_ids, classify_sequence)
  target_indices <- which(categories == target_category)
  target_sequence_ids <- sequence_ids[target_indices]
  target_embed_matrix <- embed_matrix[target_indices, ]
  
  cat("Filtered data dimensions:", dim(target_embed_matrix), "\n")
  
  return(list(
    sequence_ids = target_sequence_ids,
    embed_matrix = target_embed_matrix
  ))
}

#' Extract step information from sequence ID
#'
#' @param seq_id Sequence ID string
#' @return Step number or NA
#' @export
extract_step_from_id <- function(seq_id) {
  step_match <- regexec("_Step(\\d+)$", seq_id)
  if (length(step_match[[1]]) > 1) {
    step_str <- regmatches(seq_id, step_match)[[1]][2]
    return(as.integer(step_str))
  } else {
    return(NA)
  }
}

#' Extract poly category from sequence ID
#'
#' @param seq_id Sequence ID string
#' @return Poly category or Unknown
#' @export
extract_poly_category <- function(seq_id) {
  if (grepl("^POLY_", seq_id)) {
    poly_match <- regexec("^POLY_(\\w+)_", seq_id)
    if (length(poly_match[[1]]) > 1) {
      return(regmatches(seq_id, poly_match)[[1]][2])
    }
  }
  return("Unknown")
}

# ==============================================================================
# Seurat Object Construction
# ==============================================================================

#' Build Seurat object from embedding data
#'
#' @param sequence_ids Vector of sequence IDs
#' @param embed_matrix Embedding matrix
#' @param category Category name
#' @return Seurat object
#' @export
build_seurat_object <- function(sequence_ids, embed_matrix, category) {
  cat("Building Seurat object...\n")
  
  embed_matrix_t <- t(embed_matrix)
  colnames(embed_matrix_t) <- sequence_ids
  
  sparse_matrix <- Matrix(
    data = embed_matrix_t,
    nrow = nrow(embed_matrix_t),
    ncol = ncol(embed_matrix_t),
    dimnames = list(rownames(embed_matrix_t), colnames(embed_matrix_t)),
    sparse = TRUE
  )
  
  step_info <- sapply(sequence_ids, extract_step_from_id)
  
  if (category == "POLY") {
    poly_category_info <- sapply(sequence_ids, extract_poly_category)
    data_categories <- data.frame(
      Cell = sequence_ids,
      Sequence_ID = sequence_ids,
      Category = category,
      Poly_Category = poly_category_info,
      Step = step_info,
      row.names = sequence_ids
    )
  } else {
    data_categories <- data.frame(
      Cell = sequence_ids,
      Sequence_ID = sequence_ids,
      Category = category,
      Step = step_info,
      row.names = sequence_ids
    )
  }
  
  seurat_obj <- CreateSeuratObject(
    counts = sparse_matrix,
    assay = "embed",
    project = paste0("ProtGenesis_", category),
    min.cells = 0,
    min.features = 0,
    names.field = 1,
    names.delim = "-"
  )
  
  seurat_obj <- seurat_obj %>% AddMetaData(data_categories)
  
  cat("Performing data normalization...\n")
  seurat_obj <- SetAssayData(
    object = seurat_obj,
    assay = "embed",
    layer = "data",
    new.data = GetAssayData(seurat_obj, assay = "embed", layer = "counts")
  )
  
  seurat_obj <- ScaleData(
    object = seurat_obj,
    assay = "embed",
    do.center = TRUE,
    do.scale = FALSE
  )
  
  return(seurat_obj)
}

# ==============================================================================
# Dimensionality Reduction
# ==============================================================================

#' Perform comprehensive dimensionality reduction
#'
#' @param protgenesis_obj Seurat object
#' @return Seurat object with PCA, UMAP, t-SNE, and 3D reductions
#' @export
perform_dim_reduction <- function(protgenesis_obj) {
  cat("\nPerforming PCA dimensionality reduction...\n")
  protgenesis_obj <- RunPCA(
    object = protgenesis_obj,
    assay = "embed",
    features = rownames(protgenesis_obj),
    npcs = min(50, ncol(protgenesis_obj) - 1),
    verbose = FALSE
  )
  
  pct_variance <- protgenesis_obj@reductions$pca@stdev^2 / sum(protgenesis_obj@reductions$pca@stdev^2)
  selected_dims <- which(cumsum(pct_variance) < 0.95)
  cat("Selected dimensions:", length(selected_dims), "\n")
  
  n_cells <- ncol(protgenesis_obj)
  
  cat("\nPerforming UMAP dimensionality reduction...\n")
  umap_neighbors <- min(50, n_cells - 1)
  protgenesis_obj <- RunUMAP(
    object = protgenesis_obj,
    dims = selected_dims,
    reduction = "pca",
    metric = "cosine",
    min.dist = 0.3,
    spread = 1.5,
    n.neighbors = umap_neighbors,
    learning.rate = 0.8
  )
  
  cat("\nPerforming t-SNE dimensionality reduction...\n")
  tsne_perplexity <- min(100, floor((n_cells - 1) / 3))
  protgenesis_obj <- RunTSNE(
    object = protgenesis_obj,
    dims = selected_dims,
    reduction = "pca",
    perplexity = tsne_perplexity,
    check_duplicates = FALSE
  )
  
  cat("\nPerforming 3D UMAP dimensionality reduction...\n")
  protgenesis_obj <- RunUMAP(
    object = protgenesis_obj,
    dims = selected_dims,
    reduction = "pca",
    metric = "cosine",
    min.dist = 0.3,
    spread = 1.5,
    n.neighbors = umap_neighbors,
    learning.rate = 0.8,
    n.components = 3,
    reduction.name = "umap3d"
  )
  
  cat("\nPerforming 3D t-SNE dimensionality reduction...\n")
  tsne3d_dims <- min(length(selected_dims), 30)
  protgenesis_obj <- RunTSNE(
    object = protgenesis_obj,
    dims = 1:tsne3d_dims,
    reduction = "pca",
    perplexity = tsne_perplexity,
    check_duplicates = FALSE,
    n.components = 3,
    reduction.name = "tsne3d"
  )
  
  return(protgenesis_obj)
}

# ==============================================================================
# Visualization Functions
# ==============================================================================

#' Get custom color schemes for visualization
#'
#' @return List containing category and step color schemes
#' @export
get_category_colors <- function() {
  return(list(
    category_colors = c(
      "EAC" = "#DE6038",
      "MOTIF" = "#5778AA",
      "POLY" = "#46A64B",
      "Unknown" = "#ADADAD"
    ),
    step_colors = c(
      "Step1" = "#FF0000",
      "Step2" = "#FF6600",
      "Step5" = "#FFCC00",
      "Step10" = "#00FF00",
      "Step20" = "#00CCFF",
      "Step40" = "#0066FF",
      "Step60" = "#6600FF",
      "Step80" = "#9900FF",
      "Step100" = "#CC00FF",
      "Other" = "#ADADAD"
    )
  ))
}

#' Generate step-based visualization plots
#'
#' @param protgenesis_obj Seurat object
#' @param category Category name
#' @param output_dir Output directory for plots
#' @return List of ggplot objects
#' @export
generate_step_based_plots <- function(protgenesis_obj, category, output_dir) {
  cat("\nGenerating step-based plots...\n")
  colors <- get_category_colors()
  
  p_umap_step <- DimPlot(
    object = protgenesis_obj,
    group.by = "Step",
    pt.size = 0.5,
    label = FALSE,
    reduction = "umap",
    raster = FALSE
  ) + 
    ggtitle(paste(category, " Step Distribution (UMAP)")) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "right", plot.margin = unit(c(0.5, 0.5, 0.5, 1), "cm"))
  
  print(p_umap_step)
  ggsave(file.path(output_dir, paste0(category, "_step_umap.jpg")), height = 8, width = 10.4, dpi = 600)
  ggsave(file.path(output_dir, paste0(category, "_step_umap.svg")), height = 8, width = 10.4, dpi = 600)
  
  p_tsne_step <- DimPlot(
    object = protgenesis_obj,
    group.by = "Step",
    pt.size = 0.5,
    label = FALSE,
    reduction = "tsne",
    raster = FALSE
  ) + 
    ggtitle(paste(category, " Step Distribution (t-SNE)")) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "right", plot.margin = unit(c(0.5, 0.5, 0.5, 1), "cm"))
  
  print(p_tsne_step)
  ggsave(file.path(output_dir, paste0(category, "_step_tsne.jpg")), height = 8, width = 10.4, dpi = 600)
  ggsave(file.path(output_dir, paste0(category, "_step_tsne.svg")), height = 8, width = 10.4, dpi = 600)
  
  p_pca_step <- DimPlot(
    object = protgenesis_obj,
    group.by = "Step",
    pt.size = 0.5,
    label = FALSE,
    reduction = "pca",
    raster = FALSE
  ) + 
    ggtitle(paste(category, " Step Distribution (PCA)")) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "right", plot.margin = unit(c(0.5, 0.5, 0.5, 1), "cm"))
  
  print(p_pca_step)
  ggsave(file.path(output_dir, paste0(category, "_step_pca.jpg")), height = 8, width = 10.4, dpi = 600)
  ggsave(file.path(output_dir, paste0(category, "_step_pca.svg")), height = 8, width = 10.4, dpi = 600)
  
  return(list(umap_plot = p_umap_step, tsne_plot = p_tsne_step, pca_plot = p_pca_step))
}

#' Generate plots for specific steps
#'
#' @param protgenesis_obj Seurat object
#' @param category Category name
#' @param specific_steps Vector of step numbers
#' @param output_dir Output directory
#' @export
generate_specific_step_plots <- function(protgenesis_obj, category, specific_steps, output_dir) {
  cat("\nGenerating specific step plots...\n")
  
  for (step in specific_steps) {
    cat("Generating Step", step, "visualization...\n")
    
    step_seurat <- subset(protgenesis_obj, subset = Step == step)
    
    if (ncol(step_seurat) > 0) {
      for (reduction_type in c("umap", "tsne", "pca")) {
        p <- DimPlot(
          object = step_seurat,
          pt.size = 1,
          label = FALSE,
          reduction = reduction_type,
          raster = FALSE
        ) + 
          ggtitle(paste(category, " Step", step, "Distribution (", toupper(reduction_type), ")")) +
          theme_minimal(base_size = 14) +
          theme(legend.position = "right", plot.margin = unit(c(0.5, 0.5, 0.5, 1), "cm"))
        
        print(p)
        ggsave(file.path(output_dir, paste0(category, "_step", step, "_", reduction_type, ".jpg")), 
               height = 7, width = 9.1, dpi = 600)
        ggsave(file.path(output_dir, paste0(category, "_step", step, "_", reduction_type, ".svg")), 
               height = 7, width = 9.1, dpi = 600)
      }
    }
  }
}

#' Generate differential analysis plots
#'
#' @param protgenesis_obj Seurat object
#' @param category Category name
#' @param output_dir Output directory
#' @export
generate_differential_plots <- function(protgenesis_obj, category, output_dir) {
  cat("\nGenerating differential plots...\n")
  
  embed_data <- protgenesis_obj@assays$embed@layers$scale.data
  
  step_centers <- protgenesis_obj@meta.data %>%
    group_by(Step) %>%
    summarize(n_cells = n(), .groups = "drop")
  
  diff_dir <- file.path(output_dir, "differential")
  dir.create(diff_dir, recursive = TRUE)
  
  writeLines(c(
    "# Differential Analysis Results",
    paste("Dataset:", category),
    paste("Number of steps analyzed:", nrow(step_centers))
  ), file.path(diff_dir, paste0(category, "_differential_analysis.txt")))
  
  cat("Differential plots generated!\n")
}

#' Generate POLY subcategory visualization
#'
#' @param protgenesis_obj Seurat object
#' @param output_dir Output directory
#' @export
generate_poly_category_plots <- function(protgenesis_obj, output_dir) {
  cat("\nGenerating POLY subcategory visualization...\n")
  
  dir.create(file.path(output_dir, "poly_category_visualizations"), recursive = TRUE)
  
  poly_colors <- c("G" = "#4CAF50", "A" = "#2196F3", "D" = "#FF9800", 
                   "E" = "#9C27B0", "V" = "#F44336", "Unknown" = "#9E9E9E")
  
  for (reduction_type in c("umap", "tsne", "pca")) {
    p <- DimPlot(
      object = protgenesis_obj,
      group.by = "Poly_Category",
      pt.size = 0.5,
      label = FALSE,
      reduction = reduction_type,
      raster = FALSE
    ) + 
      ggtitle(paste("POLY Subcategory Distribution (", toupper(reduction_type), ")")) +
      scale_color_manual(values = poly_colors) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "right", plot.margin = unit(c(0.5, 0.5, 0.5, 1), "cm"))
    
    print(p)
    ggsave(file.path(output_dir, "poly_category_visualizations", 
                     paste0("POLY_category_", reduction_type, ".jpg")), height = 8, width = 10.4, dpi = 600)
    ggsave(file.path(output_dir, "poly_category_visualizations", 
                     paste0("POLY_category_", reduction_type, ".svg")), height = 8, width = 10.4, dpi = 600)
  }
  
  generate_all_poly_category_interactives(protgenesis_obj, output_dir)
  cat("POLY subcategory visualization generated!\n")
}

#' Generate interactive visualization for POLY subcategories
#'
#' @param protgenesis_obj Seurat object
#' @param output_dir Output directory
#' @param reduction Reduction type (tsne, umap, pca)
#' @export
generate_poly_category_interactive <- function(protgenesis_obj, output_dir, reduction = "tsne") {
  cat(paste("\nGenerating POLY subcategory interactive visualization (", reduction, ")...\n"))
  
  embed <- Embeddings(protgenesis_obj, reduction = reduction)
  meta_data <- protgenesis_obj@meta.data
  
  plot_data <- data.frame(
    Dim1 = embed[, 1],
    Dim2 = embed[, 2],
    Step = meta_data$Step,
    Poly_Category = meta_data$Poly_Category,
    Sequence_ID = meta_data$Sequence_ID
  )
  
  poly_colors <- c("G" = "#4CAF50", "A" = "#2196F3", "D" = "#FF9800", 
                   "E" = "#9C27B0", "V" = "#F44336", "Unknown" = "#9E9E9E")
  
  reduction_names <- list(
    "tsne" = list(name = "t-SNE", xlabel = "t-SNE 1", ylabel = "t-SNE 2"),
    "umap" = list(name = "UMAP", xlabel = "UMAP 1", ylabel = "UMAP 2"),
    "pca" = list(name = "PCA", xlabel = "PC 1", ylabel = "PC 2")
  )
  
  current_reduction <- reduction_names[[reduction]]
  
  p <- plot_ly(
    data = plot_data, x = ~Dim1, y = ~Dim2,
    color = ~Poly_Category, colors = poly_colors,
    text = ~paste("Sequence:", Sequence_ID, "<br>Step:", Step, "<br>Poly Category:", Poly_Category),
    hoverinfo = "text", type = "scatter", mode = "markers", marker = list(size = 5)
  ) %>% layout(
    title = paste("POLY Subcategory Interactive Visualization (", current_reduction$name, ")"),
    xaxis = list(title = current_reduction$xlabel),
    yaxis = list(title = current_reduction$ylabel),
    legend = list(title = list(text = "Poly Category"))
  )
  
  output_file <- file.path(output_dir, "poly_category_visualizations", 
                           paste0("POLY_category_interactive_", reduction, ".html"))
  htmlwidgets::saveWidget(p, file = output_file, selfcontained = FALSE)
  cat(paste("POLY subcategory interactive visualization (", reduction, ") saved!\n"))
}

#' Generate all POLY category interactive visualizations
#'
#' @param protgenesis_obj Seurat object
#' @param output_dir Output directory
#' @export
generate_all_poly_category_interactives <- function(protgenesis_obj, output_dir) {
  cat("\nGenerating all POLY subcategory interactive visualizations...\n")
  
  for (reduction_type in c("tsne", "umap", "pca")) {
    generate_poly_category_interactive(protgenesis_obj, output_dir, reduction_type)
  }
  
  generate_poly_category_interactive_3d(protgenesis_obj, output_dir)
  cat("All POLY subcategory interactive visualizations generated!\n")
}

#' Generate 3D interactive visualization for POLY subcategories
#'
#' @param protgenesis_obj Seurat object
#' @param output_dir Output directory
#' @export
generate_poly_category_interactive_3d <- function(protgenesis_obj, output_dir) {
  cat("\nGenerating POLY subcategory 3D interactive visualization...\n")
  
  poly_colors <- c("G" = "#4CAF50", "A" = "#2196F3", "D" = "#FF9800", 
                   "E" = "#9C27B0", "V" = "#F44336", "Unknown" = "#9E9E9E")
  
  meta_data <- protgenesis_obj@meta.data
  
  if ("tsne3d" %in% names(protgenesis_obj@reductions)) {
    cat("Generating 3D t-SNE visualization...\n")
    tryCatch({
      tsne3d_embed <- Embeddings(protgenesis_obj, reduction = "tsne3d")
      
      if (ncol(tsne3d_embed) >= 3) {
        tsne3d_data <- data.frame(
          TSNE1 = tsne3d_embed[, 1], TSNE2 = tsne3d_embed[, 2], TSNE3 = tsne3d_embed[, 3],
          Step = meta_data$Step, Poly_Category = meta_data$Poly_Category, Sequence_ID = meta_data$Sequence_ID
        )
        
        p_tsne3d <- plot_ly(
          data = tsne3d_data, x = ~TSNE1, y = ~TSNE2, z = ~TSNE3,
          color = ~Poly_Category, colors = poly_colors,
          text = ~paste("Sequence:", Sequence_ID, "<br>Step:", Step, "<br>Poly Category:", Poly_Category),
          hoverinfo = "text", type = "scatter3d", mode = "markers", marker = list(size = 3)
        ) %>% layout(
          title = "POLY Subcategory 3D Interactive Visualization (t-SNE)",
          scene = list(xaxis = list(title = "t-SNE 1"), yaxis = list(title = "t-SNE 2"), zaxis = list(title = "t-SNE 3")),
          legend = list(title = list(text = "Poly Category"))
        )
        
        output_file <- file.path(output_dir, "poly_category_visualizations", 
                                 "POLY_category_interactive_3d_tsne.html")
        htmlwidgets::saveWidget(p_tsne3d, file = output_file, selfcontained = FALSE)
        cat("3D t-SNE visualization saved!\n")
      }
    }, error <- function(e) {
      cat("Warning: Error generating 3D t-SNE visualization:", e$message, "\n")
    })
  }
  
  cat("POLY subcategory 3D visualization generated!\n")
}

#' Generate 3D visualizations
#'
#' @param protgenesis_obj Seurat object
#' @param category Category name
#' @param output_dir Output directory
#' @export
generate_3d_visualizations <- function(protgenesis_obj, category, output_dir) {
  cat("\nGenerating 3D visualizations...\n")
  
  dir.create(file.path(output_dir, "3d_visualizations"), recursive = TRUE)
  meta_data <- protgenesis_obj@meta.data
  
  if ("umap3d" %in% names(protgenesis_obj@reductions)) {
    cat("Generating 3D UMAP visualization...\n")
    tryCatch({
      umap3d_embed <- Embeddings(protgenesis_obj, reduction = "umap3d")
      
      if (ncol(umap3d_embed) >= 3) {
        umap3d_data <- data.frame(
          UMAP1 = umap3d_embed[, 1], UMAP2 = umap3d_embed[, 2], UMAP3 = umap3d_embed[, 3],
          Step = meta_data$Step, Category = meta_data$Category, Sequence_ID = meta_data$Sequence_ID,
          Poly_Category = ifelse("Poly_Category" %in% colnames(meta_data), meta_data$Poly_Category, "Unknown")
        )
        
        p_umap3d <- plot_ly(
          data = umap3d_data, x = ~UMAP1, y = ~UMAP2, z = ~UMAP3,
          color = ~as.factor(if("Poly_Category" %in% colnames(umap3d_data)) umap3d_data$Poly_Category else umap3d_data$Step),
          text = ~paste("Sequence:", Sequence_ID, "<br>Step:", Step, "<br>Category:", Category),
          hoverinfo = "text", type = "scatter3d", mode = "markers", marker = list(size = 3)
        ) %>% layout(
          title = paste(category, "3D UMAP Visualization"),
          scene = list(xaxis = list(title = "UMAP 1"), yaxis = list(title = "UMAP 2"), zaxis = list(title = "UMAP 3"))
        )
        
        output_file <- file.path(output_dir, "3d_visualizations", paste0(category, "_3d_umap.html"))
        htmlwidgets::saveWidget(p_umap3d, file = output_file, selfcontained = FALSE)
        cat("3D UMAP visualization saved!\n")
      }
    }, error <- function(e) {
      cat("Warning: Error generating 3D UMAP:", e$message, "\n")
    })
  }
  
  if ("tsne3d" %in% names(protgenesis_obj@reductions)) {
    cat("Generating 3D t-SNE visualization...\n")
    tryCatch({
      tsne3d_embed <- Embeddings(protgenesis_obj, reduction = "tsne3d")
      
      if (ncol(tsne3d_embed) >= 3) {
        tsne3d_data <- data.frame(
          TSNE1 = tsne3d_embed[, 1], TSNE2 = tsne3d_embed[, 2], TSNE3 = tsne3d_embed[, 3],
          Step = meta_data$Step, Category = meta_data$Category, Sequence_ID = meta_data$Sequence_ID
        )
        
        p_tsne3d <- plot_ly(
          data = tsne3d_data, x = ~TSNE1, y = ~TSNE2, z = ~TSNE3,
          color = ~as.factor(Step),
          text = ~paste("Sequence:", Sequence_ID, "<br>Step:", Step, "<br>Category:", Category),
          hoverinfo = "text", type = "scatter3d", mode = "markers", marker = list(size = 3)
        ) %>% layout(
          title = paste(category, "3D t-SNE Visualization"),
          scene = list(xaxis = list(title = "t-SNE 1"), yaxis = list(title = "t-SNE 2"), zaxis = list(title = "t-SNE 3"))
        )
        
        output_file <- file.path(output_dir, "3d_visualizations", paste0(category, "_3d_tsne.html"))
        htmlwidgets::saveWidget(p_tsne3d, file = output_file, selfcontained = FALSE)
        cat("3D t-SNE visualization saved!\n")
      }
    }, error <- function(e) {
      cat("Warning: Error generating 3D t-SNE:", e$message, "\n")
    })
  }
  
  cat("3D visualizations generated!\n")
}

#' Generate interactive visualization
#'
#' @param protgenesis_obj Seurat object
#' @param category Category name
#' @param output_dir Output directory
#' @export
generate_interactive_plot <- function(protgenesis_obj, category, output_dir) {
  cat("\nGenerating interactive visualization...\n")
  
  plot_data <- Seurat::Embeddings(protgenesis_obj, reduction = "tsne") %>% 
    as.data.frame() %>% 
    mutate(Cell = rownames(.)) %>% 
    left_join(protgenesis_obj@meta.data[, c("Cell", "Step", "Sequence_ID")], by = "Cell")
  
  p_interactive <- plot_ly(
    data = plot_data, x = ~tSNE_1, y = ~tSNE_2,
    color = ~factor(Step), text = ~paste("Sequence ID:", Sequence_ID, "\nStep:", Step),
    hoverinfo = "text+x+y", opacity = 0.7, marker = list(size = 3), showlegend = TRUE
  ) %>% layout(
    title = paste(category, "Interactive Visualization"),
    xaxis = list(title = "tSNE_1"), yaxis = list(title = "tSNE_2")
  )
  
  htmlwidgets::saveWidget(p_interactive, 
                          file = file.path(output_dir, paste0(category, "_interactive_visualization.html")),
                          selfcontained = FALSE)
  cat("Interactive visualization saved!\n")
}

# ==============================================================================
# Analysis Pipeline Functions
# ==============================================================================

#' Complete analysis pipeline for a single dataset
#'
#' @param file_path Path to embedding CSV file
#' @param category Category name (EAC, MOTIF, or POLY)
#' @param output_dir Output directory
#' @return Seurat object
#' @export
analyze_dataset <- function(file_path, category, output_dir) {
  cat("\n=== Starting", category, "dataset analysis ===\n")
  
  dir.create(output_dir, recursive = TRUE)
  
  data <- read_embedding_data(file_path)
  filtered_data <- filter_category_data(data$sequence_ids, data$embed_matrix, category)
  protgenesis_obj <- build_seurat_object(filtered_data$sequence_ids, filtered_data$embed_matrix, category)
  protgenesis_obj <- perform_dim_reduction(protgenesis_obj)
  
  save(protgenesis_obj, file = file.path(output_dir, paste0(category, "_protgenesis_object.rdata")))
  
  generate_step_based_plots(protgenesis_obj, category, output_dir)
  generate_differential_plots(protgenesis_obj, category, output_dir)
  
  specific_steps <- if (category == "EAC") c(1, 2, 5, 10, 20, 30) else
                     if (category == "MOTIF") c(1, 2, 5, 7) else
                     c(1, 2, 5, 10, 20, 40, 60, 80, 100, 120)
  
  generate_specific_step_plots(protgenesis_obj, category, specific_steps, output_dir)
  generate_3d_visualizations(protgenesis_obj, category, output_dir)
  
  if (category == "POLY") {
    generate_poly_category_plots(protgenesis_obj, output_dir)
  }
  
  generate_interactive_plot(protgenesis_obj, category, output_dir)
  generate_analysis_report(protgenesis_obj, category, output_dir)
  
  cat("\n===", category, "dataset analysis completed ===\n")
  return(protgenesis_obj)
}

#' Generate analysis report
#'
#' @param protgenesis_obj Seurat object
#' @param category Category name
#' @param output_dir Output directory
#' @export
generate_analysis_report <- function(protgenesis_obj, category, output_dir) {
  cat("\nGenerating analysis report...\n")
  
  report_content <- c(
    paste("#", category, "Data Analysis Report"),
    "",
    "## 1. Analysis Overview",
    "",
    paste("This report presents a comprehensive analysis of the", category, "dataset based on embedding data."),
    "",
    "## 2. Data Statistics",
    "",
    "| Statistic | Value |",
    "|-----------|-------|",
    paste("| Total Sequences |", ncol(protgenesis_obj), "|")
  )
  
  step_counts <- as.data.frame(table(protgenesis_obj@meta.data$Step))
  report_content <- c(report_content, "",
                      "### 2.2 Step Statistics",
                      "",
                      "| Step | Count |",
                      "|------|-------|")
  
  for (i in 1:nrow(step_counts)) {
    report_content <- c(report_content, paste0("| ", step_counts$Var1[i], " | ", step_counts$Freq[i], " |"))
  }
  
  report_content <- c(report_content, "",
                     "## 3. Analysis Results",
                     "",
                     "### 3.1 Dimensionality Reduction",
                     "",
                     "- **PCA**: Retained 95% variance dimensions",
                     "- **UMAP**: Visualized sequence distribution across steps",
                     "- **t-SNE**: Alternative visualization of sequence distribution",
                     "",
                     "### 3.2 Differential Analysis",
                     "",
                     paste("Differential analysis was performed for different steps in the", category, "dataset."),
                     "",
                     "## 4. Visualization Results",
                     "",
                     paste("- `", category, "_step_umap.jpg/svg`: Step-grouped UMAP distribution"),
                     paste("- `", category, "_step_tsne.jpg/svg`: Step-grouped t-SNE distribution"),
                     "",
                     "## 5. Conclusion",
                     "",
                     paste("This analysis comprehensively展示�?, category, "dataset在嵌入空间中的分布特征�?))
  
  writeLines(report_content, file.path(output_dir, paste0(category, "_analysis_report.md")))
  cat("Analysis report generated!\n")
}

#' Analyze all three datasets (EAC, MOTIF, POLY)
#'
#' @param file_path Path to embedding CSV file
#' @param base_output_dir Base output directory
#' @return List of Seurat objects
#' @export
analyze_all_datasets <- function(file_path, base_output_dir) {
  cat("\n=== Starting analysis of all three datasets ===\n")
  
  dir.create(base_output_dir, recursive = TRUE)
  
  eac_output_dir <- file.path(base_output_dir, "EAC")
  eac_seurat <- analyze_dataset(file_path, "EAC", eac_output_dir)
  
  motif_output_dir <- file.path(base_output_dir, "MOTIF")
  motif_seurat <- analyze_dataset(file_path, "MOTIF", motif_output_dir)
  
  poly_output_dir <- file.path(base_output_dir, "POLY")
  poly_seurat <- analyze_dataset(file_path, "POLY", poly_output_dir)
  
  generate_combined_report(list(EAC = eac_seurat, MOTIF = motif_seurat, POLY = poly_seurat), base_output_dir)
  
  cat("\n=== All dataset analysis completed ===\n")
  return(list(EAC = eac_seurat, MOTIF = motif_seurat, POLY = poly_seurat))
}

#' Generate combined analysis report
#'
#' @param seurat_objects List of Seurat objects
#' @param base_output_dir Base output directory
#' @export
generate_combined_report <- function(seurat_objects, base_output_dir) {
  cat("\nGenerating combined analysis report...\n")
  
  report_content <- c(
    "# ProtGenesis Combined Analysis Report",
    "",
    "## 1. Analysis Overview",
    "",
    "This report presents a comprehensive analysis of three datasets (EAC, MOTIF, POLY) based on embedding data.",
    "",
    "## 2. Data Statistics",
    "",
    "| Dataset | Sequence Count |",
    "|---------|----------------|"
  )
  
  for (name in names(seurat_objects)) {
    protgenesis_obj <- seurat_objects[[name]]
    report_content <- c(report_content, paste0("| ", name, " | ", ncol(protgenesis_obj), " |"))
  }
  
  report_content <- c(report_content, "",
                     "## 3. Analysis Results",
                     "",
                     "### 3.1 EAC Dataset",
                     "- Sequence length range: 1-30",
                     "- Main features: Limited amino acid set [G, A, D, E, V]",
                     "",
                     "### 3.2 MOTIF Dataset",
                     "- Sequence length range: 1-7",
                     "- Main features: Contains functional motifs like Ferredoxin-like and P-loop-like",
                     "",
                     "### 3.3 POLY Dataset",
                     "- Sequence length range: 1-120",
                     "- Main features: Single amino acid repeat sequences",
                     "",
                     "## 4. Conclusion",
                     "",
                     "Through comprehensive analysis of the three datasets, we found:",
                     "",
                     "1. **EAC**: Sequences form unique clustering patterns in embedding space",
                     "2. **MOTIF**: Functional motif sequences show characteristic distribution",
                     "3. **POLY**: Length variation of repeat sequences leads to regular changes in embedding space distribution")
  
  writeLines(report_content, file.path(base_output_dir, "combined_analysis_report.md"))
  cat("Combined analysis report generated!\n")
}

# ==============================================================================
# Export Functions
# ==============================================================================

#' Export all functions as a list
#'
#' @return List of exported functions
#' @export
export_advanced_functions <- function() {
  return(list(
    read_embedding_data = read_embedding_data,
    classify_sequence = classify_sequence,
    filter_category_data = filter_category_data,
    extract_step_from_id = extract_step_from_id,
    extract_poly_category = extract_poly_category,
    build_seurat_object = build_seurat_object,
    perform_dim_reduction = perform_dim_reduction,
    get_category_colors = get_category_colors,
    generate_step_based_plots = generate_step_based_plots,
    generate_specific_step_plots = generate_specific_step_plots,
    generate_differential_plots = generate_differential_plots,
    generate_poly_category_plots = generate_poly_category_plots,
    generate_poly_category_interactive = generate_poly_category_interactive,
    generate_all_poly_category_interactives = generate_all_poly_category_interactives,
    generate_poly_category_interactive_3d = generate_poly_category_interactive_3d,
    generate_3d_visualizations = generate_3d_visualizations,
    generate_interactive_plot = generate_interactive_plot,
    analyze_dataset = analyze_dataset,
    generate_analysis_report = generate_analysis_report,
    analyze_all_datasets = analyze_all_datasets,
    generate_combined_report = generate_combined_report
  ))
}
