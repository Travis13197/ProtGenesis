# ==============================================================================
# 3.Visualization.R
# ==============================================================================
#
# Purpose: Data visualization functions for protein mutation analysis
#   This script provides functions for visualizing protein mutation positions
#   in various embedding spaces (tSNE, UMAP, PCA).
#
# Author: ProtGenesis Analysis Team
# Date: 2026-02-24
# Version: 1.0.0
#
# Usage:
#   source("3.Visualization.R")
#   p <- plot_mutation_positions(seurat_obj, positions)
#
# Requirements:
#   - ggplot2 (for plotting)
#   - dplyr (for data manipulation)
#   - Seurat (for accessing embedding data)
#
# ==============================================================================

# ==============================================================================
# Package Management
# ==============================================================================

check_and_install_packages <- function() {
  required_packages <- c("ggplot2", "dplyr", "Seurat", "plotly", "htmlwidgets")
  
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (pkg %in% c("Seurat")) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager")
        }
        BiocManager::install(pkg)
      } else {
        install.packages(pkg)
      }
    }
  }
  
  library(ggplot2)
  library(dplyr)
  library(Seurat)
  library(plotly)
  library(htmlwidgets)
}

check_and_install_packages()

# ==============================================================================
# Core Functions
# ==============================================================================

#' Plot mutation positions in embedding space
#'
#' This function visualizes mutation positions in specified embedding spaces
#' (tSNE, UMAP, or PCA) with options to highlight ORI points and add arrows.
#'
#' @param seurat_obj Seurat object containing embedding data
#' @param positions Numeric vector of mutation positions to plot
#' @param colors Optional vector of colors for different positions
#' @param highlight_ori Logical indicating whether to highlight ORI points (default: TRUE)
#' @param add_arrows Logical indicating whether to add arrows between ORI points (default: TRUE)
#' @param embedding_type Type of embedding to use: "tsne", "umap", or "pca" (default: "tsne")
#'
#' @return ggplot object with mutation positions visualization
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage with tSNE embedding
#' p <- plot_mutation_positions(seurat_obj, c(10, 20, 30))
#' print(p)
#'
#' # Using UMAP embedding with custom colors
#' p <- plot_mutation_positions(
#'   seurat_obj, 
#'   c(5, 15, 25),
#'   colors = c("#FF0000", "#00FF00", "#0000FF"),
#'   embedding_type = "umap"
#' )
#' print(p)
#' }
plot_mutation_positions <- function(seurat_obj, positions, colors = NULL, highlight_ori = TRUE, add_arrows = TRUE, embedding_type = "tsne") {
  if(!is.numeric(positions) || length(positions) < 1) {
    stop("positions必须是至少包含一个数字的向量")
  }

  positions <- unique(positions)

  valid_embeddings <- c("tsne", "umap", "pca")
  if(!(embedding_type %in% valid_embeddings)) {
    stop("embedding_type必须是'tsne', 'umap'或'pca'之一")
  }

  if(is.null(colors)) {
    colors <- c("#2b8cbe", "#e41a1c", "#4daf4a", "#984ea3", "#ff7f00")
  }

  if(embedding_type == "tsne") {
    emb_coords <- as.data.frame(seurat_obj@reductions$tsne@cell.embeddings)
    x_col <- "tSNE_1"
    y_col <- "tSNE_2"
  } else if(embedding_type == "umap") {
    emb_coords <- as.data.frame(seurat_obj@reductions$umap@cell.embeddings)
    x_col <- "UMAP_1"
    y_col <- "UMAP_2"
  } else if(embedding_type == "pca") {
    emb_coords <- as.data.frame(seurat_obj@reductions$pca@cell.embeddings)
    x_col <- "PC_1"
    y_col <- "PC_2"
  }

  mut_data <- cbind(seurat_obj@meta.data, emb_coords)

  filtered_data <- mut_data %>% 
    filter(Mutation_Position %in% positions)

  position_color_map <- setNames(colors[(1:length(positions) - 1) %% length(colors) + 1], positions)

  filtered_data <- filtered_data %>% 
    mutate(
      Position_Color = position_color_map[as.character(Mutation_Position)],
      Position_Label = factor(Mutation_Position, levels = positions)
    )

  p <- ggplot(filtered_data, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(aes(color = Position_Label), size = 3, alpha = 0.8) +
    scale_color_manual(
      values = position_color_map,
      labels = function(x) paste("Position ", x),
      name = "Mutation Position"
    ) +
    ggtitle(paste("Mutation Positions", paste(positions, collapse = ", "), "- ", toupper(embedding_type), "分布")) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank())

  if(highlight_ori) {
    ori_data <- filtered_data %>% filter(grepl("ORI", Cell))
    if(nrow(ori_data) > 0) {
      p <- p +
        geom_point(data = ori_data, aes(x = .data[[x_col]], y = .data[[y_col]], fill = Position_Color),
                  color = "black", size = 5, shape = 21, stroke = 1.2) +
        geom_text(data = ori_data, aes(x = .data[[x_col]], y = .data[[y_col]], label = Mutation_Position),
                 vjust = 0.5, hjust = 0.5, color = "black", size = 3.5)

      if(add_arrows && nrow(ori_data) >= 2) {
        ori_data <- ori_data[order(match(ori_data$Mutation_Position, positions)), ]

        arrows_data <- data.frame(
          x = ori_data[[x_col]][-nrow(ori_data)],
          y = ori_data[[y_col]][-nrow(ori_data)],
          xend = ori_data[[x_col]][-1],
          yend = ori_data[[y_col]][-1]
        )

        p <- p +
          geom_segment(data = arrows_data, aes(x = x, y = y, xend = xend, yend = yend),
                      arrow = arrow(length = unit(0.5, "cm")), color = "black")
      }
    }
  }

  return(p)
}

#' Structural Map Visualization Function
#'
#' This function creates comprehensive visualizations for structural maps
#' in various embedding spaces (PCA, tSNE, UMAP) in 2D and 3D formats.
#'
#' @param protgenesis_obj Seurat/ProtGenesis object containing embedding data
#' @param group_by Character string specifying the metadata column to group by
#' @param cols Color vector for groups
#' @param pt_size Point size for plotting (default: 0.02)
#' @param order Vector specifying the order of groups
#' @param split_by Character string specifying the metadata column to split by
#' @param label Logical indicating whether to show labels (default: TRUE)
#' @param raster Logical indicating whether to use rasterization (default: TRUE)
#' @param title Plot title
#' @param output_prefix Prefix for output filenames (default: "plot")
#' @param output_dir Directory for output files (default: "Figures")
#' @param run_pca Logical indicating whether to generate PCA visualization (default: FALSE)
#' @param run_tsne Logical indicating whether to generate tSNE visualization (default: TRUE)
#' @param run_umap Logical indicating whether to generate UMAP visualization (default: TRUE)
#' @param run_pca3d Logical indicating whether to generate PCA 3D visualization (default: FALSE)
#' @param run_tsne3d Logical indicating whether to generate tSNE 3D visualization (default: FALSE)
#' @param run_umap3d Logical indicating whether to generate UMAP 3D visualization (default: FALSE)
#' @param output_svg Logical indicating whether to output SVG format (default: TRUE)
#' @param output_jpg Logical indicating whether to output JPG format (default: TRUE)
#' @param output_html Logical indicating whether to output HTML format (default: TRUE)
#' @param height Plot height in inches (default: 6)
#' @param width Plot width in inches (default: 7)
#' @param dpi Resolution for output (default: 600)
#' @param verbose Logical indicating whether to show progress messages (default: TRUE)
#'
#' @return List containing all generated ggplot objects
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage with tSNE and UMAP
#' plots <- StructuralMapVis(
#'   protgenesis_obj,
#'   group_by = "Round_Info",
#'   cols = c("Round1" = "#DE6038", "Round2" = "#D99829")
#' )
#' }
StructuralMapVis <- function(
  protgenesis_obj,
  group_by = NULL,
  cols = NULL,
  pt_size = 0.02,
  order = NULL,
  split_by = NULL,
  label = TRUE,
  raster = TRUE,
  title = NULL,
  output_prefix = "plot",
  output_dir = "Figures",
  run_pca = FALSE,
  run_tsne = TRUE,
  run_umap = TRUE,
  run_pca3d = FALSE,
  run_tsne3d = FALSE,
  run_umap3d = FALSE,
  output_svg = TRUE,
  output_jpg = TRUE,
  output_html = TRUE,
  height = 6,
  width = 7,
  dpi = 600,
  verbose = TRUE
) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  result_list <- list()
  
  check_reduction <- function(reduction_name) {
    if (!(reduction_name %in% names(protgenesis_obj@reductions))) {
      if (verbose) {
        warning(sprintf("%s reduction not available in ProtGenesis object, skipping...", 
                        toupper(reduction_name)))
      }
      return(FALSE)
    }
    return(TRUE)
  }
  
  create_2d_plot <- function(reduction_type, display_name) {
    if (!check_reduction(reduction_type)) {
      return(NULL)
    }
    
    plot_args <- list(
      object = protgenesis_obj,
      reduction = reduction_type,
      pt.size = pt_size,
      label = label,
      raster = raster
    )
    
    if (!is.null(group_by)) {
      plot_args$group.by <- group_by
    }
    
    if (!is.null(cols)) {
      plot_args$cols <- cols
    }
    
    if (!is.null(order)) {
      plot_args$order <- order
    }
    
    if (!is.null(split_by)) {
      plot_args$split.by <- split_by
    }
    
    p <- do.call(DimPlot, plot_args)
    
    if (!is.null(title)) {
      p <- p + ggtitle(title)
    } else {
      p <- p + ggtitle(sprintf("%s Visualization", display_name))
    }
    
    p <- p + theme_minimal(base_size = 14)
    
    return(p)
  }
  
  save_plot <- function(plot_obj, plot_name) {
    if (is.null(plot_obj)) {
      return(NULL)
    }
    
    if (verbose) {
      cat(sprintf("Saving %s plots...\n", plot_name))
    }
    
    if (output_jpg) {
      jpg_file <- file.path(output_dir, sprintf("%s_%s.jpg", output_prefix, plot_name))
      ggsave(jpg_file, plot = plot_obj, height = height, width = width, dpi = dpi)
      if (verbose) {
        cat(sprintf("  Saved: %s\n", jpg_file))
      }
    }
    
    if (output_svg) {
      svg_file <- file.path(output_dir, sprintf("%s_%s.svg", output_prefix, plot_name))
      ggsave(svg_file, plot = plot_obj, height = height, width = width, dpi = dpi)
      if (verbose) {
        cat(sprintf("  Saved: %s\n", svg_file))
      }
    }
    
    if (output_html) {
      html_file <- file.path(output_dir, sprintf("%s_%s.html", output_prefix, plot_name))
      plotly_obj <- ggplotly(plot_obj)
      htmlwidgets::saveWidget(plotly_obj, file = html_file, selfcontained = TRUE)
      if (verbose) {
        cat(sprintf("  Saved: %s\n", html_file))
      }
    }
    
    return(plot_obj)
  }
  
  if (verbose) {
    cat("==============================================================================\n")
    cat("StructuralMapVis - Generating Visualizations\n")
    cat("==============================================================================\n\n")
  }
  
  if (run_pca) {
    if (verbose) cat("Generating PCA visualization...\n")
    p_pca <- create_2d_plot("pca", "PCA")
    if (!is.null(p_pca)) {
      print(p_pca)
      result_list$pca <- save_plot(p_pca, "pca")
    }
  }
  
  if (run_tsne) {
    if (verbose) cat("Generating tSNE visualization...\n")
    p_tsne <- create_2d_plot("tsne", "t-SNE")
    if (!is.null(p_tsne)) {
      print(p_tsne)
      result_list$tsne <- save_plot(p_tsne, "tsne")
    }
  }
  
  if (run_umap) {
    if (verbose) cat("Generating UMAP visualization...\n")
    p_umap <- create_2d_plot("umap", "UMAP")
    if (!is.null(p_umap)) {
      print(p_umap)
      result_list$umap <- save_plot(p_umap, "umap")
    }
  }
  
  if (run_pca3d || run_tsne3d || run_umap3d) {
    if (verbose) {
      cat("\nNote: 3D visualization requires additional setup. ")
      cat("Please ensure 3D reductions are computed first.\n")
    }
  }
  
  if (verbose) {
    cat("\n==============================================================================\n")
    cat("Visualization complete!\n")
    cat("==============================================================================\n")
  }
  
  return(invisible(result_list))
}

# ==============================================================================
# Initialization
# ==============================================================================

cat("==============================================================================\n")
cat("Visualization Functions Loaded Successfully!\n")
cat("==============================================================================\n")
cat("\n")
cat("Available functions:\n")
cat("  - plot_mutation_positions()    : Visualize mutation positions in embedding space\n")
cat("  - StructuralMapVis()           : Comprehensive structural map visualization\n")
cat("\n")
cat("For help, see ?plot_mutation_positions or ?StructuralMapVis\n")
cat("==============================================================================\n")
cat("\n")
cat("Author: ProtGenesis Analysis Team\n")
cat("Date: 2026-02-24\n")
cat("Version: 1.0.0\n")
cat("By: Chuanyang Liu\n")
cat("Please cite our paper if you use this script in your research:\n")
cat("Chuanyang Liu, et al. (2026). Universal physical principles govern the deterministic genesis of protein structure.\n")
cat("BioRxiv.10.64898/2026.02.20.706798.\n")
cat("\n")
cat("For help, see ?plot_mutation_positions or ?StructuralMapVis\n")
cat("==============================================================================\n")
