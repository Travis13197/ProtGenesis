# ==============================================================================
# 3.Visualization.r | 可视化函数
# ==============================================================================
#
# Purpose: Data visualization functions for protein mutation analysis
#   This script provides functions for visualizing protein mutation positions
#   in various embedding spaces (tSNE, UMAP, PCA).
#   本脚本提供在各种嵌入空间（tSNE、UMAP、PCA）中可视化蛋白质突变位置的函数。
#
# Author: ProtGenesis Analysis Team
# Date: 2026-02-24
# Version: 1.0.0
#
# Usage: | 使用方法
#   source("3.Visualization.r")
#   p <- plot_mutation_positions(protgenesis_obj, positions)
#
# Requirements: | 依赖项
#   - ggplot2 (for plotting)
#   - dplyr (for data manipulation)
#   - Seurat (for accessing embedding data)
#
# ==============================================================================

# ==============================================================================
# Package Management | 包管理
# ==============================================================================

check_and_install_packages <- function() {
  required_packages <- c("ggplot2", "dplyr", "Seurat", "plotly", "htmlwidgets")

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
# Core Functions
# ==============================================================================

#' Plot mutation positions in embedding space
#'
#' This function visualizes mutation positions in specified embedding spaces
#' (tSNE, UMAP, or PCA) with options to highlight ORI points and add arrows.
#'
#' @param protgenesis_obj Seurat object containing embedding data
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
#' p = plot_mutation_positions(protgenesis_obj, c(10, 20, 30))
#' print(p)
#'
#' # Using UMAP embedding with custom colors
#' p = plot_mutation_positions(
#'   protgenesis_obj, 
#'   c(5, 15, 25),
#'   colors = c("#FF0000", "#00FF00", "#0000FF"),
#'   embedding_type = "umap"
#' )
#' print(p)
#' }
# Plot mutation positions function | 绘制突变位置函数
plot_mutation_positions <- function(protgenesis_obj, positions, colors = NULL, highlight_ori = TRUE, add_arrows = TRUE, embedding_type = "tsne") {
  if(!is.numeric(positions) || length(positions) < 1) {
    stop("positions必须是至少包含一个数字的向量")
  }

  positions = unique(positions)

  valid_embeddings = c("tsne", "umap", "pca")
  if(!(embedding_type %in% valid_embeddings)) {
    stop("embedding_type必须是'tsne', 'umap'或'pca'之一")
  }

  if(is.null(colors)) {
    colors = c("#2b8cbe", "#e41a1c", "#4daf4a", "#984ea3", "#ff7f00")
  }

  if(embedding_type == "tsne") {
    emb_coords = as.data.frame(protgenesis_obj@reductions$tsne@cell.embeddings)
    x_col = "tSNE_1"
    y_col = "tSNE_2"
  } else if(embedding_type == "umap") {
    emb_coords = as.data.frame(protgenesis_obj@reductions$umap@cell.embeddings)
    x_col = "UMAP_1"
    y_col = "UMAP_2"
  } else if(embedding_type == "pca") {
    emb_coords = as.data.frame(protgenesis_obj@reductions$pca@cell.embeddings)
    x_col = "PC_1"
    y_col = "PC_2"
  }

  mut_data = cbind(protgenesis_obj@meta.data, emb_coords)

  filtered_data = mut_data %>% 
    filter(Mutation_Position %in% positions)

  position_color_map = setNames(colors[(1:length(positions) - 1) %% length(colors) + 1], positions)

  filtered_data = filtered_data %>% 
    mutate(
      Position_Color = position_color_map[as.character(Mutation_Position)],
      Position_Label = factor(Mutation_Position, levels = positions)
    )

  p = ggplot(filtered_data, aes(x = .data[[x_col]], y = .data[[y_col]])) +
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
    ori_data = filtered_data %>% filter(grepl("ORI", Cell))
    if(nrow(ori_data) > 0) {
      p = p +
        geom_point(data = ori_data, aes(x = .data[[x_col]], y = .data[[y_col]], fill = Position_Color),
                  color = "black", size = 5, shape = 21, stroke = 1.2) +
        geom_text(data = ori_data, aes(x = .data[[x_col]], y = .data[[y_col]], label = Mutation_Position),
                 vjust = 0.5, hjust = 0.5, color = "black", size = 3.5)

      if(add_arrows && nrow(ori_data) >= 2) {
        ori_data = ori_data[order(match(ori_data$Mutation_Position, positions)), ]

        arrows_data = data.frame(
          x = ori_data[[x_col]][-nrow(ori_data)],
          y = ori_data[[y_col]][-nrow(ori_data)],
          xend = ori_data[[x_col]][-1],
          yend = ori_data[[y_col]][-1]
        )

        p = p +
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
#' plots = StructuralMapVis(
#'   protgenesis_obj,
#'   group_by = "Round_Info",
#'   cols = c("Round1" = "#DE6038", "Round2" = "#D99829")
#' )
#' }
# Structural map visualization function | 结构图可视化函数
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
  
  result_list = list()
  
  check_reduction = function(reduction_name) {
    if (!(reduction_name %in% names(protgenesis_obj@reductions))) {
      if (verbose) {
        warning(sprintf("%s reduction not available in ProtGenesis object, skipping...", 
                        toupper(reduction_name)))
      }
      return(FALSE)
    }
    return(TRUE)
  }
  
  create_2d_plot = function(reduction_type, display_name) {
    if (!check_reduction(reduction_type)) {
      return(NULL)
    }
    
    plot_args = list(
      object = protgenesis_obj,
      reduction = reduction_type,
      pt.size = pt_size,
      label = label,
      raster = raster
    )
    
    if (!is.null(group_by)) {
      plot_args$group.by = group_by
    }
    
    if (!is.null(cols)) {
      plot_args$cols = cols
    }
    
    if (!is.null(order)) {
      plot_args$order = order
    }
    
    if (!is.null(split_by)) {
      plot_args$split.by = split_by
    }
    
    p = do.call(DimPlot, plot_args)
    
    if (!is.null(title)) {
      p = p + ggtitle(title)
    } else {
      p = p + ggtitle(sprintf("%s Visualization", display_name))
    }
    
    p = p + theme_minimal(base_size = 14)
    
    return(p)
  }
  
  save_plot = function(plot_obj, plot_name) {
    if (is.null(plot_obj)) {
      return(NULL)
    }
    
    if (verbose) {
      cat(sprintf("Saving %s plots...\n", plot_name))
    }
    
    if (output_jpg) {
      jpg_file = file.path(output_dir, sprintf("%s_%s.jpg", output_prefix, plot_name))
      ggsave(jpg_file, plot = plot_obj, height = height, width = width, dpi = dpi)
      if (verbose) {
        cat(sprintf("  Saved: %s\n", jpg_file))
      }
    }
    
    if (output_svg) {
      svg_file = file.path(output_dir, sprintf("%s_%s.svg", output_prefix, plot_name))
      ggsave(svg_file, plot = plot_obj, height = height, width = width, dpi = dpi)
      if (verbose) {
        cat(sprintf("  Saved: %s\n", svg_file))
      }
    }
    
    if (output_html) {
      html_file = file.path(output_dir, sprintf("%s_%s.html", output_prefix, plot_name))
      plotly_obj = ggplotly(plot_obj)
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
    p_pca = create_2d_plot("pca", "PCA")
    if (!is.null(p_pca)) {
      print(p_pca)
      result_list$pca = save_plot(p_pca, "pca")
    }
  }
  
  if (run_tsne) {
    if (verbose) cat("Generating tSNE visualization...\n")
    p_tsne = create_2d_plot("tsne", "t-SNE")
    if (!is.null(p_tsne)) {
      print(p_tsne)
      result_list$tsne = save_plot(p_tsne, "tsne")
    }
  }
  
  if (run_umap) {
    if (verbose) cat("Generating UMAP visualization...\n")
    p_umap = create_2d_plot("umap", "UMAP")
    if (!is.null(p_umap)) {
      print(p_umap)
      result_list$umap = save_plot(p_umap, "umap")
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

#' Visualize differential heatmap
#'
#' This function visualizes differential embeddings using heatmap.
#'
#' @param differential_df Data frame containing differential embeddings
#' @param output_dir Directory for output files
#' @param output_prefix Prefix for output filenames
#' @param verbose Logical indicating whether to show progress messages
#'
#' @return ggplot object with heatmap
#' @export
# Differential heatmap visualization function | 差分热图可视化函数
VisualizeDifferentialHeatmap <- function(differential_df,
                                         output_dir,
                                         output_prefix = "differential_heatmap",
                                         verbose = TRUE) {
  if (verbose) {
    cat("绘制差异嵌入热图...\n")
  }
  
  library(ggplot2)
  library(reshape2)
  
  embedding_heatmap_data = differential_df[, c('Full_Sequence', 'Sequence', 'Terminal_Type', grep('^dim-', colnames(differential_df), value = TRUE))]
  embedding_heatmap_data = embedding_heatmap_data[order(embedding_heatmap_data$Terminal_Type, embedding_heatmap_data$Sequence, embedding_heatmap_data$Full_Sequence), ]
  
  long_embedding_data = reshape2::melt(embedding_heatmap_data, id.vars = c('Full_Sequence', 'Sequence', 'Terminal_Type'), variable.name = 'Dimension', value.name = 'Value')
  long_embedding_data$Dimension = factor(long_embedding_data$Dimension, levels = unique(long_embedding_data$Dimension))
  
  long_embedding_data$seq_length = nchar(long_embedding_data$Full_Sequence)
  long_embedding_data = long_embedding_data[order(long_embedding_data$seq_length, long_embedding_data$Full_Sequence), ]
  
  ordered_seq_levels = unique(long_embedding_data$Full_Sequence)
  long_embedding_data$Full_Sequence = factor(long_embedding_data$Full_Sequence, levels = ordered_seq_levels)
  
  p = ggplot(long_embedding_data, aes(x = Dimension, y = Full_Sequence, fill = Value)) +
    geom_tile(color = "white", size = 0.02) +
    scale_fill_gradient2(low = '#d60404', mid = '#f5f5f541', high = '#0091d4', midpoint = 0) +
    facet_grid(Sequence ~ Terminal_Type, scales = 'free_y', space = 'free_y') +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 6),
      strip.text.x = element_text(angle = 0, hjust = 0.5, face = "bold"),
      strip.text.y = element_text(angle = 0, hjust = 0, face = "bold"),
      legend.position = 'right',
      panel.spacing = unit(0.05, 'lines'),
      panel.border = element_rect(color = "white", fill = NA, size = 1)
    ) +
    labs(
      x = '嵌入维度 (dim-0 到 dim-1023)',
      y = '蛋白质序列',
      fill = '差异嵌入值',
      title = '不同终端类型的差异嵌入热图（按序列分组）'
    )
  
  print(p)
  
  output_file_jpg = file.path(output_dir, paste0(output_prefix, ".jpg"))
  output_file_svg = file.path(output_dir, paste0(output_prefix, ".svg"))
  ggsave(output_file_jpg, plot = p, width = 12, height = 8, dpi = 1200)
  ggsave(output_file_svg, plot = p, width = 12, height = 8, dpi = 600)
  
  if (verbose) {
    cat(sprintf("热图已保存: %s\n", output_file_jpg))
  }
  
  invisible(p)
}

#' Visualize position boxplots
#'
#' This function visualizes differential embeddings using boxplots for specific positions.
#'
#' @param differential_df Data frame containing differential embeddings
#' @param positions_to_visualize Positions to visualize (default: c(1, 2, 3))
#' @param output_dir Directory for output files
#' @param output_prefix Prefix for output filenames
#' @param verbose Logical indicating whether to show progress messages
#'
#' @return NULL
#' @export
# Position boxplots visualization function | 位置箱线图可视化函数
VisualizePositionBoxplots <- function(differential_df,
                                       positions_to_visualize = c(1, 2, 3),
                                       output_dir,
                                       output_prefix = "position_boxplot",
                                       verbose = TRUE) {
  if (verbose) {
    cat("绘制特定位点箱线图...\n")
  }
  
  library(ggplot2)
  library(reshape2)
  
  for (pos in positions_to_visualize) {
    differential_df[[paste0('Position_', pos)]] = substr(differential_df$Full_Sequence, pos, pos)
  }
  
  for (pos in positions_to_visualize) {
    pos_col = paste0('Position_', pos)
    position_data = differential_df[, c('Full_Sequence', 'Sequence', 'Terminal_Type', pos_col, grep('^dim-', colnames(differential_df), value = TRUE))]
    
    pos_long_data = reshape2::melt(position_data,
                                      id.vars = c('Full_Sequence', 'Sequence', 'Terminal_Type', pos_col),
                                      variable.name = 'Dimension',
                                      value.name = 'Value')
    
    p_box = ggplot(pos_long_data, aes(x = get(pos_col), y = Value, fill = get(pos_col))) +
      geom_boxplot() +
      facet_grid(Terminal_Type ~ Dimension, scales = 'free_y') +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = 'none'
      ) +
      labs(
        x = paste('位点', pos, '氨基酸'),
        y = '差异嵌入值',
        title = paste('特定位点', pos, '氨基酸的差异嵌入分布')
      )
    
    print(p_box)
    
    output_file_jpg = file.path(output_dir, paste0(output_prefix, "_position", pos, ".jpg"))
    output_file_svg = file.path(output_dir, paste0(output_prefix, "_position", pos, ".svg"))
    ggsave(output_file_jpg, plot = p_box, width = 15, height = 6, dpi = 600)
    ggsave(output_file_svg, plot = p_box, width = 15, height = 6, dpi = 600)
    
    if (verbose) {
      cat(sprintf("位点 %d 箱线图已保存\n", pos))
    }
  }
  
  invisible(NULL)
}

#' Visualize scatterplots
#'
#' This function visualizes differential embeddings using scatterplots.
#'
#' @param differential_df Data frame containing differential embeddings
#' @param target_sequences Target sequences to visualize
#' @param output_dir Directory for output files
#' @param output_prefix Prefix for output filenames
#' @param color_palette_terminal Color palette for terminal types
#' @param color_palette_sequences Color palette for sequences
#' @param verbose Logical indicating whether to show progress messages
#'
#' @return NULL
#' @export
# Scatterplots visualization function | 散点图可视化函数
VisualizeScatterplots <- function(differential_df,
                                   target_sequences,
                                   output_dir,
                                   output_prefix = "scatterplot",
                                   color_palette_terminal = c("Nter" = "#FF0000", "Cter" = "#0000FF"),
                                   color_palette_sequences = NULL,
                                   verbose = TRUE) {
  if (verbose) {
    cat("绘制散点图...\n")
  }
  
  library(ggplot2)
  
  if (!all(c("tSNE_1", "tSNE_2") %in% colnames(differential_df))) {
    stop("differential_df 中缺少 tSNE_1 或 tSNE_2 列，请先运行 CalculateTSNEForDifferential")
  }
  
  selected_data = differential_df[differential_df$Sequence %in% target_sequences, ]
  
  variant_patterns = paste0('^', target_sequences, '[A-Z]$')
  for (pattern in variant_patterns) {
    selected_data = rbind(selected_data, differential_df[grepl(pattern, differential_df$Sequence), ])
  }
  
  plot_df = selected_data[, c('Full_Sequence', 'Sequence', 'Terminal_Type', 'tSNE_1', 'tSNE_2')]
  
  duplicate_rows = duplicated(plot_df[, c('tSNE_1', 'tSNE_2', 'Terminal_Type')])
  if (any(duplicate_rows)) {
    plot_df = plot_df[!duplicate_rows, ]
  }
  
  all_colors = NULL
  if (!is.null(color_palette_sequences)) {
    variant_colors = setNames(sapply(color_palette_sequences, function(color) {
      substr(color, 1, 7)
    }), paste0(names(color_palette_sequences), 'X'))
    all_colors = c(color_palette_sequences, variant_colors)
    
    plot_df$Sequence_Type = ifelse(plot_df$Sequence %in% target_sequences,
                                      plot_df$Sequence,
                                      paste0(substr(plot_df$Sequence, 1, nchar(plot_df$Sequence)-1), 'X'))
  }
  
  plot_df$Terminal_Type = factor(plot_df$Terminal_Type, levels = c('Nter', 'Cter'))
  
  unique_sequences = unique(plot_df$Sequence)
  for (seq in unique_sequences) {
    seq_df = plot_df[plot_df$Sequence == seq, ]
    p_seq = ggplot(seq_df, aes(x = tSNE_1, y = tSNE_2, color = Terminal_Type, shape = Terminal_Type)) +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_manual(values = color_palette_terminal, name = '终端类型') +
      scale_shape_manual(values = c(16, 17), name = '终端类型') +
      ggtitle(paste('t-SNE可视化：序列', seq, '（按终端类型着色）')) +
      xlim(-20, 20) + ylim(-20, 20) +
      labs(x = 'tSNE_1', y = 'tSNE_2') +
      theme_minimal(base_size = 14) +
      theme(
        panel.grid = element_blank(),
        legend.position = 'right',
        plot.title = element_text(hjust = 0.5)
      )
    
    nter_count = sum(seq_df$Terminal_Type == 'Nter')
    cter_count = sum(seq_df$Terminal_Type == 'Cter')
    if (nter_count > 0 && cter_count > 0) {
      nter_mean_x = mean(seq_df$tSNE_1[seq_df$Terminal_Type == 'Nter'])
      nter_mean_y = mean(seq_df$tSNE_2[seq_df$Terminal_Type == 'Nter'])
      cter_mean_x = mean(seq_df$tSNE_1[seq_df$Terminal_Type == 'Cter'])
      cter_mean_y = mean(seq_df$tSNE_2[seq_df$Terminal_Type == 'Cter'])
      title = paste0('t-SNE可视化：序列 ', seq, '（按终端类型着色）\n',
                      'Nter均值: (', round(nter_mean_x, 2), ', ', round(nter_mean_y, 2), '), ',
                      'Cter均值: (', round(cter_mean_x, 2), ', ', round(cter_mean_y, 2), ')')
      p_seq = p_seq + ggtitle(title)
    }
    
    print(p_seq)
    
    output_file_jpg = file.path(output_dir, paste0(output_prefix, "_", seq, ".jpg"))
    output_file_svg = file.path(output_dir, paste0(output_prefix, "_", seq, ".svg"))
    ggsave(output_file_jpg, plot = p_seq, width = 5, height = 3, dpi = 600)
    ggsave(output_file_svg, plot = p_seq, width = 5, height = 3, dpi = 600)
  }
  
  if (!is.null(color_palette_sequences) && !is.null(all_colors)) {
    p_combined = ggplot(plot_df, aes(x = tSNE_1, y = tSNE_2, color = Sequence, shape = Terminal_Type, size = Terminal_Type)) +
      geom_point(alpha = 0.8) +
      scale_color_manual(values = color_palette_sequences) +
      scale_shape_manual(values = c(16, 17), name = '终端类型') +
      scale_size_manual(values = c(3, 4), name = '终端类型') +
      ggtitle('t-SNE可视化：所有目标序列（区分终端类型）') +
      xlim(-20, 20) + ylim(-20, 20) +
      labs(x = 'tSNE_1', y = 'tSNE_2') +
      theme_minimal(base_size = 14) +
      theme(
        panel.grid = element_blank(),
        legend.position = 'right',
        legend.box = 'vertical',
        plot.title = element_text(hjust = 0.5)
      )
    
    print(p_combined)
    
    output_file_jpg = file.path(output_dir, paste0(output_prefix, "_combined.jpg"))
    output_file_svg = file.path(output_dir, paste0(output_prefix, "_combined.svg"))
    ggsave(output_file_jpg, plot = p_combined, width = 6, height = 5.2, dpi = 600)
    ggsave(output_file_svg, plot = p_combined, width = 6, height = 5.2, dpi = 600)
    
    plot_df$Sequence = factor(plot_df$Sequence, levels = c("G", "GG", "GS", "GGG", "S", "SG", "GSG", "GGS"))
    p_facet = ggplot(plot_df, aes(x = tSNE_1, y = tSNE_2, color = Sequence_Type, shape = Terminal_Type)) +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_manual(values = all_colors) +
      scale_shape_manual(values = c(16, 17), name = '终端类型') +
      facet_grid(Terminal_Type ~ Sequence) +
      ggtitle('t-SNE可视化：按终端类型和序列分面') +
      xlim(-20, 20) + ylim(-20, 20) +
      labs(x = 'tSNE_1', y = 'tSNE_2') +
      theme_minimal(base_size = 10) +
      theme(
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)
      )
    
    print(p_facet)
    
    output_file_jpg = file.path(output_dir, paste0(output_prefix, "_faceted.jpg"))
    output_file_svg = file.path(output_dir, paste0(output_prefix, "_faceted.svg"))
    ggsave(output_file_jpg, plot = p_facet, width = 26, height = 5.2, dpi = 600)
    ggsave(output_file_svg, plot = p_facet, width = 26, height = 5.2, dpi = 600)
  }
  
  if (verbose) {
    cat("散点图绘制完成\n")
  }
  
  invisible(NULL)
}

#' Assembly Direction Map Visualization
#'
#' This function creates comprehensive visualizations for differential embeddings.
#'
#' @param differential_df Data frame containing differential embeddings
#' @param target_sequences Target sequences to visualize
#' @param output_dir Directory for output files
#' @param output_prefix Prefix for output filenames
#' @param color_palette_terminal Color palette for terminal types
#' @param color_palette_sequences Color palette for sequences
#' @param heatmap_vis Logical indicating whether to generate heatmap (default: FALSE)
#' @param boxplot_vis Logical indicating whether to generate boxplots (default: FALSE)
#' @param scatterplot_vis Logical indicating whether to generate scatterplots (default: TRUE)
#' @param verbose Logical indicating whether to show progress messages (default: TRUE)
#'
#' @return List containing all generated visualizations
#' @export
# Assembly direction map function | 组装方向图函数
AssemblyDirectionMap <- function(
  differential_df,
  target_sequences,
  output_dir,
  output_prefix = "assembly_map",
  color_palette_terminal = c("Nter" = "#FF0000", "Cter" = "#0000FF"),
  color_palette_sequences = NULL,
  heatmap_vis = FALSE,
  boxplot_vis = FALSE,
  scatterplot_vis = TRUE,
  verbose = TRUE
) {
  result_list = list()
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  if (scatterplot_vis) {
    scatterplot_dir = file.path(output_dir, "scatterplot")
    if (!dir.exists(scatterplot_dir)) {
      dir.create(scatterplot_dir, recursive = TRUE)
    }
    result_list$scatterplot = VisualizeScatterplots(
      differential_df = differential_df,
      target_sequences = target_sequences,
      output_dir = scatterplot_dir,
      output_prefix = output_prefix,
      color_palette_terminal = color_palette_terminal,
      color_palette_sequences = color_palette_sequences,
      verbose = verbose
    )
  }
  
  if (heatmap_vis) {
    heatmap_dir = file.path(output_dir, "heatmap")
    if (!dir.exists(heatmap_dir)) {
      dir.create(heatmap_dir, recursive = TRUE)
    }
    result_list$heatmap = VisualizeDifferentialHeatmap(
      differential_df = differential_df,
      output_dir = heatmap_dir,
      output_prefix = output_prefix,
      verbose = verbose
    )
  }
  
  if (boxplot_vis) {
    boxplot_dir = file.path(output_dir, "boxplot")
    if (!dir.exists(boxplot_dir)) {
      dir.create(boxplot_dir, recursive = TRUE)
    }
    result_list$boxplot = VisualizePositionBoxplots(
      differential_df = differential_df,
      positions_to_visualize = c(1, 2, 3),
      output_dir = boxplot_dir,
      output_prefix = output_prefix,
      verbose = verbose
    )
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
cat("  - VisualizeDifferentialHeatmap() : Visualize differential embeddings using heatmap\n")
cat("  - VisualizePositionBoxplots()  : Visualize differential embeddings using boxplots\n")
cat("  - VisualizeScatterplots()      : Visualize differential embeddings using scatterplots\n")
cat("  - AssemblyDirectionMap()       : Comprehensive differential embeddings visualization\n")
cat("\n")
cat("For help, see ?plot_mutation_positions or ?StructuralMapVis or ?AssemblyDirectionMap\n")
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
