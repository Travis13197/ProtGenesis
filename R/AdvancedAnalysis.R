# ==============================================================================
# AdvancedAnalysis.R
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
# ==============================================================================

check_and_install_packages_adv = function() {
  required_packages = c("Seurat", "data.table", "dplyr", "plotly", 
                        "ggplot2", "htmlwidgets", "reshape2")
  
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
  }
  
  library(Seurat, quietly = TRUE)
  library(data.table, quietly = TRUE)
  library(dplyr, quietly = TRUE)
  library(plotly, quietly = TRUE)
  library(ggplot2, quietly = TRUE)
  library(htmlwidgets, quietly = TRUE)
  library(reshape2, quietly = TRUE)
}

check_and_install_packages_adv()

read_embedding_data = function(file_path) {
  raw_data <- data.table::fread(file_path) %>% as.data.frame()
  sequence_ids <- raw_data[, 1]
  embed_matrix <- as.matrix(raw_data[, -1])
  return(list(sequence_ids = sequence_ids, embed_matrix = embed_matrix))
}

classify_sequence = function(seq_id) {
  if (grepl("^EAC_", seq_id)) return("EAC")
  else if (grepl("^MOTIF_", seq_id)) return("MOTIF")
  else if (grepl("^POLY_", seq_id)) return("POLY")
  else return("Unknown")
}

filter_category_data = function(sequence_ids, embed_matrix, target_category) {
  categories <- sapply(sequence_ids, classify_sequence)
  target_indices <- which(categories == target_category)
  return(list(
    sequence_ids = sequence_ids[target_indices],
    embed_matrix = embed_matrix[target_indices, ]
  ))
}

extract_step_from_id = function(seq_id) {
  step_match <- regexec("_Step(\\d+)$", seq_id)
  if (length(step_match[[1]]) > 1) {
    return(as.integer(regmatches(seq_id, step_match)[[1]][2]))
  }
  return(NA)
}

extract_poly_category = function(seq_id) {
  if (grepl("^POLY_", seq_id)) {
    poly_match <- regexec("^POLY_(\\w+)_", seq_id)
    if (length(poly_match[[1]]) > 1) {
      return(regmatches(seq_id, poly_match)[[1]][2])
    }
  }
  return("Unknown")
}

build_seurat_object = function(sequence_ids, embed_matrix, category) {
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
  
  protgenesis_obj <- CreateSeuratObject(
    counts = sparse_matrix,
    assay = "embed",
    project = paste0("ProtGenesis_", category),
    min.cells = 0,
    min.features = 0
  )
  
  protgenesis_obj <- protgenesis_obj %>% AddMetaData(data_categories)
  
  protgenesis_obj <- SetAssayData(
    object = protgenesis_obj,
    assay = "embed",
    layer = "data",
    new.data = GetAssayData(protgenesis_obj, assay = "embed", layer = "counts")
  )
  
  protgenesis_obj <- ScaleData(
    object = protgenesis_obj,
    assay = "embed",
    do.center = TRUE,
    do.scale = FALSE
  )
  
  return(protgenesis_obj)
}

perform_dim_reduction = function(protgenesis_obj) {
  protgenesis_obj <- RunPCA(
    object = protgenesis_obj,
    assay = "embed",
    features = rownames(protgenesis_obj),
    npcs = min(50, ncol(protgenesis_obj) - 1),
    verbose = FALSE
  )
  
  pct_variance <- protgenesis_obj@reductions$pca@stdev^2 / sum(protgenesis_obj@reductions$pca@stdev^2)
  selected_dims <- which(cumsum(pct_variance) < 0.95)
  
  n_cells <- ncol(protgenesis_obj)
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
  
  tsne_perplexity <- min(100, floor((n_cells - 1) / 3))
  protgenesis_obj <- RunTSNE(
    object = protgenesis_obj,
    dims = selected_dims,
    reduction = "pca",
    perplexity = tsne_perplexity,
    check_duplicates = FALSE
  )
  
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

get_category_colors = function() {
  list(
    category_colors = c("EAC" = "#DE6038", "MOTIF" = "#5778AA", "POLY" = "#46A64B", "Unknown" = "#ADADAD"),
    step_colors = c("Step1" = "#FF0000", "Step2" = "#FF6600", "Step5" = "#FFCC00", "Step10" = "#00FF00",
                    "Step20" = "#00CCFF", "Step40" = "#0066FF", "Step60" = "#6600FF", "Step80" = "#9900FF",
                    "Step100" = "#CC00FF", "Other" = "#ADADAD")
  )
}

generate_step_based_plots = function(protgenesis_obj, category, output_dir) {
  p_umap <- DimPlot(object = protgenesis_obj, group.by = "Step", pt.size = 0.5, label = FALSE, reduction = "umap") + 
    ggtitle(paste(category, "Step Distribution (UMAP)")) + theme_minimal(base_size = 14)
  ggsave(file.path(output_dir, paste0(category, "_step_umap.jpg")), height = 8, width = 10.4, dpi = 600)
  
  p_tsne <- DimPlot(object = protgenesis_obj, group.by = "Step", pt.size = 0.5, label = FALSE, reduction = "tsne") + 
    ggtitle(paste(category, "Step Distribution (t-SNE)")) + theme_minimal(base_size = 14)
  ggsave(file.path(output_dir, paste0(category, "_step_tsne.jpg")), height = 8, width = 10.4, dpi = 600)
  
  return(list(umap = p_umap, tsne = p_tsne))
}

generate_poly_category_plots = function(protgenesis_obj, output_dir) {
  dir.create(file.path(output_dir, "poly_category_visualizations"), recursive = TRUE, showWarnings = FALSE)
  
  poly_colors <- c("G" = "#4CAF50", "A" = "#2196F3", "D" = "#FF9800", "E" = "#9C27B0", "V" = "#F44336", "Unknown" = "#9E9E9E")
  
  p <- DimPlot(object = protgenesis_obj, group.by = "Poly_Category", pt.size = 0.5, label = FALSE, reduction = "umap") + 
    scale_color_manual(values = poly_colors) + ggtitle("POLY Subcategory Distribution (UMAP)") + theme_minimal(base_size = 14)
  ggsave(file.path(output_dir, "poly_category_visualizations", "POLY_category_umap.jpg"), height = 8, width = 10.4, dpi = 600)
  
  return(p)
}

generate_3d_visualizations = function(protgenesis_obj, category, output_dir) {
  dir.create(file.path(output_dir, "3d_visualizations"), recursive = TRUE, showWarnings = FALSE)
  meta_data <- protgenesis_obj@meta.data
  
  if ("umap3d" %in% names(protgenesis_obj@reductions)) {
    umap3d_embed <- Embeddings(protgenesis_obj, reduction = "umap3d")
    if (ncol(umap3d_embed) >= 3) {
      umap3d_data <- data.frame(
        UMAP1 = umap3d_embed[, 1], UMAP2 = umap3d_embed[, 2], UMAP3 = umap3d_embed[, 3],
        Step = meta_data$Step, Sequence_ID = meta_data$Sequence_ID
      )
      
      p <- plot_ly(data = umap3d_data, x = ~UMAP1, y = ~UMAP2, z = ~UMAP3, color = ~as.factor(Step),
                   text = ~paste("Sequence:", Sequence_ID, "<br>Step:", Step),
                   hoverinfo = "text", type = "scatter3d", mode = "markers", marker = list(size = 3)) %>%
        layout(title = paste(category, "3D UMAP"), scene = list(xaxis = list(title = "UMAP 1"), yaxis = list(title = "UMAP 2"), zaxis = list(title = "UMAP 3")))
      
      htmlwidgets::saveWidget(p, file = file.path(output_dir, "3d_visualizations", paste0(category, "_3d_umap.html")), selfcontained = FALSE)
    }
  }
}

generate_interactive_plot = function(protgenesis_obj, category, output_dir) {
  plot_data <- Seurat::Embeddings(protgenesis_obj, reduction = "tsne") %>% as.data.frame() %>% 
    mutate(Cell = rownames(.)) %>% left_join(protgenesis_obj@meta.data[, c("Cell", "Step", "Sequence_ID")], by = "Cell")
  
  p <- plot_ly(data = plot_data, x = ~tSNE_1, y = ~tSNE_2, color = ~factor(Step),
               text = ~paste("Sequence ID:", Sequence_ID, "\nStep:", Step), hoverinfo = "text+x+y",
               opacity = 0.7, marker = list(size = 3)) %>%
    layout(title = paste(category, "Interactive Visualization"), xaxis = list(title = "tSNE_1"), yaxis = list(title = "tSNE_2"))
  
  htmlwidgets::saveWidget(p, file = file.path(output_dir, paste0(category, "_interactive.html")), selfcontained = FALSE)
}

analyze_dataset = function(file_path, category, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  data <- read_embedding_data(file_path)
  filtered_data <- filter_category_data(data$sequence_ids, data$embed_matrix, category)
  protgenesis_obj <- build_seurat_object(filtered_data$sequence_ids, filtered_data$embed_matrix, category)
  protgenesis_obj <- perform_dim_reduction(protgenesis_obj)
  
  save(protgenesis_obj, file = file.path(output_dir, paste0(category, "_protgenesis_object.rdata")))
  
  generate_step_based_plots(protgenesis_obj, category, output_dir)
  generate_3d_visualizations(protgenesis_obj, category, output_dir)
  generate_interactive_plot(protgenesis_obj, category, output_dir)
  
  if (category == "POLY") {
    generate_poly_category_plots(protgenesis_obj, output_dir)
  }
  
  return(protgenesis_obj)
}

analyze_all_datasets = function(file_path, base_output_dir) {
  dir.create(base_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  eac_obj <- analyze_dataset(file_path, "EAC", file.path(base_output_dir, "EAC"))
  motif_obj <- analyze_dataset(file_path, "MOTIF", file.path(base_output_dir, "MOTIF"))
  poly_obj <- analyze_dataset(file_path, "POLY", file.path(base_output_dir, "POLY"))
  
  return(list(EAC = eac_obj, MOTIF = motif_obj, POLY = poly_obj))
}
