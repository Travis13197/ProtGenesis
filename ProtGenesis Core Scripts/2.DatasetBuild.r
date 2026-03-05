# ==============================================================================
# 2.DatasetBuild.r | 数据集构建函数
# ==============================================================================
#
# Purpose: Dimensionality reduction and dataset construction functions
#   This script provides functions for dimensionality reduction of protein
#   embedding data, including PCA, UMAP, and tSNE methods.
#   本脚本提供蛋白质嵌入数据降维的函数，包括 PCA、UMAP 和 tSNE 方法。
#
# Author: ProtGenesis Analysis Team
# Date: 2026-02-24
# Version: 1.0.0
# By: Chuanyang Liu
# Please cite our paper if you use this script in your research:
#   Chuanyang Liu, et al. (2026). Universal physical principles govern the deterministic genesis of protein structure.
#   BioRxiv.10.64898/2026.02.20.706798.
#
# Usage: | 使用方法
#   source("2.DatasetBuild.r")
#   protgenesis_obj <- ReduceDim(protgenesis_obj)
#
# Requirements: | 依赖项
#   - Seurat (for single-cell analysis and dimensionality reduction)
#
# ==============================================================================

# ==============================================================================
# Package Management | 包管理
# ==============================================================================

check_and_install_packages <- function() {
  required_packages <- c("Seurat", "Matrix", "data.table", "dplyr", "stringr", "plotly", "htmlwidgets")

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

#' Perform dimensionality reduction on ProtGenesis object
#'
#' This function performs comprehensive dimensionality reduction on a ProtGenesis
#' object, including PCA, UMAP, and tSNE. It includes intelligent dimension
#' selection and robust error handling.
#'
#' @param protgenesis_obj ProtGenesis object containing protein embedding data
#' @param pca_npcs Number of PCA components to compute (default: 50)
#' @param variance_cutoff Cumulative variance cutoff for dimension selection (default: 0.95)
#' @param umap_dist Minimum distance for UMAP (default: 0.3)
#' @param umap_spread Spread parameter for UMAP (default: 1.5)
#' @param umap_neighbors Number of neighbors for UMAP (default: 50)
#' @param tsne_perplexity Perplexity parameter for tSNE (default: 100)
#'
#' @return ProtGenesis object with dimensionality reductions added
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage with default parameters
#' protgenesis_obj = ReduceDim(protgenesis_obj)
#'
#' # Custom parameters
#' protgenesis_obj <- ReduceDim(protgenesis_obj, pca_npcs = 50, umap_dist = 0.5)
#' }
# Main function: Dimensionality reduction | 主函数：降维
ReduceDim <- function(protgenesis_obj,
                      pca_npcs = 50,
                      variance_cutoff = 0.95,
                      umap_dist = 0.3,
                      umap_spread = 1.5,
                      umap_neighbors = 50,
                      tsne_perplexity = 100) {

  protgenesis_obj = FindVariableFeatures(protgenesis_obj, selection.method = "vst", nfeatures = 1000)

  protgenesis_obj = RunPCA(
    object = protgenesis_obj,
    assay = "embed",
    features = rownames(protgenesis_obj),
    npcs = min(pca_npcs, ncol(protgenesis_obj) - 1),
    verbose = FALSE
  )

  tryCatch({
    pct_variance = protgenesis_obj@reductions$pca@stdev^2 / sum(protgenesis_obj@reductions$pca@stdev^2)
    selected_dims = which(cumsum(pct_variance) < variance_cutoff)

    if (length(selected_dims) == 0) {
      warning("Cannot select dimensions using cumulative variance, using default dimensions")
      selected_dims = 1:min(20, pca_npcs, ncol(protgenesis_obj) - 1)
    }
  }, error = function(e) {
    warning(paste("Dimension selection error:", e$message, "using default dimensions"))
    selected_dims = 1:min(20, pca_npcs, ncol(protgenesis_obj) - 1)
  })

  if (length(selected_dims) < 2) {
    warning("Insufficient dimensions, expanding to at least 2 dimensions")
    selected_dims = 1:min(2, pca_npcs, ncol(protgenesis_obj) - 1)
  }

  tryCatch({
    protgenesis_obj = RunUMAP(
      object = protgenesis_obj,
      dims = selected_dims,
      reduction = "pca",
      metric = "cosine",
      min.dist = umap_dist,
      spread = umap_spread,
      n.neighbors = umap_neighbors,
      learning.rate = 0.8,
      n.components = 3
    )
  }, error = function(e) {
    warning(paste("UMAP error:", e$message, "trying alternative parameters"))
    protgenesis_obj = RunUMAP(
      object = protgenesis_obj,
      dims = selected_dims,
      reduction = "pca",
      metric = "euclidean",
      min.dist = 0.5,
      spread = 1.0,
      n.neighbors = 50,
      learning.rate = 0.5,
      verbose = FALSE,
      n.components = 3
    )
  })

  tryCatch({
    protgenesis_obj = RunTSNE(
      object = protgenesis_obj,
      dims = selected_dims,
      reduction = "pca",
      perplexity = min(tsne_perplexity, length(selected_dims) - 1),
      check_duplicates = FALSE,
      dim.embed = 3
    )
  }, error = function(e) {
    warning(paste("tSNE error:", e$message))
  })

  return(protgenesis_obj)
}

#' Build ProtGenesis Dataset from protein embedding matrix
#'
#' This function creates a ProtGenesis object (Seurat-based) from a protein
#' embedding matrix file. It handles reading the matrix, transposing,
#' converting to sparse matrix, creating ProtGenesis object, and normalization.
#'
#' @param embedding_file Path to embedding matrix file (csv.gz format)
#' @param assay_name Name of the assay in ProtGenesis object (default: "embed")
#' @param project_name Name of the project (default: "ProteinEmbed")
#' @param do_center Logical indicating whether to center the data (default: TRUE)
#' @param do_scale Logical indicating whether to scale the data (default: FALSE)
#' @param min.cells Minimum cells threshold for CreateSeuratObject (default: 0)
#' @param min.features Minimum features threshold for CreateSeuratObject (default: 0)
#' @param verbose Logical indicating whether to print progress messages (default: TRUE)
#'
#' @return A ProtGenesis object containing the protein embedding data
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' protgenesis_obj = BuildProtGenesisDataset("protein_embedding_matrix_16wHuman.csv.gz")
#'
#' # Custom parameters
#' protgenesis_obj = BuildProtGenesisDataset(
#'   "protein_embedding_matrix_16wHuman.csv.gz",
#'   assay_name = "protein",
#'   project_name = "MyProject"
#' )
#' }
# Build dataset function | 构建数据集函数
BuildProtGenesisDataset <- function(embedding_file,
                                      assay_name = "embed",
                                      project_name = "ProteinEmbed",
                                      do_center = TRUE,
                                      do_scale = FALSE,
                                      min.cells = 0,
                                      min.features = 0,
                                      verbose = TRUE) {

  if (!file.exists(embedding_file)) {
    stop(sprintf("Embedding file not found: %s", embedding_file))
  }

  raw_matrix = fread(embedding_file) %>%
    as.data.frame() %>%
    as.matrix()

  dat = as.data.frame(t(raw_matrix))
  colnames(dat) = dat[1, ]
  dat = dat[-1, ]

  dat_matrix = as.matrix(dat)
  storage.mode(dat_matrix) = "numeric"

  sparse_matrix = Matrix(
    data = dat_matrix,
    nrow = nrow(dat),
    ncol = ncol(dat),
    dimnames = list(rownames(dat), colnames(dat)),
    sparse = TRUE
  )

  protgenesis_obj = CreateSeuratObject(
    counts = sparse_matrix,
    assay = assay_name,
    project = project_name,
    min.cells = min.cells,
    min.features = min.features,
    names.field = 1,
    names.delim = "-"
  )

  protgenesis_obj = SetAssayData(
    object = protgenesis_obj,
    assay = assay_name,
    layer = "data",
    new.data = GetAssayData(protgenesis_obj, assay = assay_name, layer = "counts")
  )

  protgenesis_obj = ScaleData(
    object = protgenesis_obj,
    assay = assay_name,
    do.center = do_center,
    do.scale = do_scale
  )

  if (verbose) {
    cat("=== ProtGenesis Dataset Complete ===\n")
    cat(sprintf("Number of proteins: %d\n", ncol(protgenesis_obj)))
    cat(sprintf("Number of embedding dimensions: %d\n", nrow(protgenesis_obj)))
  }

  return(protgenesis_obj)
}

#' Extract metadata from ProtGenesis sequence IDs
#'
#' This function extracts comprehensive metadata from sequence IDs generated
#' by the DataGen functions (DataGen_Mutation and DataGen_Stepwise). It
#' parses the ID format to extract protein name, sequence index, segment type,
#' position, amino acid, direction, and segment information.
#'
#' @param ids Vector of sequence IDs to parse
#' @param verbose Logical indicating whether to print progress messages (default: TRUE)
#'
#' @return A data.frame containing extracted metadata with columns:
#'   \item{Protein_Name}{Original protein name}
#'   \item{Seq_Index}{Sequence index}
#'   \item{Seg_Type}{Segment type (ORIG, MUT, STEP)}
#'   \item{Position}{Position in sequence (if applicable)}
#'   \item{Amino_Acid}{Amino acid at position (if applicable)}
#'   \item{Direction}{Direction (Fwd or Rev)}
#'   \item{Seg_Start}{Segment start position}
#'   \item{Seg_End}{Segment end position}
#'   \item{Round}{Round information (from legacy format, if applicable)}
#' @export
#'
#' @examples
#' \dontrun{
#' # Example IDs from generate_sequence_name()
#' ids = c(
#'   "TestProtein_Seq1_ORIG_Fwd_Seg1-50",
#'   "TestProtein_Seq1_MUT_Pos15_A_Fwd_Seg1-50",
#'   "TestProtein_Seq1_STEP_Pos20_K_Fwd_Seg1-50"
#' )
#' metadata <- MetadataExtract(ids)
#' }
# Metadata extraction function | 元数据提取函数
MetadataExtract <- function(ids, verbose = TRUE) {

  metadata_list = lapply(ids, function(id) {
    result = list(
      Protein_Name = NA,
      Seq_Index = NA,
      Seg_Type = NA,
      Position = NA,
      Amino_Acid = NA,
      Direction = NA,
      Seg_Start = NA,
      Seg_End = NA,
      Round = NA,
      Original_ID = id
    )

    parts = strsplit(id, "_")[[1]]

    if (grepl("^([1-4])_", id)) {
      round_match = regmatches(id, regexec("^([1-4])_", id))
      if (length(round_match[[1]]) > 1) {
        result$Round = paste0("Round", round_match[[1]][2])
      }
    }

    if (grepl("_GFP_([1-4])_", id)) {
      gfp_round_match = regmatches(id, regexec("_GFP_([1-4])_", id))
      if (length(gfp_round_match[[1]]) > 1) {
        result$Round = paste0("Round", gfp_round_match[[1]][2])
      }
    }

    seq_idx = grep("^Seq\\d+$", parts)
    if (length(seq_idx) > 0) {
      result$Seq_Index = as.integer(sub("Seq", "", parts[seq_idx[1]]))

      if (seq_idx[1] > 1) {
        result$Protein_Name = paste(parts[1:(seq_idx[1]-1)], collapse = "_")
      }

      if (seq_idx[1] < length(parts)) {
        seg_type_idx = seq_idx[1] + 1
        if (seg_type_idx <= length(parts)) {
          result$Seg_Type = parts[seg_type_idx]

          pos_idx = grep("^Pos\\d+$", parts)
          if (length(pos_idx) > 0 && pos_idx[1] > seg_type_idx) {
            result$Position = as.integer(sub("Pos", "", parts[pos_idx[1]]))

            if (pos_idx[1] + 1 <= length(parts)) {
              aa_candidate = parts[pos_idx[1] + 1]
              valid_aa = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                            "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
              if (aa_candidate %in% valid_aa) {
                result$Amino_Acid = aa_candidate
              }
            }
          }

          dir_idx = which(parts %in% c("Fwd", "Rev"))
          if (length(dir_idx) > 0) {
            result$Direction = parts[dir_idx[1]]
          }

          seg_idx = grep("^Seg\\d+-\\d+$", parts)
          if (length(seg_idx) > 0) {
            seg_parts = strsplit(sub("Seg", "", parts[seg_idx[1]]), "-")[[1]]
            result$Seg_Start = as.integer(seg_parts[1])
            result$Seg_End = as.integer(seg_parts[2])
          }
        }
      }
    }

    as.data.frame(result, stringsAsFactors = FALSE)
  })

  metadata_df = do.call(rbind, metadata_list)
  rownames(metadata_df) = ids

  if (verbose) {
    cat("=== Metadata Extraction Complete ===\n")
    cat(sprintf("Protein names found: %d unique\n", length(unique(metadata_df$Protein_Name[!is.na(metadata_df$Protein_Name)]))))
    cat(sprintf("Sequence types: %s\n", paste(unique(metadata_df$Seg_Type[!is.na(metadata_df$Seg_Type)]), collapse = ", ")))
    cat(sprintf("Round info found: %d\n", sum(!is.na(metadata_df$Round))))
  }

  return(metadata_df)
}

#' Perform clustering analysis on ProtGenesis object
#'
#' This function performs comprehensive clustering analysis including
#' FindNeighbors, FindClusters, and optional visualization.
#'
#' @param protgenesis_obj ProtGenesis object containing protein embedding data
#' @param dims Dimension range for PCA (default: 1:50)
#' @param k_param Number of neighbors for FindNeighbors (default: 50)
#' @param annoy_metric Distance metric for FindNeighbors (default: "cosine")
#' @param prune_snn Pruning parameter for SNN graph (default: 1/15)
#' @param resolution Clustering resolution (default: 3)
#' @param algorithm Clustering algorithm (1-4, default: 4)
#' @param random_seed Random seed for reproducibility (default: 42)
#' @param run_clustering Logical indicating whether to perform clustering (default: TRUE)
#' @param run_visualization Logical indicating whether to generate visualization (default: TRUE)
#' @param group_by Metadata column to group by for visualization
#' @param cols Color vector for groups
#' @param pt_size Point size for visualization (default: 0.02)
#' @param label Logical indicating whether to show labels (default: TRUE)
#' @param raster Logical indicating whether to use rasterization (default: TRUE)
#' @param reduction Reduction type for visualization (default: "tsne")
#' @param title Plot title
#' @param output_prefix Prefix for output filenames (default: "cluster_plot")
#' @param output_dir Directory for output files (default: "Figures")
#' @param output_svg Logical indicating whether to output SVG format (default: TRUE)
#' @param output_jpg Logical indicating whether to output JPG format (default: TRUE)
#' @param output_html Logical indicating whether to output HTML format (default: TRUE)
#' @param height Plot height in inches (default: 6)
#' @param width Plot width in inches (default: 7)
#' @param dpi Resolution for output (default: 600)
#' @param verbose Logical indicating whether to print progress messages (default: TRUE)
#'
#' @return List containing updated protgenesis_obj and cluster_plot (if generated)
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' result = ClusterFind(protgenesis_obj)
#'
#' # With visualization
#' result = ClusterFind(
#'   protgenesis_obj,
#'   group_by = "Round_Info",
#'   cols = c("Round1" = "#DE6038", "Round2" = "#D99829")
#' )
#' }
# Clustering function | 聚类函数
ClusterFind <- function(
  protgenesis_obj,
  dims = 1:50,
  k_param = 50,
  annoy_metric = "cosine",
  prune_snn = 1/15,
  resolution = 3,
  algorithm = 4,
  random_seed = 42,
  run_clustering = TRUE,
  run_visualization = TRUE,
  group_by = NULL,
  cols = NULL,
  pt_size = 0.02,
  label = FALSE,
  raster = FALSE,
  reduction = "tsne",
  title = NULL,
  output_prefix = "cluster_plot",
  output_dir = "Figures",
  output_svg = TRUE,
  output_jpg = TRUE,
  output_html = TRUE,
  height = 6,
  width = 7,
  dpi = 600,
  verbose = TRUE
) {

  result_list = list(
    protgenesis_obj = protgenesis_obj,
    cluster_plot = NULL
  )

  if (run_clustering) {
    protgenesis_obj = FindNeighbors(
      object = protgenesis_obj,
      dims = dims,
      k.param = k_param,
      annoy.metric = annoy_metric,
      prune.SNN = prune_snn
    )

    protgenesis_obj = FindClusters(
      protgenesis_obj,
      resolution = resolution,
      algorithm = algorithm,
      random.seed = random_seed
    )

    result_list$protgenesis_obj = protgenesis_obj
  }

  if (run_visualization) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }

    plot_args = list(
      object = protgenesis_obj,
      reduction = reduction,
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

    p_cluster = do.call(DimPlot, plot_args)

    if (!is.null(title)) {
      p_cluster = p_cluster + ggtitle(title)
    } else {
      p_cluster = p_cluster + ggtitle("Clustering Analysis")
    }

    p_cluster = p_cluster + theme_minimal(base_size = 14)

    print(p_cluster)
    result_list$cluster_plot = p_cluster

    if (output_jpg) {
      jpg_file = file.path(output_dir, sprintf("%s.jpg", output_prefix))
      ggsave(jpg_file, plot = p_cluster, height = height, width = width, dpi = dpi)
    }

    if (output_svg) {
      svg_file = file.path(output_dir, sprintf("%s.svg", output_prefix))
      ggsave(svg_file, plot = p_cluster, height = height, width = width, dpi = dpi)
    }

    if (output_html) {
      html_file = file.path(output_dir, sprintf("%s.html", output_prefix))
      plotly_obj = ggplotly(p_cluster)
      htmlwidgets::saveWidget(plotly_obj, file = html_file, selfcontained = TRUE)
    }
  }

  if (verbose) {
    cat("=== ClusterFind Complete ===\n")
  }

  return(invisible(result_list))
}

#' Extract metadata for GFP Variants sequence IDs
#'
#' This function extracts metadata from GFP Variants sequence IDs,
#' including main category, mutation position, linker information, etc.
#'
#' @param ids Vector of sequence IDs to parse
#' @param verbose Logical indicating whether to print progress messages (default: TRUE)
#'
#' @return A data.frame containing extracted metadata
#' @export
# GFP variants metadata extraction | GFP 变体元数据提取
MetadataExtract_GFPVariants <- function(ids, verbose = TRUE) {

  metadata_list = lapply(ids, function(id) {
    result = list(
      Cell = id,
      Main_Category = NA,
      Is_Original = NA,
      Structure_Type = NA,
      Mutation_Domain = NA,
      Linker_Length = NA,
      Linker_Composition = NA,
      Original_ID = id
    )

    if (grepl("GFPMut", id)) {
      result$Main_Category = "GFPMut"
    } else if (grepl("^GFP_L", id)) {
      result$Main_Category = "GFPlinker"
    } else if (grepl("^RevGFP_L", id)) {
      result$Main_Category = "RevGFPlinker"
    } else if (grepl("StepGFP", id)) {
      result$Main_Category = "GFPstepwisemut"
    } else if (grepl("Trim21", id)) {
      result$Main_Category = "TRIM21mut"
    } else {
      result$Main_Category = "Other"
    }

    result$Is_Original = grepl("^ORI_", id)

    if (grepl("^(GFP|RevGFP)_(GFP|RevGFP)(_Rev)?$", id)) {
      result$Structure_Type = "Basic-vari"
    } else if (grepl("^(GFP|RevGFP)_L([1-5])_", id)) {
      prefix_match = regmatches(id, regexec("^(GFP|RevGFP)", id))
      linker_match = regmatches(id, regexec("_L([1-5])_", id))
      if (length(prefix_match[[1]]) > 1 && length(linker_match[[1]]) > 1) {
        result$Structure_Type = paste0(prefix_match[[1]][2], "-L", linker_match[[1]][2])
      }
    } else {
      result$Structure_Type = "Others"
    }

    pos_match = regmatches(id, regexec("(?<=Mut_)\\d{3}", id, perl = TRUE))
    if (length(pos_match[[1]]) > 0) {
      pos_num = as.integer(pos_match[[1]][1])
      result$Mutation_Domain = sprintf("MutPos_%03d", pos_num)
    }

    linker_match = regmatches(id, regexec("(?<=_L)[1-5](?=_)", id, perl = TRUE))
    if (length(linker_match[[1]]) > 0) {
      result$Linker_Length = linker_match[[1]][1]
    }

    aa_match = regmatches(id, regexec("L[1-5]_(\\w+)_", id))
    if (length(aa_match[[1]]) > 1) {
      result$Linker_Composition = aa_match[[1]][2]
    }

    as.data.frame(result, stringsAsFactors = FALSE)
  })

  metadata_df = do.call(rbind, metadata_list)
  rownames(metadata_df) = ids

  if (verbose) {
    cat("=== MetadataExtraction_GFPVariants Complete ===\n")
    cat(sprintf("Main categories found: %s\n", paste(unique(metadata_df$Main_Category[!is.na(metadata_df$Main_Category)]), collapse = ", ")))
    cat(sprintf("Structure types: %s\n", paste(unique(metadata_df$Structure_Type[!is.na(metadata_df$Structure_Type)]), collapse = ", ")))
    cat(sprintf("Linker lengths found: %s\n", paste(unique(metadata_df$Linker_Length[!is.na(metadata_df$Linker_Length)]), collapse = ", ")))
  }

  return(metadata_df)
}

#' Extract metadata for TRIM21 sequence IDs
#'
#' This function extracts metadata from TRIM21 sequence IDs,
#' including species category, ORI status, and position number.
#'
#' @param ids Vector of sequence IDs to parse
#' @param verbose Logical indicating whether to print progress messages (default: TRUE)
#'
#' @return A data.frame containing extracted metadata
#' @export
# TRIM21 metadata extraction | TRIM21 元数据提取
MetadataExtract_TRIM21 <- function(ids, verbose = TRUE) {

  metadata_list = lapply(ids, function(id) {
    result = list(
      Cell = id,
      Species_Category = NA,
      ORI_Status = NA,
      Position_Number = NA,
      Original_ID = id
    )

    if (grepl("^(ORI_)?Trim21_\\d+", id)) {
      result$Species_Category = "HumanTRIM21"
    } else if (grepl("^(ORI_)?Trim21_mouse_", id)) {
      result$Species_Category = "MouseTRIM21"
    } else if (grepl("^(ORI_)?Trim5_buman_", id)) {
      result$Species_Category = "HumanTRIM5"
    } else {
      result$Species_Category = "Other"
    }

    if (grepl("^ORI_", id)) {
      result$ORI_Status = "ORI"
    } else {
      result$ORI_Status = "Mut"
    }

    pos_match = regmatches(id, regexec("(?<=_)\\d+(?=_)", id, perl = TRUE))
    if (length(pos_match[[1]]) > 0) {
      result$Position_Number = pos_match[[1]][1]
    }

    as.data.frame(result, stringsAsFactors = FALSE)
  })

  metadata_df = do.call(rbind, metadata_list)
  rownames(metadata_df) = ids

  if (verbose) {
    cat("=== MetadataExtract_TRIM21 Complete ===\n")
    cat(sprintf("Species categories: %s\n", paste(unique(metadata_df$Species_Category[!is.na(metadata_df$Species_Category)]), collapse = ", ")))
    cat(sprintf("ORI status: %s\n", paste(unique(metadata_df$ORI_Status[!is.na(metadata_df$ORI_Status)]), collapse = ", ")))
  }

  return(metadata_df)
}

#' Extract Metadata for TetR-rtTA Sequence IDs
#'
#' This function extracts comprehensive metadata from TetR-rtTA sequence IDs.
#'
#' @param ids Vector of sequence IDs to parse
#' @param verbose Logical indicating whether to print progress messages (default: TRUE)
#'
#' @return A data.frame containing extracted metadata
#' @export
# TetR metadata extraction | TetR 元数据提取
MetadataExtract_TetR <- function(ids, verbose = TRUE) {

  metadata_list = lapply(ids, function(id) {
    result = list(
      Cell = id,
      Original_Name = id,
      Category = NA,
      Protein_Type = NA,
      ON_OFF_Status = NA,
      Original_Sample = NA,
      Position_Number = NA,
      Amino_Acid = NA,
      Is_Original = NA,
      Original_ID = id
    )

    # Category: Mutations or Stepwise
    if (grepl("GFPMut", id)) {
      result$Category = "Mutations"
    } else if (grepl("ORI_|Step|ORI_", id)) {
      result$Category = "Stepwise"
    } else {
      result$Category = "Other"
    }

    # Protein_Type
    if (grepl("TetR_OFF", id)) {
      result$Protein_Type = "TetR_OFF"
    } else if (grepl("rtTA", id)) {
      result$Protein_Type = "rtTA"
    } else {
      result$Protein_Type = "Other"
    }

    # ON/OFF Status
    if (grepl("TetR_OFF|WT1|WT2|WT3", id)) {
      result$ON_OFF_Status = "OFF"
    } else if (grepl("rtTA|ON1|ON2|ON3", id)) {
      result$ON_OFF_Status = "ON"
    } else {
      result$ON_OFF_Status = "Unknown"
    }

    # Original Sample
    if (grepl("WT1", id)) {
      result$Original_Sample = "WT1"
    } else if (grepl("WT2", id)) {
      result$Original_Sample = "WT2"
    } else if (grepl("WT3", id)) {
      result$Original_Sample = "WT3"
    } else if (grepl("ON1", id)) {
      result$Original_Sample = "ON1"
    } else if (grepl("ON2", id)) {
      result$Original_Sample = "ON2"
    } else if (grepl("ON3", id)) {
      result$Original_Sample = "ON3"
    } else {
      result$Original_Sample = "Other"
    }

    # Is Original
    result$Is_Original = grepl("^ORI_", id)

    # Position Number and Amino Acid (for Mutations)
    if (grepl("GFPMut", id)) {
      pos_match = regmatches(id, regexec("GFPMut_(\\d+)_(\\w+)$", id))
      if (length(pos_match[[1]]) > 2) {
        result$Position_Number = as.integer(pos_match[[1]][2])
        result$Amino_Acid = pos_match[[1]][3]
      }
    }

    # Position Number (for Stepwise)
    if (grepl("ORI_|Step", id)) {
      pos_match = regmatches(id, regexec("_(\\d+)_(\\w+)$", id))
      if (length(pos_match[[1]]) > 2) {
        result$Position_Number = as.integer(pos_match[[1]][2])
        result$Amino_Acid = pos_match[[1]][3]
      }
    }

    as.data.frame(result, stringsAsFactors = FALSE)
  })

  metadata_df = do.call(rbind, metadata_list)
  rownames(metadata_df) = ids

  if (verbose) {
    cat("=== MetadataExtract_TetR Complete ===\n")
    cat(sprintf("Categories found: %s\n", paste(unique(metadata_df$Category[!is.na(metadata_df$Category)]), collapse = ", ")))
    cat(sprintf("Protein types: %s\n", paste(unique(metadata_df$Protein_Type[!is.na(metadata_df$Protein_Type)]), collapse = ", ")))
    cat(sprintf("ON/OFF status: %s\n", paste(unique(metadata_df$ON_OFF_Status[!is.na(metadata_df$ON_OFF_Status)]), collapse = ", ")))
    cat(sprintf("Original samples: %s\n", paste(unique(metadata_df$Original_Sample[!is.na(metadata_df$Original_Sample)]), collapse = ", ")))
  }

  return(metadata_df)
}

#' Build TetR-rtTA Dataset from protein embedding matrix
#'
#' This function creates a ProtGenesis object (Seurat-based) from a TetR-rtTA
#' protein embedding matrix file, with TetR-specific metadata extraction.
#'
#' @param embedding_file Path to embedding matrix file (csv.gz format)
#' @param assay_name Name of the assay in ProtGenesis object (default: "embed")
#' @param project_name Name of the project (default: "TetR_ProteinEmbed")
#' @param do_center Logical indicating whether to center the data (default: TRUE)
#' @param do_scale Logical indicating whether to scale the data (default: FALSE)
#' @param min.cells Minimum cells threshold for CreateSeuratObject (default: 0)
#' @param min.features Minimum features threshold for CreateSeuratObject (default: 0)
#' @param verbose Logical indicating whether to print progress messages (default: TRUE)
#'
#' @return A ProtGenesis object containing the protein embedding data
#' @export
# Build TetR dataset function | 构建 TetR 数据集函数
BuildTetRDataset <- function(embedding_file,
                                  assay_name = "embed",
                                  project_name = "TetR_ProteinEmbed",
                                  do_center = TRUE,
                                  do_scale = FALSE,
                                  min.cells = 0,
                                  min.features = 0,
                                  verbose = TRUE) {

  if (!file.exists(embedding_file)) {
    stop(sprintf("Embedding file not found: %s", embedding_file))
  }

  raw_matrix = fread(embedding_file) %>%
    as.data.frame()

  protein_names = raw_matrix[[1]]
  raw_matrix = raw_matrix[, -1]
  dat = as.matrix(raw_matrix)
  storage.mode(dat) = "numeric"
  rownames(dat) = protein_names

  dat_matrix_t = t(dat)
  protein_names_t = colnames(dat_matrix_t)
  simple_names = paste0("Protein_", 1:length(protein_names_t))
  rownames(dat_matrix_t) = simple_names

  dat_matrix = as.matrix(dat_matrix_t)
  storage.mode(dat_matrix) = "numeric"

  sparse_matrix = Matrix(
    data = dat_matrix,
    nrow = nrow(dat_matrix),
    ncol = ncol(dat_matrix),
    dimnames = list(rownames(dat_matrix), colnames(dat_matrix)),
    sparse = TRUE
  )

  protgenesis_obj = CreateSeuratObject(
    counts = sparse_matrix,
    assay = assay_name,
    project = project_name,
    min.cells = min.cells,
    min.features = min.features
  )

  protgenesis_obj = SetAssayData(
    object = protgenesis_obj,
    assay = assay_name,
    layer = "data",
    new.data = GetAssayData(protgenesis_obj, assay = assay_name, layer = "counts")
  )

  protgenesis_obj = ScaleData(
    object = protgenesis_obj,
    assay = assay_name,
    do.center = do_center,
    do.scale = do_scale
  )

  metadata = MetadataExtract_TetR(protein_names_t, verbose = verbose)
  rownames(metadata) = simple_names

  protgenesis_obj = AddMetaData(protgenesis_obj, metadata = metadata)

  if (verbose) {
    cat("=== TetR-rtTA Dataset Complete ===\n")
    cat(sprintf("Number of proteins: %d\n", ncol(protgenesis_obj)))
    cat(sprintf("Number of embedding dimensions: %d\n", nrow(protgenesis_obj)))
  }

  return(protgenesis_obj)
}

# ==============================================================================
# Initialization
# ==============================================================================

cat("==============================================================================\n")
cat("DatasetBuild Functions Loaded Successfully!\n")
cat("==============================================================================\n")
cat("\n")
cat("Available functions:\n")
cat("  - ReduceDim()              : Perform dimensionality reduction (PCA, UMAP, tSNE)\n")
cat("  - BuildProtGenesisDataset(): Build ProtGenesis object from embedding matrix\n")
cat("  - BuildTetRDataset()      : Build TetR-rtTA dataset from embedding matrix\n")
cat("  - MetadataExtract()        : Extract metadata from sequence IDs (Round format)\n")
cat("  - MetadataExtract_GFPVariants(): Extract metadata for GFP Variants\n")
cat("  - MetadataExtract_TRIM21()  : Extract metadata for TRIM21 sequences\n")
cat("  - MetadataExtract_TetR()    : Extract metadata for TetR-rtTA sequences\n")
cat("  - ClusterFind()            : Perform clustering analysis and visualization\n")
cat("\n")
cat("Author: ProtGenesis Analysis Team\n")
cat("Date: 2026-02-24\n")
cat("Version: 1.0.0\n")
cat("By: Chuanyang Liu\n")
cat("Please cite our paper if you use this script in your research:\n")
cat("Chuanyang Liu, et al. (2026). Universal physical principles govern the deterministic genesis of protein structure.\n")
cat("BioRxiv.10.64898/2026.02.20.706798.\n")
cat("\n")
cat("For help, see ?ReduceDim, ?BuildProtGenesisDataset, ?MetadataExtract, or ?ClusterFind\n")
cat("==============================================================================\n")
