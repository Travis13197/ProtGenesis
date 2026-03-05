# ==============================================================================
# 5.DifferentialEmbeddingCal.R
# ==============================================================================
#
# Purpose: Comprehensive differential embedding calculation functions
#   This script provides a unified framework for calculating differential
#   embeddings from protein sequence data, supporting multiple embedding types
#   including raw 1024D, PCA, t-SNE, and UMAP.
#
# Author: ProtGenesis Analysis Team
# Date: 2026-02-26
# Version: 1.0.0
# By: Chuanyang Liu
# Please cite our paper if you use this script in your research:
#  Chuanyang Liu, et al. (2026). Universal physical principles govern the deterministic genesis of protein structure.
#  BioRxiv.10.64898/2026.02.20.706798.
#
# Usage:
#   source("5.DifferentialEmbeddingCal.R")
#   result = CalculateDifferentialEmbedding(
#     protgenesis_obj = protgenesis_obj,
#     target_sequences = c("G", "GG", "GGG"),
#     embedding_type = "raw"
#   )
#
# Requirements:
#   - Seurat (for single-cell analysis)
#   - Matrix (for sparse matrix operations)
#   - dplyr (for data manipulation)
#   - data.table (for fast data reading)
#   - tidyverse (for tidy data operations)
#
# ==============================================================================

# ==============================================================================
# Package Management
# ==============================================================================



# ==============================================================================
# Global Constants
# ==============================================================================

standard_amino_acids = c('A','C','D','E','F',
                         'G','H','I','K','L',
                         'M','N','P','Q','R',
                         'S','T','V','W','Y')

# ==============================================================================
# Core Functions
# ==============================================================================

#' Convert Cell ID to Amino Acid Sequence
#'
#' This function converts a numeric Cell ID to its corresponding amino acid
#' sequence using the standard amino acid mapping.
#'
#' @param cell_id Character string representing the cell ID (e.g., "1_7" for "G")
#' @return Amino acid sequence as character string, or NA if conversion fails
#' @export
#'
#' @examples
#' \dontrun{
#' get_amino_acid_sequence("1_7")  # Returns "G"
#' get_amino_acid_sequence("2_7_1")  # Returns "GA"
#' }
GetAminoAcidSequence <- function(cell_id) {
  all_parts = strsplit(cell_id, "_")[[1]]
  
  numeric_start_index = which(sapply(all_parts, function(x) !is.na(suppressWarnings(as.numeric(x)))))[1]
  
  if (!is.na(numeric_start_index)) {
    numeric_parts = all_parts[numeric_start_index:length(all_parts)]
    parts = as.numeric(numeric_parts)
  } else {
    parts = as.numeric(all_parts)
  }
  
  if (any(is.na(parts))) {
    return(NA)
  }
  
  len = parts[1]
  indices = parts[-1]
  
  if (len == 1) {
    sequence = standard_amino_acids[indices[1]]
  } else if (len == 2) {
    sequence = paste(standard_amino_acids[indices[c(1,2)]], collapse="")
  } else if (len == 3) {
    sequence = paste(standard_amino_acids[indices[c(1,2,3)]], collapse="")
  } else if (len == 4) {
    sequence = paste(standard_amino_acids[indices[c(1,2,3,4)]], collapse="")
  } else {
    sequence = NA
  }
  
  return(sequence)
}

#' Convert Amino Acid Sequence to Cell ID
#'
#' This function converts an amino acid sequence to its corresponding Cell ID
#' using the standard amino acid mapping.
#'
#' @param sequences Character vector of amino acid sequences
#' @return Character vector of corresponding Cell IDs
#' @export
#'
#' @examples
#' \dontrun{
#' SeqToCellId(c("G", "GA"))  # Returns c("1_7", "2_7_1")
#' }
SeqToCellId <- function(sequences) {
  amino_acid_codes = 1:length(standard_amino_acids)
  names(amino_acid_codes) = standard_amino_acids
  
  result = c()
  
  for(seq in sequences) {
    if(nchar(seq) == 0) {
      result = c(result, NA)
      next
    }
    
    seq = toupper(seq)
    valid_seq = all(strsplit(seq, "")[[1]] %in% standard_amino_acids)
    if(!valid_seq) {
      result = c(result, paste0("Invalid_", seq))
      next
    }
    
    seq_length = nchar(seq)
    aa_codes = sapply(strsplit(seq, "")[[1]], function(aa) amino_acid_codes[aa])
    cell_id = paste(c(seq_length, aa_codes), collapse = "_")
    result = c(result, cell_id)
  }
  
  return(result)
}

#' Extract Sequence Embeddings from Seurat Object
#'
#' This function extracts embeddings for target sequences and their next-position
#' variants from a Seurat object.
#'
#' @param sequence Target amino acid sequence to extract
#' @param seq_data Embedding matrix from Seurat object
#' @param meta_data Metadata dataframe from Seurat object
#' @param verbose Logical indicating whether to print progress information
#' @return List containing exact and next-variant embeddings
#' @export
#'
#' @examples
#' \dontrun{
#' result = ExtractSequenceEmbeddingsSingle(
#'   sequence = "G",
#'   seq_data = protgenesis_obj@assays$embed@layers$scale.data,
#'   meta_data = protgenesis_obj@meta.data
#' )
#' }
ExtractSequenceEmbeddingsSingle <- function(sequence, seq_data, meta_data, verbose = TRUE) {
  exact_match_indices = which(meta_data$Amino_Acid_Sequence == sequence)
  exact_embeddings = seq_data[, exact_match_indices, drop = FALSE]
  exact_sequences = meta_data$Amino_Acid_Sequence[exact_match_indices]
  
  next_variants = paste0(sequence, standard_amino_acids)
  next_var_indices = which(meta_data$Amino_Acid_Sequence %in% next_variants)
  next_var_embeddings = seq_data[, next_var_indices, drop = FALSE]
  next_var_sequences = meta_data$Amino_Acid_Sequence[next_var_indices]
  
  result = list(
    exact = exact_embeddings,
    exact_sequences = exact_sequences,
    next_variants = next_var_embeddings,
    next_var_sequences = next_var_sequences
  )
  
  if (verbose) {
    cat(sprintf("序列 '%s' 提取结果:\n", sequence))
    cat(sprintf("- 精确匹配: %d 个\n", length(exact_match_indices)))
    cat(sprintf("- 下一位变体: %d 个 (预期20个)\n", length(next_var_indices)))
  }
  
  return(result)
}

#' Get Embedding Data from Seurat Object
#'
#' This function extracts different types of embedding data from a Seurat object.
#'
#' @param protgenesis_obj Seurat object containing protein embeddings
#' @param embedding_type Type of embedding to extract: "raw", "pca", "tsne", or "umap"
#' @param pca_npcs Number of PCA components to use (default: 50)
#' @param verbose Logical indicating whether to print progress information
#' @return Matrix of embeddings
#' @export
#'
#' @examples
#' \dontrun{
#' raw_embeddings = GetEmbeddingData(protgenesis_obj, "raw")
#' pca_embeddings = GetEmbeddingData(protgenesis_obj, "pca", pca_npcs = 50)
#' }
GetEmbeddingData <- function(protgenesis_obj, embedding_type = c("raw", "pca", "tsne", "umap"), 
                             pca_npcs = 50, verbose = TRUE) {
  embedding_type = match.arg(embedding_type)
  
  if (verbose) {
    cat(sprintf("提取 %s 嵌入数据...\n", toupper(embedding_type)))
  }
  
  switch(embedding_type,
         "raw" = {
           if ("embed" %in% names(protgenesis_obj@assays)) {
             if ("scale.data" %in% names(protgenesis_obj@assays$embed@layers)) {
               embeddings = t(protgenesis_obj@assays$embed@layers$scale.data)
             } else if ("data" %in% names(protgenesis_obj@assays$embed@layers)) {
               embeddings = t(protgenesis_obj@assays$embed@layers$data)
             } else {
               stop("在Seurat对象中未找到嵌入数据")
             }
             if (is.null(rownames(embeddings))) {
               rownames(embeddings) = colnames(protgenesis_obj)
             }
           } else {
             stop("Seurat对象中未找到'embed' assay")
           }
         },
         "pca" = {
           if ("pca" %in% names(protgenesis_obj@reductions)) {
             npcs_use = min(pca_npcs, ncol(Embeddings(protgenesis_obj, "pca")))
             embeddings = Embeddings(protgenesis_obj, "pca")[, 1:npcs_use, drop = FALSE]
           } else {
             stop("Seurat对象中未找到PCA降维结果")
           }
         },
         "tsne" = {
           if ("tsne" %in% names(protgenesis_obj@reductions)) {
             embeddings = Embeddings(protgenesis_obj, "tsne")
           } else {
             stop("Seurat对象中未找到t-SNE降维结果")
           }
         },
         "umap" = {
           if ("umap" %in% names(protgenesis_obj@reductions)) {
             embeddings = Embeddings(protgenesis_obj, "umap")
           } else {
             stop("Seurat对象中未找到UMAP降维结果")
           }
         })
  
  if (verbose) {
    cat(sprintf("%s 嵌入数据维度: %d x %d\n", toupper(embedding_type), nrow(embeddings), ncol(embeddings)))
  }
  
  return(embeddings)
}

#' Calculate Differential Embeddings
#'
#' This is the main function for calculating differential embeddings from a
#' Seurat object. It supports multiple embedding types and provides flexible
#' configuration options.
#'
#' @param protgenesis_obj Seurat object containing protein embeddings
#' @param target_sequences Character vector of target sequences to analyze
#' @param embedding_type Type of embedding to use: "raw", "pca", "tsne", or "umap"
#' @param split_by Column name in meta.data to split data by (e.g., "Terminal_Type")
#' @param split_values Values to use for splitting (if NULL, uses all unique values)
#' @param pca_npcs Number of PCA components to use (default: 50)
#' @param add_amino_acid_sequence Logical indicating whether to add amino acid sequences to metadata
#' @param verbose Logical indicating whether to print progress information
#' @return List containing differential embedding results and metadata
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage with raw embeddings
#' result = CalculateDifferentialEmbedding(
#'   protgenesis_obj = protgenesis_obj,
#'   target_sequences = c("G", "GG", "GGG"),
#'   embedding_type = "raw"
#' )
#'
#' # With Terminal_Type splitting
#' result = CalculateDifferentialEmbedding(
#'   protgenesis_obj = protgenesis_obj,
#'   target_sequences = c("G", "GG"),
#'   embedding_type = "pca",
#'   split_by = "Terminal_Type",
#'   split_values = c("Nter", "Cter")
#' )
#' }
CalculateDifferentialEmbedding <- function(protgenesis_obj,
                                           target_sequences,
                                           embedding_type = c("raw", "pca", "tsne", "umap"),
                                           split_by = NULL,
                                           split_values = NULL,
                                           pca_npcs = 50,
                                           add_amino_acid_sequence = TRUE,
                                           verbose = TRUE) {
  
  if (verbose) {
    cat("==============================================================================\n")
    cat("CalculateDifferentialEmbedding - 差异嵌入计算\n")
    cat("==============================================================================\n\n")
  }
  
  embedding_type = match.arg(embedding_type)
  
  if (!inherits(protgenesis_obj, "Seurat")) {
    stop("protgenesis_obj 必须是一个Seurat对象")
  }
  
  if (missing(target_sequences) || length(target_sequences) == 0) {
    stop("必须提供 target_sequences 参数")
  }
  
  if (add_amino_acid_sequence && !"Amino_Acid_Sequence" %in% colnames(protgenesis_obj@meta.data)) {
    if (verbose) {
      cat("添加氨基酸序列到meta.data...\n")
    }
    protgenesis_obj@meta.data$Amino_Acid_Sequence = sapply(
      protgenesis_obj@meta.data$Cell,
      GetAminoAcidSequence
    )
  }
  
  if (!is.null(split_by)) {
    if (!split_by %in% colnames(protgenesis_obj@meta.data)) {
      stop(sprintf("在meta.data中未找到列 '%s'", split_by))
    }
    
    if (is.null(split_values)) {
      split_values = unique(protgenesis_obj@meta.data[[split_by]])
    }
    
    if (verbose) {
      cat(sprintf("按 '%s' 拆分数据，值: %s\n", split_by, paste(split_values, collapse = ", ")))
    }
  } else {
    split_values = NULL
  }
  
  embeddings = GetEmbeddingData(protgenesis_obj, embedding_type, pca_npcs, verbose)
  
  result_list = list()
  
  process_split = function(split_value = NULL) {
    if (verbose && !is.null(split_value)) {
      cat(sprintf("\n处理拆分值: %s\n", split_value))
    }
    
    meta_subset = protgenesis_obj@meta.data
    if (!is.null(split_value)) {
      meta_subset = meta_subset[meta_subset[[split_by]] == split_value, , drop = FALSE]
    }
    
    common_cells = intersect(rownames(meta_subset), rownames(embeddings))
    if (length(common_cells) == 0) {
      warning(sprintf("拆分值 '%s' 没有匹配的细胞", split_value))
      return(NULL)
    }
    
    embed_subset = embeddings[common_cells, , drop = FALSE]
    meta_subset = meta_subset[common_cells, , drop = FALSE]
    
    combined_df_list = list()
    
    for (seq in target_sequences) {
      if (verbose) {
        cat(sprintf("处理序列: %s\n", seq))
      }
      
      seq_result = ExtractSequenceEmbeddingsSingle(seq, t(embed_subset), meta_subset, verbose = FALSE)
      
      if (ncol(seq_result$exact) > 0) {
        exact_df = as.data.frame(t(seq_result$exact))
        colnames(exact_df) = paste0('dim-', 0:(ncol(exact_df)-1))
        exact_df$Sequence = seq
        exact_df$Full_Sequence = seq_result$exact_sequences
        exact_df$Type = 'Exact'
        if (!is.null(split_value)) {
          exact_df[[split_by]] = split_value
        }
        combined_df_list[[length(combined_df_list) + 1]] = exact_df
      }
      
      if (ncol(seq_result$next_variants) > 0) {
        next_var_df = as.data.frame(t(seq_result$next_variants))
        colnames(next_var_df) = paste0('dim-', 0:(ncol(next_var_df)-1))
        next_var_df$Sequence = seq
        next_var_df$Full_Sequence = seq_result$next_var_sequences
        next_var_df$Type = 'Next_Variant'
        if (!is.null(split_value)) {
          next_var_df[[split_by]] = split_value
        }
        combined_df_list[[length(combined_df_list) + 1]] = next_var_df
      }
    }
    
    if (length(combined_df_list) == 0) {
      return(NULL)
    }
    
    do.call(rbind, combined_df_list)
  }
  
  if (is.null(split_values)) {
    combined_df = process_split()
  } else {
    split_results = lapply(split_values, process_split)
    split_results = split_results[!sapply(split_results, is.null)]
    if (length(split_results) == 0) {
      stop("没有找到有效的拆分结果")
    }
    combined_df = do.call(rbind, split_results)
  }
  
  if (verbose) {
    cat(sprintf("\n合并数据框维度: %d x %d\n", nrow(combined_df), ncol(combined_df)))
  }
  
  calculate_differential = function(df) {
    unique_sequences = unique(df$Sequence)
    result_list = list()
    
    for (base_seq in unique_sequences) {
      group_data = df[df$Sequence == base_seq, ]
      base_data = group_data[group_data$Type == 'Exact', ]
      
      if (nrow(base_data) > 0) {
        dim_cols = grep('^dim-', colnames(base_data), value = TRUE)
        base_embedding = as.numeric(unlist(base_data[1, dim_cols]))
        
        for (i in 1:nrow(group_data)) {
          seq_data = group_data[i, ]
          seq_data_diff = seq_data
          seq_data_diff[, dim_cols] = as.numeric(seq_data[, dim_cols]) - base_embedding
          result_list[[length(result_list) + 1]] = seq_data_diff
        }
      }
    }
    
    if (length(result_list) == 0) {
      return(NULL)
    }
    
    do.call(rbind, result_list)
  }
  
  if (!is.null(split_by)) {
    differential_list = list()
    for (sv in split_values) {
      split_df = combined_df[combined_df[[split_by]] == sv, ]
      if (nrow(split_df) > 0) {
        diff_df = calculate_differential(split_df)
        if (!is.null(diff_df)) {
          differential_list[[length(differential_list) + 1]] = diff_df
        }
      }
    }
    differential_df = do.call(rbind, differential_list)
  } else {
    differential_df = calculate_differential(combined_df)
  }
  
  if (verbose && !is.null(differential_df)) {
    cat(sprintf("差异嵌入数据框维度: %d x %d\n", nrow(differential_df), ncol(differential_df)))
  }
  
  result = list(
    combined_df = combined_df,
    differential_df = differential_df,
    protgenesis_obj = protgenesis_obj,
    embedding_type = embedding_type,
    target_sequences = target_sequences,
    split_by = split_by,
    split_values = split_values,
    timestamp = Sys.time()
  )
  
  if (verbose) {
    cat("\n==============================================================================\n")
    cat("差异嵌入计算完成!\n")
    cat("==============================================================================\n")
  }
  
  return(result)
}

#' Save Differential Embedding Results
#'
#' This function saves the differential embedding results to Rdata and CSV files.
#'
#' @param result Result list from CalculateDifferentialEmbedding
#' @param output_prefix Prefix for output files
#' @param output_dir Directory for output files
#' @param save_csv Logical indicating whether to save CSV files
#' @param verbose Logical indicating whether to print progress information
#' @export
#'
#' @examples
#' \dontrun{
#' SaveDifferentialEmbeddingResults(
#'   result = result,
#'   output_prefix = "GFPnc",
#'   output_dir = "Rdata"
#' )
#' }
SaveDifferentialEmbeddingResults <- function(result,
                                             output_prefix = "differential_embedding",
                                             output_dir = ".",
                                             save_csv = TRUE,
                                             verbose = TRUE) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  combined_df = result$combined_df
  differential_df = result$differential_df
  protgenesis_obj = result$protgenesis_obj
  
  rdata_file = file.path(output_dir, paste0(output_prefix, "_differential_embeddings.rdata"))
  save(combined_df, differential_df, protgenesis_obj, file = rdata_file)
  
  if (verbose) {
    cat(sprintf("保存Rdata文件: %s\n", rdata_file))
  }
  
  if (save_csv) {
    if (!is.null(differential_df)) {
      csv_file = file.path(output_dir, paste0(output_prefix, "_differential_embeddings.csv"))
      write.csv(differential_df, file = csv_file, row.names = FALSE)
      if (verbose) {
        cat(sprintf("保存CSV文件: %s\n", csv_file))
      }
    }
  }
  
  invisible(rdata_file)
}

#' Calculate t-SNE for differential embeddings
#'
#' This function calculates t-SNE for differential embeddings data frame.
#'
#' @param differential_df Data frame containing differential embeddings
#' @param perplexity Perplexity parameter for t-SNE (default: 30)
#' @param theta Theta parameter for t-SNE (default: 0.5)
#' @param seed Random seed for reproducibility (default: 42)
#' @param verbose Logical indicating whether to show progress messages (default: TRUE)
#'
#' @return Data frame with added tSNE_1 and tSNE_2 columns
#' @export
#'
#' @examples
#' \dontrun{
#' differential_df_with_tsne = CalculateTSNEForDifferential(differential_df)
#' }
CalculateTSNEForDifferential <- function(differential_df,
                                           perplexity = 30,
                                           theta = 0.5,
                                           seed = 42,
                                           verbose = TRUE) {
  if (verbose) {
    cat("为差异嵌入计算t-SNE...\n")
  }
  
  library(Rtsne)
  
  dim_cols = grep('^dim-', colnames(differential_df), value = TRUE)
  diff_embeddings = differential_df[, dim_cols]
  
  unique_indices = !duplicated(diff_embeddings)
  differential_df = differential_df[unique_indices, ]
  diff_embeddings = differential_df[, dim_cols]
  
  set.seed(seed)
  tsne_result = Rtsne::Rtsne(diff_embeddings, perplexity = perplexity, theta = theta, dims = 2)
  
  differential_df$tSNE_1 = tsne_result$Y[, 1]
  differential_df$tSNE_2 = tsne_result$Y[, 2]
  
  if (verbose) {
    cat("t-SNE计算完成\n")
  }
  
  return(differential_df)
}

cat("==============================================================================\n")
cat("Differential Embedding Calculation Functions Loaded Successfully!\n")
cat("==============================================================================\n\n")
cat("Available functions:\n")
cat("  - GetAminoAcidSequence()   : Convert Cell ID to amino acid sequence\n")
cat("  - SeqToCellId()            : Convert amino acid sequence to Cell ID\n")
cat("  - ExtractSequenceEmbeddingsSingle(): Extract embeddings for target sequences\n")
cat("  - GetEmbeddingData()       : Extract embedding data from Seurat object\n")
cat("  - CalculateDifferentialEmbedding(): Main differential embedding calculation\n")
cat("  - CalculateTSNEForDifferential(): Calculate t-SNE for differential embeddings\n")
cat("  - SaveDifferentialEmbeddingResults(): Save results to files\n\n")
cat("For help, see ?GetAminoAcidSequence, ?CalculateDifferentialEmbedding, etc.\n")
cat("==============================================================================\n")
cat("\n")
cat("Author: ProtGenesis Analysis Team\n")
cat("Date: 2026-02-26\n")
cat("Version: 1.0.0\n")
cat("By: Chuanyang Liu\n")
cat("Please cite our paper if you use this script in your research:\n")
cat("Chuanyang Liu, et al. (2026). Universal physical principles govern the deterministic genesis of protein structure.\n")
cat("BioRxiv.10.64898/2026.02.20.706798.\n")
cat("\n")
cat("For help, see the function documentation above.\n")
cat("==============================================================================\n")
