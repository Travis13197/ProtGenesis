# ==============================================================================
# 4.SpatialAnalysis.r | 空间分析函数
# ==============================================================================
#
# Purpose: Spatial analysis functions for protein embedding data
#   This script provides comprehensive functions for analyzing protein embedding
#   spaces, including density analysis, clustering, and visualization.
#
# Author: ProtGenesis Analysis Team
# Date: 2026-02-24
# Version: 1.0.0
#
# Usage:
#   source("4.SpatialAnalysis.R")
#   result = protein_analysis(data_matrix)
#
# Requirements:
#   - ggplot2 (for visualization)
#   - dplyr (for data manipulation)
#   - proxy (for distance calculations)
#   - cluster (for clustering analysis)
#   - umap (for dimensionality reduction)
#   - doParallel (for parallel processing)
#
# ==============================================================================

# ==============================================================================
# Package Management | 包管理
# ==============================================================================

check_and_install_packages <- function() {
  required_packages <- c("ggplot2", "dplyr", "proxy", "cluster", "umap", "doParallel")

  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, repos = "https://cloud.r-project.org")
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

#' Validate protein embedding data
#'
#' This function validates protein embedding data to ensure it meets
#' the required format and quality standards.
#'
#' @param data_matrix Numeric matrix of protein embeddings (samples × features)
#'
#' @return NULL (throws error if validation fails)
#' @noRd
validate_embedding <- function(data_matrix) {
  if(ncol(data_matrix)!=1024) stop("需1024维ProstT5嵌入向量")
  
  # 同时检查中心化和标准化
  col_means = colMeans(data_matrix)
  col_sd = apply(data_matrix, 2, sd)
  
  if(any(abs(col_means) > 0.1)) 
    warning("检测到非中心化嵌入值 (平均偏移量 > 0.1)")
  if(any(abs(col_sd - 1) > 0.2)) 
    warning("检测到非标准化嵌入值 (标准差偏离1 > 0.2)")
}

#' Calculate dispersion of protein embeddings
#'
#' This function calculates the dispersion of protein embeddings
#' using the specified distance method.
#'
#' @param data_matrix Numeric matrix of protein embeddings
#' @param method Distance method to use (default: "cosine")
#'
#' @return Numeric value representing dispersion
#' @noRd
calculate_dispersion <- function(data_matrix, method="cosine") {
  dist_matrix = as.matrix(dist(data_matrix, method=method, diag=TRUE))
  mean(apply(dist_matrix, 1, sd))
}

#' Calculate density of protein embeddings
#'
#' This function calculates the density of protein embeddings
#' using k-nearest neighbors approach with dynamic k value.
#'
#' @param data_matrix Numeric matrix of protein embeddings
#' @param method Distance method to use (default: "cosine")
#'
#' @return Numeric value representing mean density
#' @noRd
calculate_density <- function(data_matrix, method="cosine") {
  n_samples = nrow(data_matrix)
  k = max(5, floor(0.01 * n_samples))  # 自适应k值
  
  dist_matrix = as.matrix(dist(data_matrix, method=method, diag=TRUE))
  diag(dist_matrix) = NA
  
  density_values = apply(dist_matrix, 1, function(x) {
    sorted = sort(x, na.last=NA)
    if(length(sorted) < k) return(NA)
    1/mean(sorted[1:k])
  })
  
  mean(density_values, na.rm=TRUE)
}

#' Calculate consecutive distance between protein embeddings
#'
#' This function calculates the cosine distance between consecutive
#' protein sequences in the embedding space.
#'
#' @param embedding_matrix Numeric matrix of protein embeddings
#'
#' @return Named vector of consecutive distances
#' @export
#'
#' @examples
#' \dontrun{
#' # Calculate consecutive distances
#' distances = calculate_consecutive_distance(embedding_matrix)
#' print(distances)
#' }
calculate_consecutive_distance <- function(embedding_matrix) {
  # 验证输入矩阵
  if(nrow(embedding_matrix) < 2) stop("嵌入矩阵至少需要包含2个蛋白质序列")
  if(ncol(embedding_matrix)!=1024) stop("需1024维ProstT5嵌入向量")

  # 确保输入是矩阵格式
  if(!is.matrix(embedding_matrix)) {
    embedding_matrix = as.matrix(embedding_matrix)
  }

  # 计算连续样本间的余弦距离
  distances = numeric(nrow(embedding_matrix) - 1)
  for (i in 1:(nrow(embedding_matrix)-1)) {
    x = embedding_matrix[i, ]
    y = embedding_matrix[i+1, ]
    # 计算余弦相似度
    cos_sim = sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2)))
    # 转换为余弦距离
    distances[i] = 1 - cos_sim
  }

  # 返回包含连续距离的命名向量
  names(distances) = paste0("N", 1:(nrow(embedding_matrix)-1), "-N", 2:nrow(embedding_matrix))
  distances
}

#' Calculate Euclidean distance between protein embeddings
#'
#' This function calculates the Euclidean distance between consecutive
#' protein sequences in the embedding space.
#'
#' @param embedding_matrix Numeric matrix of protein embeddings
#'
#' @return Named vector of Euclidean distances
#' @export
#'
#' @examples
#' \dontrun{
#' # Calculate Euclidean distances
#' distances = calculate_euclidean_distance(embedding_matrix)
#' print(distances)
#' }
calculate_euclidean_distance <- function(embedding_matrix) {
  # 验证输入矩阵
  if(nrow(embedding_matrix) < 2) stop("嵌入矩阵至少需要包含2个蛋白质序列")
  if(ncol(embedding_matrix)!=1024) stop("需1024维ProstT5嵌入向量")

  # 计算连续样本间的欧氏距离
  distances = sapply(1:(nrow(embedding_matrix)-1), function(i) {
    x = embedding_matrix[i, ]
    y = embedding_matrix[i+1, ]
    as.numeric(proxy::dist(x, y, method = "euclidean"))
  })

  # 返回包含连续距离的命名向量
  names(distances) = paste0("N", 1:(nrow(embedding_matrix)-1), "-N", 2:nrow(embedding_matrix))
  distances
}

#' Calculate cluster significance
#'
#' This function calculates cluster significance using gap statistic
#' with optional parallel processing.
#'
#' @param data_matrix Numeric matrix of protein embeddings
#' @param bootstrap Logical indicating whether to use bootstrap (default: FALSE)
#' @param B Number of bootstrap samples (default: 10)
#' @param parallel Logical indicating whether to use parallel processing (default: FALSE)
#' @param workers Number of parallel workers (default: 4)
#'
#' @return Numeric value representing cluster significance
#' @noRd
calculate_cluster_significance <- function(data_matrix, bootstrap=FALSE, 
                                          B=10, parallel=FALSE, workers=4) {
  if(!bootstrap) return(NA)
  
  if(parallel) {
    cl = makeCluster(workers)
    registerDoParallel(cl)
    on.exit(stopCluster(cl))
  }
  
  # 动态调整参数，确保k_max不超过样本数
  n_samples = nrow(data_matrix)
  k_max = min(if(n_samples > 1e5) 6 else 8, max(2, floor(n_samples * 0.5)))
  
  gap_stat = clusGap(data_matrix, 
                      FUN = function(x,k) {
                        kmeans(x, k, nstart=3, iter.max=15)
                      },
                      K.max = k_max,
                      B = B,
                      verbose = FALSE)
  
  max(gap_stat$Tab[,"gap"] - gap_stat$Tab[,"SE.sim"])
}

#' Plot protein embeddings using UMAP
#'
#' This function visualizes protein embeddings using UMAP dimensionality reduction.
#'
#' @param data_matrix Numeric matrix of protein embeddings
#'
#' @return ggplot object with UMAP visualization
#' @noRd
plot_protein_umap <- function(data_matrix) {
  config = umap.defaults
  config$n_neighbors = 15
  layout = umap(data_matrix, config)
  ggplot(data.frame(layout$layout), aes(X1,X2)) + 
    geom_point(alpha=0.6, color="#2b8cbe") +
    ggtitle("UMAP蛋白质嵌入可视化")
}

#' Plot local density distribution
#'
#' This function plots the distribution of local density values
#' for protein embeddings.
#'
#' @param data_matrix Numeric matrix of protein embeddings
#'
#' @return ggplot object with density distribution
#' @noRd
plot_local_density <- function(data_matrix) {
  dist_matrix = as.matrix(dist(data_matrix, method="cosine"))
  n_samples = nrow(data_matrix)
  k = max(5, floor(0.01 * n_samples))
  density_values = apply(dist_matrix, 1, function(x) {
    sorted = sort(x, na.last=NA)
    if(length(sorted) < k) return(NA)
    1/mean(sorted[1:k])
  })
  
  ggplot(data.frame(density=density_values), aes(x=density)) +
    geom_density(fill="#2b8cbe", alpha=0.6) +
    ggtitle("局部密度分布曲线") +
    xlab("密度值 (1/平均最近邻距离)") +
    ylab("概率密度")
}

#' Comprehensive protein embedding analysis
#'
#' This function performs comprehensive analysis of protein embedding data,
#' including dispersion calculation, cluster significance testing,
#' density estimation, and visualization.
#'
#' @param data_matrix Numeric matrix of protein embeddings (samples × 1024)
#' @param method Distance method to use (default: "cosine")
#' @param bootstrap Logical indicating whether to use bootstrap for cluster significance (default: TRUE)
#' @param parallel Logical indicating whether to use parallel processing (default: FALSE)
#' @param workers Number of parallel workers (default: 4)
#' @param fast_mode Logical indicating whether to use fast mode with data scaling (default: FALSE)
#'
#' @return List containing analysis results and visualizations
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic analysis
#' result = protein_analysis(data_matrix)
#' print(result$stats)
#' 
#' # Analysis with parallel processing
#' result = protein_analysis(data_matrix, parallel = TRUE, workers = 4)
#' print(result$stats)
#' 
#' # View visualizations
#' print(result$plots$umap)
#' print(result$plots$density)
#' }
# Main protein analysis function | 主蛋白质分析函数
protein_analysis <- function(data_matrix, method="cosine", bootstrap=TRUE, 
                            parallel=FALSE, workers=4, fast_mode=FALSE) {
  # 数据预处理
  if(fast_mode) {
    data_matrix = scale(data_matrix, center=TRUE, scale=apply(data_matrix, 2, sd))
    message("启用快速模式：执行维度缩放优化")
  }
  
  # 增强验证
  validate_embedding(data_matrix)
  
  # 核心计算流程
  metrics = list(
    dispersion = calculate_dispersion(data_matrix, method),
    cluster_sig = calculate_cluster_significance(
      data_matrix,
      bootstrap = bootstrap,
      parallel = parallel,
      workers = workers
    ),
    density = calculate_density(data_matrix, method)
  )
  
  # 可视化生成
  visualization = list(
    umap = plot_protein_umap(data_matrix),
    density = plot_local_density(data_matrix)
  )
  
  # 资源清理
  if(parallel) {
    message("正在释放并行计算资源...")
    stopImplicitCluster()
    gc()
  }
  
  # 返回标准化结果
  list(stats = metrics, plots = visualization)
}

# ==============================================================================
# Initialization
# ==============================================================================

cat("==============================================================================\n")
cat("SpatialAnalysis Functions Loaded Successfully!\n")
cat("==============================================================================\n")
cat("\n")
cat("Available functions:\n")
cat("  - protein_analysis()                : Comprehensive protein embedding analysis\n")
cat("  - calculate_consecutive_distance()  : Calculate consecutive cosine distances\n")
cat("  - calculate_euclidean_distance()    : Calculate consecutive Euclidean distances\n")
cat("\n")
cat("For help, see ?protein_analysis\n")
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
cat("For help, see ?protein_analysis\n")
cat("==============================================================================\n")
