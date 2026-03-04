# ==============================================================================
# ProtGenesis-package.R | ProtGenesis 包定义文件
# ==============================================================================
#
# Purpose: Universal Protein Structure Genesis Analysis Toolkit
#   A comprehensive toolkit for protein structure genesis analysis based on
#   the ProtGenesis framework (Liu et al., 2026).
#
# Author: ProtGenesis Analysis Team
# Maintainer: Chuanyang Liu <example@email.com>
# Version: 1.0.0
# License: GPL-3
#
# ==============================================================================

#' @keywords internal
#' @aliases ProtGenesis-package
#' @docType package
#' @name ProtGenesis-package
#' @title ProtGenesis: Universal Protein Structure Genesis Analysis Toolkit
#'
#' @description
#' ProtGenesis is a unified computational framework that recasts protein genesis
#' as a structured, deterministic navigation within a discrete structural space.
#' This package provides functions for:
#' \itemize{
#'   \item Protein sequence generation (single-point mutations, stepwise truncations)
#'   \item Building ProtGenesis objects from protein language model embeddings
#'   \item Dimensionality reduction (PCA, UMAP, tSNE)
#'   \item Clustering analysis using graph-based methods
#'   \item Spatial analysis with three novel metrics (Local Density, Spatial Dispersion,
#'         Differential Embedding Distance)
#'   \item Differential embedding analysis for critical site identification
#'   \item Multi-modal visualization (Assembly Direction Map, Genesis Path Map,
#'         Status Transition Map)
#' }
#'
#' @details
#' The ProtGenesis framework is based on the discovery that protein structural space
#' follows deterministic physical laws. By embedding protein sequences using protein
#' language models (e.g., ProstT5, ProtT5), the framework enables systematic analysis
#' of protein genesis trajectories.
#'
#' @section Main Functions:
#' \describe{
#'   \item{\code{DataGen_Mutation}}{Generate single-point mutation sequences}
#'   \item{\code{DataGen_Stepwise}}{Generate stepwise/truncated sequences}
#'   \item{\code{BuildProtGenesisDataset}}{Build ProtGenesis object from embeddings}
#'   \item{\code{ReduceDim}}{Perform dimensionality reduction}
#'   \item{\code{ClusterFind}}{Perform clustering analysis}
#'   \item{\code{StructuralMapVis}}{Generate structural visualizations}
#'   \item{\code{protein_analysis}}{Perform spatial analysis with Tripartite metrics}
#'   \item{\code{CalculateDifferentialEmbedding}}{Calculate differential embeddings}
#' }
#'
#' @section Citation:
#' If you use ProtGenesis in your research, please cite:
#' Chuanyang Liu, et al. (2026). Universal physical principles govern the
#' deterministic genesis of protein structure. BioRxiv.
#' https://doi.org/10.64898/2026.02.20.706798
#'
#' @section Acknowledgments:
#' ProtGenesis builds upon the Seurat framework and various Bioconductor packages.
#' We thank the authors of ProstT5 for making the protein language model available.
#'
#' @references
#' Liu, C. et al. (2026). Universal physical principles govern the deterministic
#' genesis of protein structure. BioRxiv.
#'
#' @examples
#' \dontrun{
#' # Load the package
#' library(ProtGenesis)
#'
#' # Build dataset from embeddings
#' prot_obj <- BuildProtGenesisDataset("protein_embeddings.csv.gz")
#'
#' # Perform dimensionality reduction
#' prot_obj <- ReduceDim(prot_obj)
#'
#' # Clustering
#' cluster_result <- ClusterFind(prot_obj)
#'
#' # Visualization
#' StructuralMapVis(prot_obj, group_by = "seurat_clusters")
#'
#' # Spatial analysis
#' spatial_result <- protein_analysis(Embedding(prot_obj, "pca"))
#' }
#'
NULL
