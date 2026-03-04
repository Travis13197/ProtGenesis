# ProtGenesis R Package
## ProtGenesis: Universal Physical Principles Govern the Deterministic Genesis of Protein Structure

===============================================================================

## Table of Contents

1. [Overview](#1-overview)
2. [Installation](#2-installation)
3. [Core Concepts](#3-core-concepts)
4. [Quick Start](#4-quick-start)
5. [Function Reference](#5-function-reference)
6. [Data Structures](#6-data-structures)
7. [Workflow](#7-workflow)
8. [Examples](#8-examples)
9. [Citation](#9-citation)

===============================================================================

## 1. Overview

### 1.1 Breif Introduction of ProtGenesis Concept

ProtGenesis is a unified computational framework that recasts protein genesis as a structured, deterministic navigation within a discrete structural space | ProtGenesis 是基于R的蛋白质结构空间中Protein Genesis过程统一计算框架，将蛋白质起源重新定义为离散结构空间内有结构的、确定性的导航过程。

This framework hypothesizes that protein structural space is not a stochastic continuum but a hierarchically organized physical system | 该框架假设蛋白质结构空间不是随机的连续体，而是分层组织的物理系统。

ProtGenesis quantitatively reflects deterministic biological processes of protein genesis including amino acid assembly, translation and folding processes, and evolutionary adaptation | ProtGenesis 定量反映蛋白质起源的确定性生物学过程，包括氨基酸组装、翻译和折叠过程以及进化适应。



### 1.2 Key Features

| Feature | Description |
|---------|-------------|
| Protein Language Model-based embedding dataset extraction and standardization | Supports ProstT5, ProtT5 and other protein language model embeddings based on customized python scripts(Jupyter Notebook) | 支持 ProstT5、ProtT5 等面向蛋白质Genesis的蛋白质语言模型嵌入数据获取（基于自定义Python脚本）与蛋白质生成过程的标准化分析（ProtGenesis R包） |
| Dimensionality Reduction | PCA, UMAP, tSNE and diverse visualizations (heatmap, scatterplot, 3d plot) for 1024D embedding dataset | 1024D 嵌入可视化的 PCA、UMAP、tSNE，包含2d和3d可视化 |
| Tripartite Mathematical Metrics of Spatial Feature | Local Density (ρ), Spatial Dispersion (D), Differential Embedding Distance (δ) | 三种分析蛋白质结构空间特征的数学指标：局部密度(ρ)、空间离散度(D)、差分嵌入距离(δ) |
| Visualization Modes | Assembly direction map, Genesis path map, Status transition map | 三种ProtGenesis的可视化方式：组装方向图、起源路径图、状态转换图 |
| Four-tiered Framework | Models protein genesis from amino acids to full architectures | 四层框架：从氨基酸到完整蛋白质结构建模（数据集生成到完整分析与可视化过程） |

### 1.3 Core Scientific Principles

Based on the ProtGenesis paper (Liu et al., 2026), this framework uncovers three first principles | 基于 ProtGenesis 论文（Liu et al., 2026），该框架揭示了三个蛋白质生成的第一性原理及工程应用（具体请参考ProtGenesis文章）：

1. **Principle I：Hierarchical short-range order governed by directional assembly** | 蛋白质生成过程受到方向性组装的分层短程有序控制 （氨基酸按照特定规律组装成为多肽）

2. **Principle II：Differential embedding analysis identifies critical sites for protein engineering** | 差分嵌入分析，经过数学统计关键转变点，能用于识别蛋白质工程的关键操作位点（如蛋白质工程中氨基酸替换、插入或删除的关键位点）

3. **Principle III：Phase transition analysis enables principle-driven de novo protein design** | 相变分析实现原理驱动的全新蛋白质设计（类比生物分子的间断平衡进化过程，通过分析蛋白质结构空间中的突变点，实现基于第一性原理的全新蛋白质设计）




===============================================================================

## 2. Installation

### 2.1 Install from GitHub (Recommended)

```r
# Install devtools if not already installed
install.packages("devtools")

# Install ProtGenesis from GitHub
devtools::install_github("Travis131/ProtGenesis")
```

### 2.2 Local Installation

```r
# After cloning or downloading the repository
install.packages("path/to/ProtGenesis", repos = NULL, type = "source")
```

### 2.3 Dependencies

ProtGenesis requires the following R packages | ProtGenesis 需要以下 R 包：

**CRAN Packages:**
- `ggplot2` - Data visualization | 数据可视化
- `dplyr` - Data manipulation | 数据操作
- `stringr` - String processing | 字符串处理
- `data.table` - Efficient data handling | 高效数据处理
- `plotly` - Interactive visualization | 交互式可视化
- `htmlwidgets` - Interactive HTML widgets | 交互式 HTML 部件
- `Matrix` - Sparse matrices | 稀疏矩阵
- `proxy` - Distance calculations | 距离计算
- `cluster` - Clustering analysis | 聚类分析
- `umap` - UMAP dimensionality reduction | UMAP 降维
- `doParallel` - Parallel computing | 并行计算

**Bioconductor Packages:**
- `Seurat` - Single-cell analysis framework | 单细胞数据及基本分析框架
- `Biostrings` - Biological sequence handling | 生物序列处理

### 2.4 Auto-installation

ProtGenesis automatically checks and installs missing dependencies when loaded | ProtGenesis 加载时自动检查并安装缺失的依赖项：

```r
library(ProtGenesis)
# Automatically installs missing packages
```
✨ 您也可以仅调用 ProtGenesis 的五大核心脚本，快速集成全部通用功能： | ✨ You can also directly source the five core ProtGenesis scripts to unlock all built-in utilities:


| ✨ Core scripts核心脚本 | 功能 | English Function |
|-------------|----------|------------------|
| 1.Gen_funs.r| 序列生成 | Sequence generation |
| 2.DatasetBuild.r| 数据构建、标准化、降维与聚类 | Data construction, normalization, dimensionality reduction & clustering |
| 3.Visualization.r | 多模态可视化 | Multi-modal visualization |
| 4.DifferentialEmbeddingCal.r | 差分嵌入分析 | Differential-embedding analysis |
| 5.SpatialAnalysis.r | 空间特征统计 | Spatial-feature statistics |

===============================================================================


## 3. Core Concepts

### 3.1 Protein Language Models and Embedding Generation

The ProtGenesis framework uses state-of-the-art protein language models to generate high-dimensional structural embeddings. | ProtGenesis 框架使用先进的蛋白质语言模型生成高维结构嵌入。

**ProstT5** (Protein structure-sequence T5) is a bilingual language model optimized for structure-aware embedding generation. ProstT5 finetunes ProtT5-XL-U50 on translating between protein sequence and structure using 17M proteins with high-quality 3D structure predictions from the AlphaFoldDB. Unlike standard sequence-only models, ProstT5 captures both evolutionary semantics and implicit structural information, allowing mapping sequence space directly onto a continuous structural space.  
ProstT5 (Protein structure-sequence T5) 是针对结构感知嵌入生成优化的双语语言模型，基于 ProtT5-XL-U50，使用 17M 蛋白质与 AlphaFoldDB 高质量 3D 结构预测进行微调。与标准仅序列模型不同，ProstT5 同时捕获进化语义和隐式结构信息，允许将序列空间直接映射到连续的结构空间。

For each protein sequence, ProtGenesis extracts a 1024-dimensional vector representation from the model's latent space (encoder output). | 对于每个蛋白质序列，ProtGenesis 从模型的潜在空间（编码器输出）中提取 1024 维向量表示。



✨ We recommend using the 1024-dimensional ProstT5 embeddings showcased in our paper. For comparative studies, you may also employ ProtT5 (1024D) or a 2048-dimensional chimera embedding that concatenates ProtT5 and ProstT5. | ✨我们建议使用论文中展示的 1024 维 ProstT5 嵌入。如需对比研究，亦可选用 ProtT5（1024D）或将 ProtT5 与 ProstT5 拼接而成的 2048 维嵌合嵌入，我们论文中发现其实 1024 维ProstT5在某些应用场景中效果足够好。

**Installation** (see mheinzinger/ProstT5 for details): | **安装** (详见 mheinzinger/ProstT5):

```bash
pip install torch
pip install transformers
pip install sentencepiece
pip install protobuf
```

**Jupyter Notebooks for Embedding Extraction:** | **用于嵌入提取的 Jupyter Notebook:**

We provide custom Jupyter notebooks for embedding extraction. These notebooks can be found in the project repository. | 我们提供自定义的 Jupyter notebook 用于嵌入提取。这些 notebook 可以在项目仓库中找到。

| Notebook | Description |
|----------|-------------|
| `ProstT5_embedding_extract.ipynb` | Extract ProstT5 embeddings from protein sequences | 提取 ProstT5 嵌入 |
| `ProtT5_embedding_extract.ipynb` | Extract ProtT5 embeddings from protein sequences | 提取 ProtT5 嵌入 |

✨ ✨ For detailed ProtGenesis workflow scripts for principle discovery and protein engineering applications, please refer to the scripts in the Travis131/ProtGenesis_Workflow directory. | ✨ ✨ 具体的 ProtGenesis 规律发现与工程应用的 Workflow 请参考：Travis131/ProtGenesis_Workflow 目录下的脚本文件

### 3.2 Three Novel Metrics

The framework introduces three intrinsic metrics for quantifying protein genesis trajectory in the structural space | 该框架引入三个内在指标来量化蛋白质结构空间中的起源Genesis轨迹：

| Metric | Symbol | Description |
|--------|--------|-------------|
| Local Density | ρ (rho) | Measures local neighborhood density in embedding space | 测量嵌入空间中局部邻域密度 |
| Spatial Dispersion | D | Measures overall spread of protein variants | 测量蛋白质变体的整体扩散 |
| Differential Embedding Distance | δ (delta) | Measures phase transition between structural states | 测量结构状态之间的相变 |


### 3.3 Three Visualization Modalities

| Visualization | Description |
|--------------|-------------|
| Assembly Direction Map | Highlights key organizational direction of genesis | 突出显示起源的关键组织方向 |
| Genesis Path Map | Shows possible trajectories of genesis process | 显示起源过程的可能轨迹 |
| Status Transition Map | Displays mutation-driven shifts of structural status | 显示突变驱动的结构状态变化 |

### 3.4 Four-tiered Framework

ProtGenesis models protein genesis from amino acids to full architectures through four tiers | ProtGenesis 通过四层从氨基酸到完整结构建模蛋白质起源：

1. **Amino Acid Assembly** - Primitive rules of protein genesis | 氨基酸组装 - 蛋白质起源的原始规则
2. **Peptide Formation** - Short-range order establishment | 肽形成 - 短程有序建立
3. **Topology/Domain Development** - Intermediate structural motifs | 拓扑/域开发 - 中间结构 motifs
4. **Full Protein Architecture** - Complete folded structure | 完整蛋白质结构 - 完整折叠结构

===============================================================================

## 4. Quick Start

### 4.1 Basic Analysis Pipeline

```r
# Step 1: Load the ProtGenesis package | 加载 ProtGenesis 包
library(ProtGenesis)

# Step 2: Build dataset from embedding matrix | 从蛋白质语言模型的嵌入矩阵构建数据集
# Input: Protein embeddings from ProstT5/ProtT5 (1024-dimensional vectors)
protgenesis_obj <- BuildProtGenesisDataset(
  embedding_file = "protein_embeddings.csv.gz",
  project_name = "MyProject"
)

# Step 3: Perform dimensionality reduction | 执行降维分析
# Supports PCA, UMAP, and tSNE for 1024D embedding visualization
protgenesis_obj <- ReduceDim(
  protgenesis_obj,
  pca_npcs = 50,
  umap_dist = 0.3
)

# Step 4: Clustering analysis | 聚类分析
cluster_result <- ClusterFind(
  protgenesis_obj,
  dims = 1:50,
  resolution = 3
)

# Step 5: Customized Visualization | 各类别可视化，可自定义
# Generate structural map with Assembly direction map
StructuralMapVis(
  protgenesis_obj,
  group_by = "seurat_clusters"
)

# Step 6: Spatial analysis | 空间特征统计分析
# Calculate Tripartite metrics: Local Density, Spatial Dispersion, Differential Distance
spatial_result <- protein_analysis(
  data_matrix = Embedding(protgenesis_obj, "pca")
)
```

### 4.2 Sequence Generation

```r
# Generate single-point mutation sequences | 生成单点突变序列
DataGen_Mutation(
  fasta_path = "protein.fasta",
  output_dir = "output/"
)

# Generate stepwise/truncated sequences | 生成渐进式截短序列
DataGen_Stepwise(
  fasta_path = "protein.fasta",
  direction = TRUE,
  output_dir = "output/"
)
```

### 4.3 Differential Embedding Analysis

```r
# Calculate differential embedding for identifying critical mutation sites
# This reveals phase transition boundaries in structural space
differential_result <- CalculateDifferentialEmbedding(
  protgenesis_obj,
  target_sequences = c("seq1", "seq2"),
  group_by = "origident"
)

# Visualize differential results | 可视化差分结果
VisualizeDifferentialHeatmap(
  differential_result,
  output_dir = "figures/"
)
```

==============================================================================

## 5. Function Reference

### 5.1 Sequence Generation Module (1.Gen_funs.r)

#### DataGen_Mutation()

**Description**: Generates comprehensive single-site mutation sequences from input protein sequences | 从输入蛋白序列生成全面单点突变序列

**Parameters**:
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `fasta_path` | character | - | Input FASTA file path |
| `segments` | list | NULL | Optional segment specification |
| `direction` | logical | TRUE | TRUE=N→C, FALSE=C→N |
| `output_dir` | character | NULL | Output directory |
| `output_prefix` | character | "DataGen" | Output file prefix |
| `include_original` | logical | TRUE | Include original sequences |

**Returns**: Named vector of all generated sequences (invisibly) | 返回：所有生成序列的命名向量（invisibly）

**Example**:
```r
# Generate mutations for full-length sequences in forward direction
DataGen_Mutation("protein.fasta")

# Generate mutations for specific segments
DataGen_Mutation(
  "protein.fasta",
  segments = list(1:30),
  direction = FALSE,
  output_dir = "mutations/"
)
```

#### DataGen_Stepwise()

**Description**: Generates stepwise/truncated sequences mimicking translation-folding process | 生成模拟翻译-折叠过程的渐进式组装/截短序列

**Parameters**:
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `fasta_path` | character | - | Input FASTA file path |
| `segments` | list | NULL | Optional segment specification |
| `direction` | logical | TRUE | TRUE=N→C, FALSE=C→N |
| `include_all_residues` | logical | TRUE | Include all residues at each step |

**Returns**: Named vector of all generated sequences | 返回：所有生成序列的命名向量

### 5.2 Data Construction Module (2.DatasetBuild.r)

#### BuildProtGenesisDataset()

**Description**: Creates ProtGenesis object (Seurat-based) from protein embedding matrix | 从蛋白嵌入矩阵创建 ProtGenesis 对象（基于 Seurat）

**Parameters**:
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `embedding_file` | character | - | Embedding matrix file (csv.gz format) |
| `assay_name` | character | "embed" | Assay name in object |
| `project_name` | character | "ProteinEmbed" | Project name |
| `do_center` | logical | TRUE | Center the data (mean-centering) |
| `do_scale` | logical | FALSE | Scale the data |

**Note**: The input embedding matrix should contain 1024-dimensional vectors from ProstT5 or similar protein language models | 注意：输入嵌入矩阵应包含来自 ProstT5 或类似蛋白质语言模型的 1024 维向量

**Returns**: ProtGenesis object containing protein embedding data | 返回：包含蛋白嵌入数据的 ProtGenesis 对象

#### ReduceDim()

**Description**: Performs comprehensive dimensionality reduction (PCA, UMAP, tSNE) | 执行综合降维分析（PCA、UMAP、tSNE）

**Parameters**:
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `protgenesis_obj` | object | - | ProtGenesis object |
| `pca_npcs` | numeric | 50 | Number of PCA components |
| `variance_cutoff` | numeric | 0.95 | Cumulative variance cutoff |
| `umap_dist` | numeric | 0.3 | UMAP minimum distance |
| `umap_spread` | numeric | 1.5 | UMAP spread parameter |
| `tsne_perplexity` | numeric | 100 | tSNE perplexity |

**Returns**: ProtGenesis object with dimensionality reductions added | 返回：添加了降维结果的 ProtGenesis 对象

**Example**:
```r
# Standard dimensionality reduction
protgenesis_obj <- ReduceDim(protgenesis_obj)

# Custom parameters for large datasets
protgenesis_obj <- ReduceDim(
  protgenesis_obj,
  pca_npcs = 50,
  umap_dist = 0.3,
  umap_neighbors = 50
)
```

#### ClusterFind()

**Description**: Performs clustering analysis including FindNeighbors and FindClusters | 执行包括 FindNeighbors 和 FindClusters 的聚类分析

**Parameters**:
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `protgenesis_obj` | object | - | ProtGenesis object |
| `dims` | numeric | 1:50 | PCA dimension range |
| `k_param` | numeric | 50 | Number of neighbors |
| `resolution` | numeric | 3 | Clustering resolution |
| `run_visualization` | logical | TRUE | Generate visualization |

**Returns**: Clustering results list | 返回：聚类结果列表

#### MetadataExtract()

**Description**: Extracts comprehensive metadata from sequence IDs | 从序列 ID 中提取全面元数据

**Parameters**:
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `ids` | character | - | Vector of sequence IDs |
| `verbose` | logical | TRUE | Print progress messages |

**Returns**: Data frame with columns: Protein_Name, Seq_Index, Seg_Type, Position, Amino_Acid, Direction, Seg_Start, Seg_End | 返回：包含以下列的数据框：Protein_Name, Seq_Index, Seg_Type, Position, Amino_Acid, Direction, Seg_Start, Seg_End

### 5.3 Visualization Module (3.Visualization.r)

#### StructuralMapVis()

**Description**: Generates comprehensive structural map visualization | 生成综合结构图可视化

**Parameters**:
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `protgenesis_obj` | object | - | ProtGenesis object |
| `group_by` | character | NULL | Grouping variable |
| `reduction` | character | "tsne" | Reduction method |
| `pt_size` | numeric | 0.02 | Point size |
| `output_dir` | character | "Figures" | Output directory |

**Returns**: Visualization results list including Assembly Direction Map | 返回：可视化结果列表，包括组装方向图

#### AssemblyDirectionMap()

**Description**: Generates assembly direction map highlighting key organizational directions | 生成组装方向图，突出显示关键组织方向

**Parameters**:
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `differential_df` | data.frame | - | Differential embedding data |
| `target_sequences` | character | - | Target sequences |
| `output_dir` | character | "." | Output directory |

#### VisualizeDifferentialHeatmap()

**Description**: Visualizes differential embedding as heatmap | 将差分嵌入可视化为热图

**Parameters**:
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `differential_df` | data.frame | - | Differential embedding data |
| `output_dir` | character | "." | Output directory |

#### VisualizePositionBoxplots()

**Description**: Generates boxplots for position-specific analysis | 生成位置特异性分析的箱线图

**Parameters**:
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `differential_df` | data.frame | - | Differential embedding data |
| `positions_to_visualize` | numeric | c(1,2,3) | Positions to visualize |
| `output_dir` | character | "." | Output directory |

### 5.4 Spatial Analysis Module (4.SpatialAnalysis.r)

#### protein_analysis()

**Description**: Performs comprehensive protein embedding space analysis including Tripartite metrics | 执行综合蛋白嵌入空间分析，包括三大指标

**Parameters**:
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `data_matrix` | matrix | - | Embedding matrix (1024D from ProstT5) |
| `method` | character | "cosine" | Distance calculation method |
| `bootstrap` | logical | TRUE | Perform bootstrap analysis |
| `parallel` | logical | FALSE | Enable parallel computing |
| `workers` | numeric | 4 | Number of parallel workers |
| `fast_mode` | logical | FALSE | Fast computation mode |

**Returns**: List containing:
- `dispersion` - Spatial Dispersion (D) metric | 空间离散度(D)指标
- `density` - Local Density (ρ) metric | 局部密度(ρ)指标
- `cluster_significance` - Clustering significance | 聚类显著性
- `consecutive_distance` - Consecutive distance between variants | 变体之间的连续距离
- `euclidean_distance` - Euclidean distance matrix | 欧氏距离矩阵
- `plots` - Visualization plots | 可视化图

**Example**:
```r
# Calculate Tripartite metrics for protein embedding space
spatial_result <- protein_analysis(
  data_matrix = Embedding(protgenesis_obj, "pca")[, 1:50],
  method = "cosine",
  bootstrap = TRUE,
  parallel = TRUE,
  workers = 4
)

# Access individual metrics
rho <- spatial_result$density      # Local Density
D <- spatial_result$dispersion    # Spatial Dispersion
```

### 5.5 Differential Analysis Module (5.DifferentialEmbeddingCal.r)

#### CalculateDifferentialEmbedding()

**Description**: Calculates differential embedding vectors for identifying critical mutation sites | 计算差分嵌入向量以识别关键突变位点

This function reveals phase transition boundaries in structural space | 此功能揭示结构空间中的相变边界

**Parameters**:
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `protgenesis_obj` | object | - | ProtGenesis object |
| `target_sequences` | character | - | Target sequence IDs |
| `embedding_type` | character | "pca" | Embedding type |
| `group_by` | character | "ident" | Grouping variable |

**Returns**: Differential analysis results list | 返回：差分分析结果列表

===============================================================================

## 6. Data Structures

### 6.1 ProtGenesis Object

The ProtGenesis object is built on Seurat with the following components | ProtGenesis 对象基于 Seurat 构建，包含以下组件：

```
ProtGenesis Object
├── assay: "embed"
├── embeddings
│   ├── raw (1024D embeddings from ProstT5/ProtT5)
│   ├── pca (PCA reduced dimensions)
│   ├── umap (UMAP projection)
│   └── tsne (tSNE projection)
├── metadata (元数据)
│   ├── origident (original sequence ID)
│   ├── seq_type (ORIG/MUT/STEP)
│   ├── position (mutation position)
│   ├── amino_acid (mutated amino acid)
│   └── ... (custom metadata columns)
└── graphs (SNN graph for clustering)
```

### 6.2 Sequence ID Format

Standard sequence ID format | 标准序列 ID 格式：

```
{ProteinName}_Seq{Index}_{Type}_{Position}_{AA}_{Direction}_Seg{Start}-{End}
```

**Examples**:
- `GFP_Seq1_ORIG_Fwd_Seg1-238` - Original GFP sequence
- `GFP_Seq1_MUT_Pos65_V_Fwd_Seg1-238` - Position 65 mutated to V
- `GFP_Seq1_STEP_Pos50_K_Fwd_Seg1-238` - Stepwise at position 50

### 6.3 Embedding Data Format

Input embedding files should contain | 输入嵌入文件应包含：

- Rows: Individual protein sequences | 行：单个蛋白质序列
- Columns: 1024-dimensional embedding vectors | 列：1024 维嵌入向量

Supported file formats: `.csv.gz`, `.csv`, `.tsv.gz` | 支持文件格式：.csv.gz, .csv, .tsv.gz

===============================================================================

## 7. Workflow

### 7.1 Standard Analysis Workflow

```r
# Step 1: Set working directory | 设置工作目录
setwd("your/project/directory")

# Step 2: Load ProtGenesis package | 加载 ProtGenesis 包
library(ProtGenesis)

# Step 3: Build dataset from protein embeddings
# Input: 1024D embeddings from ProstT5 or similar protein language model
protgenesis_obj <- BuildProtGenesisDataset(
  embedding_file = "protein_embeddings.csv.gz"
)

# Step 4: Extract metadata (if applicable) | 提取元数据
metadata <- MetadataExtract(Idents(protgenesis_obj))
protgenesis_obj <- AddMetaData(protgenesis_obj, metadata)

# Step 5: Dimensionality reduction | 降维分析
# PCA for variance extraction, UMAP for topology, tSNE for clustering
protgenesis_obj <- ReduceDim(protgenesis_obj)

# Step 6: Clustering analysis | 聚类分析
cluster_result <- ClusterFind(protgenesis_obj)

# Step 7: Visualization | 可视化
# Generate Assembly Direction Map, Genesis Path Map
StructuralMapVis(
  protgenesis_obj,
  group_by = "seurat_clusters"
)

# Step 8: Spatial analysis with Tripartite metrics
# Calculate Local Density (ρ), Spatial Dispersion (D), Differential Distance (δ)
spatial_result <- protein_analysis(
  data_matrix = Embedding(protgenesis_obj, "pca")
)

# Step 9: Save results | 保存结果
saveRDS(protgenesis_obj, "protgenesis_analysis.rds")
```

### 7.2 Differential Analysis Workflow

```r
# Step 1: Load previous analysis results | 加载之前的分析结果
protgenesis_obj <- readRDS("protgenesis_analysis.rds")

# Step 2: Select target sequences for differential analysis
# Identify critical mutation sites through differential embedding
target_seqs <- c(
  "GFP_Seq1_MUT_Pos65_V_Fwd_Seg1-238",
  "GFP_Seq1_MUT_Pos65_A_Fwd_Seg1-238"
)

# Step 3: Calculate differential embedding
# Reveals phase transition between structural states
diff_result <- CalculateDifferentialEmbedding(
  protgenesis_obj,
  target_sequences = target_seqs,
  group_by = "amino_acid"
)

# Step 4: Visualize differential results | 可视化差分结果
# Heatmap showing differential embedding patterns
VisualizeDifferentialHeatmap(
  diff_result,
  output_dir = "differential_figures/"
)

# Boxplots for position-specific analysis
VisualizePositionBoxplots(
  diff_result,
  positions_to_visualize = c(1, 2, 3),
  output_dir = "differential_figures/"
)
```

===============================================================================

## 8. Examples

### 8.1 GFP Variant Analysis

```r
# GFP protein variant analysis | GFP 蛋白变体分析
library(ProtGenesis)

# Build GFP dataset from ProstT5 embeddings
gfp_obj <- BuildProtGenesisDataset(
  embedding_file = "GFP_embeddings.csv.gz",
  project_name = "GFP_Analysis"
)

# Dimensionality reduction
gfp_obj <- ReduceDim(gfp_obj, pca_npcs = 50)

# Clustering
gfp_clusters <- ClusterFind(gfp_obj, resolution = 2.5)

# Visualization with Assembly Direction Map
StructuralMapVis(
  gfp_obj,
  group_by = "origident",
  reduction = "umap"
)

# Spatial analysis - Calculate Tripartite metrics
gfp_spatial <- protein_analysis(
  Embedding(gfp_obj, "pca")[, 1:50],
  method = "cosine"
)

# Access metrics
cat("Local Density (rho):", gfp_spatial$density, "\n")
cat("Spatial Dispersion (D):", gfp_spatial$dispersion, "\n")
```

### 8.2 TRIM21 Variant Analysis

```r
# TRIM21 protein analysis | TRIM21 蛋白分析
trim21_obj <- BuildProtGenesisDataset(
  embedding_file = "TRIM21_embeddings.csv.gz",
  project_name = "TRIM21"
)

trim21_obj <- ReduceDim(trim21_obj)
trim21_clusters <- ClusterFind(trim21_obj)

# Differential analysis for critical site identification
diff_trim21 <- CalculateDifferentialEmbedding(
  trim21_obj,
  target_sequences = NULL,
  group_by = "position"
)

# Visualize differential patterns
VisualizeScatterplots(
  diff_trim21,
  target_sequences = unique(diff_trim21$group),
  output_dir = "TRIM21_differential/"
)
```

===============================================================================

## 9. Citation

If you use ProtGenesis in your research, please cite | 如果您在研究中使用 ProtGenesis，请引用：

```
Chuanyang Liu, et al. (2026). Universal physical principles govern the deterministic genesis of protein structure.
BioRxiv. https://doi.org/10.64898/2026.02.20.706798
```

### Corresponding Paper Abstract

This framework introduces a unified methodological framework ProtGenesis that recasts genesis of protein as a structured, deterministic navigation within a discrete structural space | 该框架引入统一的方法学框架 ProtGenesis，将蛋白质起源重新定义为离散结构空间内有结构的、确定性导航。

ProtGenesis identifies three universal principles governing this hierarchical organization: the Assembly Principle directs amino acids condensation into multilayer fractal-like architectures; the Emergence Principle ensures nascent peptides’ emergence follow deterministic spatial trajectories; and the Phase-Transition Principle describes wherein incremental residue accrual or mutations drives precise topological phase shifts from short-range to long-range order | ProtGenesis 揭示支配这种层级组织的三大普适原理：组装原理指导氨基酸缩合形成多层分形结构；涌现原理确保新生肽段沿确定性空间轨迹出现；相变原理阐明残基逐步累积或突变如何驱动从短程到长程有序的精确拓扑相移。

By quantifying these trajectories with novel tripartite spatial metrics, ProtGenesis provides an universal interpretable mathematical foundation for decoding “black-box” of deep learning models and establishes a rigorous basis for exploring, understanding, and engineering the molecular blueprint of life | 凭借三项新型空间指标对这些轨迹进行量化，ProtGenesis 为解码深度学习模型的“黑箱”提供了可解释的普适数学基础，并为探索、理解和工程化生命的分子蓝图奠定了严谨根基。

===============================================================================

## Appendix

### A. Package Information

```r
# Check package version
packageVersion("ProtGenesis")

# View all exported functions
ls("package:ProtGenesis")
```

### B. Contact

- Authors: ProtGenesis Analysis Team
- Maintainer: Chuanyang Liu
- Email: liuchuanyang13@email.com

### C. License

GPL-3

===============================================================================

Document Version: 1.0.0
Last Updated: 2026-03-03

Based on: Liu et al. (2026). Universal physical principles govern the deterministic genesis of protein structure. BioRxiv.
===============================================================================
