# scRNA-Seq. VPS 6 weeks - Immune Cells
# Script 03: Clustering
# written by Gideon JL Schaefer
# 2026-08-14
# ============================================================================

### define paths
setwd("your_wd")
input_rds  <- "data/rds/harmony_integrated_RDS.rds"
output_rds <- "data/rds/clustered_RDS.rds"

### activate venv
source("renv/activate.R")

## load libraries
library(Seurat)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(dplyr)
library(patchwork)
library(harmony)
library(writexl)

# ============================================================================
# Load integrated object
# ============================================================================

data <- readRDS(input_rds)
DefaultAssay(data) <- "RNA"

# ============================================================================
# Normalization, Scaling, Variable Features, PCA
# ============================================================================
data <- NormalizeData(data, normalization.method = "LogNormalize",scale.factor = 10000, verbose = FALSE)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 4000, verbose = FALSE)
data <- ScaleData(data, verbose = FALSE)
data <- RunPCA(data, npcs = 70, verbose = FALSE)

# ============================================================================
# UMAP reduction
# ============================================================================
data <- FindNeighbors(data, reduction = "harmony", dims = 1:20)
data <- RunUMAP(data, reduction = "harmony", dims = 1:20)
 
# ============================================================================
# Clustering
# ============================================================================
data <- FindClusters(data, resolution = 0.1, verbose = FALSE)

# Final UMAP plot
DimPlot(data, reduction = "umap", label = TRUE, label.size = 10, pt.size = 0.1)

# ============================================================================
# check QC per Cluster 
# ============================================================================
VlnPlot(data, group.by = "seurat_clusters", features = "percent.mt", pt.size = 0.1) 
VlnPlot(data, group.by = "seurat_clusters", features = "nFeature_RNA", pt.size = 0.1) 
VlnPlot(data, group.by = "seurat_clusters", features = "nCount_RNA", pt.size = 0.1)
VlnPlot(data, group.by = "seurat_clusters", features = "doublet_score_scDblFinder", pt.size = 0.1) 

# ============================================================================
# get Marker Genes (only positives)
# ============================================================================
Idents(data) <- "seurat_clusters"
all.markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)


### Exclude cluster 8+9 based on QC-plots and marker genes
data <-  subset(data, subset = c(0:7))

# ============================================================================
# Save final clustered object
# ============================================================================
saveRDS(data, file = output_rds)


