# scRNA-Seq. Injury 6 weeks - Immune Cells
# Script 04: Cluster_annotation
# written by Gideon JL Schaefer
# 2026-08-14
# ============================================================================

### define paths
setwd("your_wd")
input_rds  <- "data/rds/clustered_RDS.rds"
output_rds <- "data/rds/annotated_RDS.rds"

### activate venv
source("renv/activate.R")

## load libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(writexl)
library(openxlsx)
library(dittoSeq)
library(patchwork)
library(scales)
library(tidyr)

# ============================================================================
# Load integrated object
# ============================================================================

data <- readRDS(input_rds)
DefaultAssay(data) <- "RNA"

# 
# !! NOTE: After exclusion of cluster 8+9, harmony integration and clustering were reperformed with the same parameterts set as in script 02_QC_HarmonyIntegration and 03_clustering !!
# 

# ============================================================================
# Manual final annotation
# ============================================================================

# Define final cluster annotations
cluster_annotations <- c(
  "0" = "Granulocytes",
  "1" = "naive_CD4_TCs",
  "2" = "BCs",
  "3" = "CD8_TCs",
  "4" = "Ly6c2low_Monos",
  "5" = "Ly6c2high_Monos",
  "6" = "NKCs",
  "7" = "activated_CD4_TCs"
)

# Map cluster numbers to annotations
cluster_vec <- as.character(data$seurat_clusters)
annotation_vec <- unname(cluster_annotations[cluster_vec])

cluster_order <- c(
  "Granulocytes",
  "Ly6c2high_Monos",
  "Ly6c2low_Monos",
  "BCs",
  "naive_CD4_TCs",
  "activated_CD4_TCs",
  "CD8_TCs",
  "NKCs"
)

data@meta.data$cluster_annotation <- factor(annotation_vec, levels = cluster_order)

# Set as identity
Idents(data) <- data$cluster_annotation

# ============================================================================
# Save annotated object
# ============================================================================

saveRDS(data, file = output_rds)

