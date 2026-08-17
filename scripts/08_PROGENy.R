# scRNA-Seq. VPS 6 weeks - Immune Cells
# Script 08: PROGENy pathway activity
# written by Gideon JL Schaefer
# 2026-08-14
# ============================================================================

### define paths
setwd("your_wd")
input_rds <- "data/rds/annotated_RDS.rds"

### activate venv
source("renv/activate.R")

# load libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(progeny)
library(decoupleR)
library(pheatmap)
library(tidyr)
library(writexl)
library(SeuratObject)
library(Matrix)
library(scales)
library(patchwork)

# ============================================================================
# Load object
# ============================================================================

data <- readRDS(input_rds)
DefaultAssay(data) <- "RNA"
Idents(data) <- data$cluster_annotation

# Metadata from script 02: HTO_classification = Mouse1-Mouse11; intervention = VPS/Sham
data$intervention <- factor(as.character(data$intervention), levels = c("Sham", "VPS"))

# ============================================================================
# PROGENy (Mouse) 
# ============================================================================

# extract count-matrix
expr_data <- LayerData(data, assay = "RNA", layer = "data")

# PROGENy model: wide gene (rownames) x pathway (columns) weights
prog_model <- progeny::getModel(organism = "Mouse", top = 500)
genes_present <- intersect(rownames(expr_data), rownames(prog_model))
prog_model <- prog_model[genes_present, , drop = FALSE]

W <- Matrix(as.matrix(prog_model), sparse = TRUE)

# pathway x cells (ensure rowname alignment)
expr_sub <- expr_data[rownames(W), , drop = FALSE]
progeny_scores <- t(W) %*% expr_sub
progeny_scores <- as(progeny_scores, "dgCMatrix")

data[["PROGENy"]] <- CreateAssayObject(data = progeny_scores)
DefaultAssay(data) <- "PROGENy"
data <- ScaleData(data, assay = "PROGENy", verbose = TRUE)

# ============================================================================
# Calculate difference in progeny activity (VPS - Sham)
# ============================================================================

calc_delta <- function(seurat_obj, assay_name) {
  mat <- LayerData(seurat_obj, assay = assay_name, layer = "scale.data")
  md <- seurat_obj@meta.data
  out <- list()
  for (cl in unique(as.character(md$cluster_annotation))) {
    cells <- rownames(md)[as.character(md$cluster_annotation) == cl]
    c1 <- cells[as.character(md[cells, "intervention"]) == "VPS"]
    c2 <- cells[as.character(md[cells, "intervention"]) == "Sham"]
    c1 <- intersect(c1, colnames(mat))
    c2 <- intersect(c2, colnames(mat))
    if (!length(c1) || !length(c2)) next
    out[[cl]] <- rowMeans(mat[, c1, drop = FALSE], na.rm = TRUE) - rowMeans(mat[, c2, drop = FALSE], na.rm = TRUE)
  }
  m <- do.call(cbind, out)
  colnames(m) <- names(out)
}

prog_delta <- calc_delta(data, "PROGENy")




