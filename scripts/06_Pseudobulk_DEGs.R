# scRNA-Seq. Injury 6 weeks - Immune Cells
# Script 06: Pseudobulk DEGs (DESeq2)
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
library(dplyr)
library(tidyr)
library(writexl)
library(EnhancedVolcano)
library(ggplot2)
library(grid)

data <- readRDS(input_rds)
DefaultAssay(data) <- "RNA" 
Idents(data) <- data$cluster_annotation

# Metadata from script 02: HTO_classification = Mouse1-Mouse11; intervention = Injury/Sham
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

# ============================================================================
# Pseudobulking
# ============================================================================

####Aggregate cells per mouse (HTO_classification = Mouse1-Mouse11) into pseudobulk
data$intervention <- factor(as.character(data$intervention), levels = c("Sham", "Injury"))
data$cluster_annotation <- factor(as.character(data$cluster_annotation), levels = cluster_order)

pseudo <- AggregateExpression(
  data,
  assays = "RNA",
  return.seurat = TRUE,
  group.by = c("intervention", "HTO_classification", "cluster_annotation"),
  verbose = TRUE
)

#generate a new meta-col
pseudo$newIdent <- paste0(pseudo$cluster_annotation, "_", pseudo$intervention)
Idents(pseudo) <- pseudo$newIdent

# ============================================================================
# DESeq2 DEGs per cluster (Injury vs Sham)
# ============================================================================

clusters <- as.character(unique(pseudo$cluster_annotation))
deg_list_all <- list()

for (cl in clusters) {
  ident1 <- paste0(cl, "_Injury")
  ident2 <- paste0(cl, "_Sham")

   de <- FindMarkers(
    pseudo,
    ident.1 = ident1,
    ident.2 = ident2,
    test.use = "DESeq2",
    verbose = TRUE
  )

  de$gene <- rownames(de)
  de$cluster <- cl
  de$comparison <- "Injury_vs_Sham"

  deg_list_all[[cl]] <- de

}

deg_all <- dplyr::bind_rows(deg_list_all)


