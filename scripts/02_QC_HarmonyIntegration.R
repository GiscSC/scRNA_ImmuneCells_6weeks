# scRNA-Seq. Injury 6 weeks - Immune Cells
# Script 02: Doublet Exclusion, QC Filtering, and Harmony Integration
# written by Gideon JL Schaefer
# 2026-08-14

# ============================================================================
### define paths
setwd("your_wd")
output_rds <- "data/rds/harmony_integrated_RDS.rds"

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
library(scDblFinder)
library(SingleCellExperiment)

# ============================================================================
# Sample labels + intervention mapping
# ============================================================================
# library1: Mouse1-3  -> Injury
# library2: Mouse4-5  -> Injury
# library3: Mouse6-8  -> Sham
# library4: Mouse9-11 -> Sham
# HTO_classification (Mouse1-Mouse11) is set in script 01 during demultiplexing.

library_labels <- c("library1", "library2", "library3", "library4")

mouse_to_intervention <- c(
  "Mouse1" = "Injury", "Mouse2" = "Injury", "Mouse3" = "Injury",           
  "Mouse4" = "Injury", "Mouse5" = "Injury",
  "Mouse6" = "Sham", "Mouse7" = "Sham", "Mouse8" = "Sham",
  "Mouse9" = "Sham", "Mouse10" = "Sham", "Mouse11" = "Sham"
)

# ============================================================================
# Load + annotate demultiplexed RDS files
# ============================================================================
sample1 <- readRDS("data/rds/library1_demultiplexed.rds")
sample2 <- readRDS("data/rds/library2_demultiplexed.rds")
sample3 <- readRDS("data/rds/library3_demultiplexed.rds")
sample4 <- readRDS("data/rds/library4_demultiplexed.rds")

DefaultAssay(sample1) <- "RNA"
DefaultAssay(sample2) <- "RNA"
DefaultAssay(sample3) <- "RNA"
DefaultAssay(sample4) <- "RNA"

sample1@meta.data$multiplexed_library <- library_labels[1]
sample2@meta.data$multiplexed_library <- library_labels[2]
sample3@meta.data$multiplexed_library <- library_labels[3]
sample4@meta.data$multiplexed_library <- library_labels[4]

# ============================================================================
# Merge all samples into one object
# ============================================================================
combined <- merge(
  sample1,
  y = c(sample2, sample3, sample4),
  add.cell.ids = library_labels,
  project = "Injury_6weeks_ICs"
)
combined <- JoinLayers(combined, overwrite = TRUE)

combined$HTO_classification <- as.character(combined$HTO_classification)
combined$intervention <- unname(mouse_to_intervention[combined$HTO_classification])

if (any(is.na(combined$intervention))) {
  stop(
    "Unmapped HTO_classification values: ",
    paste(unique(combined$HTO_classification[is.na(combined$intervention)]), collapse = ", ")
  )
}

Idents(combined) <- combined$multiplexed_library

# ============================================================================
# Run ScDblFinder for doublet detection
# ============================================================================
counts_matrix <- LayerData(combined, assay = "RNA", layer = "counts")
sce <- SingleCellExperiment(list(counts = counts_matrix))
colData(sce) <- DataFrame(combined@meta.data)

sce <- scDblFinder(sce, samples = combined$multiplexed_library)

doublet_score <- sce$scDblFinder.score
names(doublet_score) <- colnames(sce)

doublet_class <- sce$scDblFinder.class
names(doublet_class) <- colnames(sce)

common_cells <- intersect(colnames(combined), names(doublet_score))

doublet_score_matched <- rep(NA_real_, ncol(combined))
names(doublet_score_matched) <- colnames(combined)
doublet_score_matched[common_cells] <- as.numeric(doublet_score[common_cells])

doublet_class_matched <- rep(NA_character_, ncol(combined))
names(doublet_class_matched) <- colnames(combined)
doublet_class_matched[common_cells] <- as.character(doublet_class[common_cells])

combined$doublet_score_scDblFinder <- doublet_score_matched
combined$doublet_class_scDblFinder <- doublet_class_matched

combined <- combined[, combined$doublet_class_scDblFinder == "singlet"]

# ============================================================================
# QC Filtering
# ============================================================================
combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^mt-")

percent_mt_cutoff <- 10
nFeature_RNA_min <- 200
nFeature_RNA_max <- 6000
nCount_RNA_max <- 40000

combined <- combined[, combined$percent.mt < percent_mt_cutoff]
combined <- combined[, combined$nFeature_RNA >= nFeature_RNA_min]
combined <- combined[, combined$nFeature_RNA <= nFeature_RNA_max]
combined <- combined[, combined$nCount_RNA < nCount_RNA_max]

# ============================================================================
# Harmony Integration
# ============================================================================
combined <- RunHarmony(combined, group.by.vars = "multiplexed_library", plot_convergence = TRUE, verbose = FALSE)

# ============================================================================
# Save final integrated object
# ============================================================================
saveRDS(combined, file = output_rds)
