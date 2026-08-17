# scRNA-Seq. Injury 6 weeks - Immune Cells
# Script 01: Demultiplexing with Hashtag Oligos (HTO)
# written by Gideon JL Schaefer
# 2026-08-14
# ============================================================================

### define paths
setwd("your_wd")

multiplexed_library <- "library1"  # library1-4
cDNA_path     <- "path_to_filtered_10X_matrix"
hashtag_path  <- "path_to_HTO_library"
output_rds    <- "data/rds/library1_demultiplexed.rds"

### activate venv
source("renv/activate.R")

## load libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(Matrix)
library(ggpubr)
library(writexl)
library(dittoSeq)
library(openxlsx)
library(stringi)
library(R.utils)

# ============================================================================
# Annotation of #-barcode to Mouse-ID (per multiplexed library)
# ============================================================================
# library1: Mouse1-3  
# library2: Mouse4-5 (demultiplexing for #AB "189-GGTCGAGAGCATTCA" failed, mouse is therefore excluded!)
# library3: Mouse6-8
# library4: Mouse9-11

hto_barcode_to_mouse <- switch(
  multiplexed_library,
  library1 = c(
    "189-GGTCGAGAGCATTCA" = "Mouse1",
    "190-CTTGCCGCATGTCAT" = "Mouse2",
    "216-AAAGCATTCTTCACG" = "Mouse3"
  ),
  library2 = c(
    "190-CTTGCCGCATGTCAT" = "Mouse4",
    "216-AAAGCATTCTTCACG" = "Mouse5"
  ),
  library3 = c(
    "189-GGTCGAGAGCATTCA" = "Mouse6",
    "190-CTTGCCGCATGTCAT" = "Mouse7",
    "216-AAAGCATTCTTCACG" = "Mouse8"
  ),
  library4 = c(
    "189-GGTCGAGAGCATTCA" = "Mouse9",
    "190-CTTGCCGCATGTCAT" = "Mouse10",
    "216-AAAGCATTCTTCACG" = "Mouse11"
  ),
  stop("Unknown multiplexed_library: ", multiplexed_library)
)

failed_hashtag_barcode <- "189-GGTCGAGAGCATTCA"  # excluded from demux in library2 only

# ============================================================================
# LOAD cDNA-Matrix
# ============================================================================
cDNA_data <- Read10X(data.dir = cDNA_path)
cDNA_data <- CreateSeuratObject(
  counts = cDNA_data,
  project = paste0(multiplexed_library, "_demultiplexed"),
  min.cells = 3,
  min.features = 200
)

# ============================================================================
# Create HTO-Matrix
# ============================================================================
barcode.path  <- file.path(hashtag_path, "barcodes.tsv.gz")
features.path <- file.path(hashtag_path, "features.tsv.gz")
matrix.path   <- file.path(hashtag_path, "matrix.mtx.gz")

hto_data <- readMM(file = matrix.path)

feature.names <- read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
barcode.names <- read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)

colnames(hto_data) <- barcode.names$V1
rownames(hto_data) <- feature.names$V1
rownames(hto_data) <- gsub("_", "-", rownames(hto_data))

for (barcode in names(hto_barcode_to_mouse)) {
  rownames(hto_data) <- gsub(barcode, hto_barcode_to_mouse[[barcode]], rownames(hto_data), fixed = TRUE)
}

hto_data <- hto_data[rownames(hto_data) != "unmapped", ]

# ============================================================================
# merge HTO- and cDNA-matrices
# ============================================================================
cDNA_data <- RenameCells(object = cDNA_data, new.names = gsub("-1", "", as.list(colnames(cDNA_data))))
shared_barcodes <- intersect(colnames(cDNA_data), colnames(hto_data))

cDNA_data <- cDNA_data[, shared_barcodes]
sample.htos <- as.matrix(hto_data[, shared_barcodes])

cDNA_HTOs <- CreateSeuratObject(
  counts = cDNA_data@assays[["RNA"]],
  project = multiplexed_library,
  min.cells = 3,
  min.features = 200
)

cDNA_HTOs[["HTO"]] <- CreateAssayObject(counts = sample.htos[, colnames(cDNA_HTOs)])
cDNA_HTOs <- NormalizeData(cDNA_HTOs, assay = "HTO", normalization.method = "CLR")

# library2: remove failed hashtag before demultiplexing
if (multiplexed_library == "library2") {
  hto_counts <- LayerData(cDNA_HTOs, assay = "HTO", layer = "counts")
  failed_feature <- gsub("_", "-", failed_hashtag_barcode)
  if (failed_feature %in% rownames(hto_counts)) {
    hto_counts <- hto_counts[!rownames(hto_counts) %in% failed_feature, , drop = FALSE]
    cDNA_HTOs[["HTO"]] <- CreateAssayObject(counts = hto_counts)
    cDNA_HTOs <- NormalizeData(cDNA_HTOs, assay = "HTO", normalization.method = "CLR")
  }
}

# ============================================================================
# HTODemux + singlet selection
# ============================================================================
cDNA_HTOs <- HTODemux(
  cDNA_HTOs,
  assay = "HTO",
  kfunc = "clara",
  positive.quantile = 0.95,
  verbose = TRUE
)

Idents(cDNA_HTOs) <- "HTO_classification.global"
cDNA_HTOs <- subset(cDNA_HTOs, idents = "Singlet", invert = FALSE)

cDNA_HTOs$multiplexed_library <- multiplexed_library
cDNA_HTOs$HTO_classification <- as.character(cDNA_HTOs$HTO_classification)

saveRDS(object = cDNA_HTOs, file = output_rds)
