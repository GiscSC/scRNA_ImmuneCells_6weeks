# scRNA-Seq. Injury 6 weeks - Immune Cells
# Script 05: Compositional_analysis
# written by Gideon JL Schaefer
# 2026-08-14
# ============================================================================

### define paths
setwd("your_wd")
input_rds <- "data/rds/annotated_RDS.rds"

### activate venv
source("renv/activate.R")

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(writexl)
library(openxlsx)
library(dittoSeq)
library(speckle)
library(RColorBrewer)

data <- readRDS(input_rds)
Idents(data) <- data$cluster_annotation

# Metadata from script 02: HTO_classification = Mouse1-Mouse11; multiplexed_library = library1-library4; intervention = Injury/Sham
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
data$cluster_annotation <- factor(as.character(data$cluster_annotation), levels = cluster_order)
data$intervention <- factor(as.character(data$intervention), levels = c("Sham", "Injury"))

cluster_colors <- setNames(
  RColorBrewer::brewer.pal(length(cluster_order), "Set2"),
  cluster_order
)

# ============================================================================
# Propeller statistics (speckle)
# ============================================================================

prop <- speckle::propeller(
  clusters = data$cluster_annotation,
  sample = data$HTO_classification,
  group = data$intervention,
  transform = "asin",
  robust = TRUE,
  trend = TRUE
)

# extract statistics of propeller (already FDR-corrected)
prop_res <- as.data.frame(prop)

# ============================================================================
# DittoSeq
# ============================================================================

p_ditto <- dittoSeq::dittoBarPlot(
  data,
  var = "cluster_annotation",
  group.by = "intervention",
  scale = "percent",
  color.panel = cluster_colors[levels(data$cluster_annotation)],
  ylab = "Proportion of cells (%)",
  main = "Compositional analysis by intervention (propeller)",
  legend.title = "Cluster"
)


