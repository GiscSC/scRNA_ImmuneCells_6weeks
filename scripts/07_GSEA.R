# scRNA-Seq. Injury 6 weeks - Immune Cells
# Script 07: Gene set enrichement analysis
# written by Gideon JL Schaefer
# 2026-08-14
# ============================================================================

### define paths
setwd("your_wd")

### activate venv
source("renv/activate.R")

library(dplyr)
library(ggplot2)
library(openxlsx)
library(gprofiler2)
library(stringr)

## load DEGs (from DESeq2-DEG analysis) and filter for only significantly reg. genes
deg <- load("DESEQ2_DEGs")
deg <- deg %>% filter(p_val_adj < 0.045)


## ============================================================================
## Run gProfiler per cluster
## ============================================================================

clusters <- unique(deg$cluster)
all_rows <- list()

# for-loop for gprofiler2
for (cl in clusters) {
  deg_cl <- deg[deg$cluster == cl, , drop = FALSE]
  
  genes_up <- deg_cl$gene[deg_cl$avg_log2FC > 0]
  genes_down <- deg_cl$gene[deg_cl$avg_log2FC < 0]
  
  gost_up <- gost(query = genes_up,
                  organism = "mmusculus",
                  ordered_query = FALSE,
                  significant = TRUE,
                  measure_underrepresentation = FALSE,
                  user_threshold = 0.05,
                  evcodes = TRUE)
  
  
  gost_down <- gost(query = genes_down,
                    organism = "mmusculus",
                    ordered_query = FALSE,
                    significant = TRUE,
                    measure_underrepresentation = FALSE,
                    user_threshold = 0.05,
                    evcodes = TRUE)
  
  gost_both <- dplyr::bind_rows(gost_up, gost_down)
  all_rows[[cl]] <- gost_both
}

gost_all <- dplyr::bind_rows(all_rows)



