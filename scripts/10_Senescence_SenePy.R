# scRNA-Seq. VPS 6 weeks - Immune Cells
# Script 10: SenePy
# written by Gideon JL Schaefer
# 2026-08-14
# ============================================================================

### define paths
setwd("your_wd")
input_rds <- "data/rds/annotated_RDS.rds"

### activate venv
source("renv/activate.R")

#libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(writexl)
library(tidyr)
library(readxl)
library(reticulate)

# ============================================================================
# load RDS
# ============================================================================
seurat_obj <- readRDS(input_rds)
DefaultAssay(seurat_obj) <- "RNA"
Idents(seurat_obj) <- seurat_obj$cluster_annotation

# Metadata from script 02: HTO_classification = Mouse1-Mouse11; intervention = VPS/Sham
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
seurat_obj$cluster_annotation <- factor(as.character(seurat_obj$cluster_annotation), levels = cluster_order)
seurat_obj$intervention <- factor(as.character(seurat_obj$intervention), levels = c("Sham", "VPS"))

# ============================================================================
# STEP 5: SenePy module score - weighted (Python senepy via reticulate)
# ============================================================================
# Python deps in reticulate env: pip install senepy anndata scipy
# SenePy R tutorial: https://github.com/jaleesr/SenePy/blob/main/EXAMPLE_OIS_hepatocytes.md

sp <- import("senepy")
anndata <- import("anndata")
scipy_sparse <- import("scipy.sparse")
pd <- import("pandas")

hubs <- sp$load_hubs(species = "Mouse", sig_type = "cell_type")
senepy_merge_results <- hubs$merge_hubs(
  hubs$metadata,
  calculate_thresh = TRUE,
  p_thres = 0.01,
  new_name = "Universal"
)

universal_hub <- hubs$hubs[["Universal"]]
senepy_universal_genes <- vapply(universal_hub, function(x) x[[1]], character(1))
senepy_universal_weights <- vapply(universal_hub, function(x) as.numeric(x[[2]]), numeric(1))

senepy_weight_tbl <- data.frame(
  gene_symbol = senepy_universal_genes,
  weight = senepy_universal_weights,
  stringsAsFactors = FALSE
)

## check for missing genes
senepy_genes_present_weighted <- intersect(senepy_universal_genes, rownames(seurat_obj))
senepy_genes_missing_weighted <- setdiff(senepy_universal_genes, senepy_genes_present_weighted)
senepy_merge_results_tbl <- py_to_r(senepy_merge_results)


## extract scRNA-IC count matrix 
counts_mat <- LayerData(seurat_obj, assay = "RNA", layer = "counts")

# transfrom to python format
adata <- anndata$AnnData(
  X = scipy_sparse$csr_matrix(t(counts_mat)),
  obs = seurat_obj@meta.data,
  var = pd$DataFrame(index = rownames(seurat_obj))
)

# calculate module expression 
senepy_weighted_scores <- unlist(sp$score_hub(adata, universal_hub))
# add module expression as coldata to seurat_obj
seurat_obj$SenePy_weighted <- senepy_weighted_scores

###Density plot with SenePy-Score by cell cluster 
intervention_colors <- c("Sham" = "#1F4E79", "VPS" = "#E64A19")

senepy_density_df <- seurat_obj@meta.data %>%
  dplyr::mutate(
    cluster_annotation = factor(as.character(cluster_annotation), levels = cluster_order),
    intervention = factor(as.character(intervention), levels = c("Sham", "VPS"))
  )

p_senepy_density <- ggplot(senepy_density_df, aes(x = SenePy_weighted, color = intervention, fill = intervention)) +
  geom_density(alpha = 0.25, linewidth = 0.6) +
  facet_wrap(~ cluster_annotation, scales = "free_y", ncol = 4) +
  scale_color_manual(values = intervention_colors) +
  scale_fill_manual(values = intervention_colors) +
  labs(
    x = "SenePy weighted score",
    y = "Density",
    title = "SenePy weighted score density by cell type"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

