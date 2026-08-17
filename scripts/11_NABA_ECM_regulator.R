# scRNA-Seq. Injury 6 weeks - Immune Cells
# Script 11: NABA_ECM
# written by Gideon JL Schaefer
# 2026-08-14
# ============================================================================

### define paths
setwd("your_wd")
input_rds <- "data/rds/annotated_RDS.rds"

### activate venv
source("renv/activate.R")

library(Seurat)
library(dplyr)
library(ggplot2)
library(msigdbr)
library(writexl)


# ============================================================================
# Load annotated object
# ============================================================================

seurat_obj <- readRDS(input_rds)
DefaultAssay(seurat_obj) <- "RNA"
Idents(seurat_obj) <- seurat_obj$cluster_annotation

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
seurat_obj$cluster_annotation <- factor(as.character(seurat_obj$cluster_annotation), levels = cluster_order)
seurat_obj$intervention <- factor(as.character(seurat_obj$intervention), levels = c("Sham", "Injury"))

# ============================================================================
# Download NABA_ECM_REGULATORS gene set (MSigDB mouse via msigdbr)
# ============================================================================

naba_genes_df <- msigdbr::msigdbr(species = "Mus musculus") %>%
  dplyr::filter(gs_name == "NABA_ECM_REGULATORS")

naba_genes_all <- unique(naba_genes_df$gene_symbol)
naba_genes_present <- intersect(naba_genes_all, rownames(seurat_obj))

# ============================================================================
# Module score
# ============================================================================

seurat_obj <- AddModuleScore(
  object = seurat_obj,
  features = list(naba_genes_present),
  name = "NABA_ECM_Regulator"
)

# ============================================================================
# Density plot by cell cluster (Sham vs Injury)
# ============================================================================

intervention_colors <- c("Sham" = "#1F4E79", "Injury" = "#E64A19")

naba_density_df <- seurat_obj@meta.data %>%
  dplyr::mutate(
    cluster_annotation = factor(as.character(cluster_annotation), levels = cluster_order),
    intervention = factor(as.character(intervention), levels = c("Sham", "Injury"))
  )

p_naba_density <- ggplot(naba_density_df, aes(x = NABA_ECM_Regulator1, color = intervention, fill = intervention)) +
  geom_density(alpha = 0.25, linewidth = 0.6) +
  facet_wrap(~ cluster_annotation, scales = "free_y", ncol = 4) +
  scale_color_manual(values = intervention_colors) +
  scale_fill_manual(values = intervention_colors) +
  labs(
    x = "NABA ECM Regulator ModuleScore",
    y = "Density",
    title = "NABA ECM Regulator ModuleScore density by cell type"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )


