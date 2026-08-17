# scRNA-Seq. Injury 6 weeks - Immune Cells
# Script 09: MultiNicheNet
# written by Gideon JL Schaefer
# 2026-08-14
# ============================================================================

### define paths
setwd("your_wd")
input_rds <- "data/rds/annotated_RDS.rds"
ligand_target_matrix_path <- "path_to_ligand_target_matrix.rds"
output_dir <- "your_output_dir"

### activate venv
source("renv/activate.R")

### load libraries
library(Seurat)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(RColorBrewer)
library(circlize)
library(writexl)
library(nichenetr)
library(multinichenetr)

# ============================================================================
# Load Seurat and generate SCE
# ============================================================================
seurat_obj <- readRDS(input_rds)
DefaultAssay(seurat_obj) <- "RNA"

# Metadata from script 02:
# HTO_classification = Mouse1-Mouse11 (sample_id for MultiNicheNet pseudobulk)
# multiplexed_library = library1-library4
# intervention = Injury / Sham

sce <- Seurat::as.SingleCellExperiment(seurat_obj, assay = "RNA")
sce <- alias_to_symbol_SCE(sce, organism = "mouse") %>% makenames_SCE()

cd <- SummarizedExperiment::colData(sce)
cd$intervention <- factor(make.names(as.character(cd$intervention)))
cd$HTO_classification <- make.names(as.character(cd$HTO_classification)) ### Mouse1-Mouse11
cd$cluster_annotation <- make.names(as.character(cd$cluster_annotation))
if ("multiplexed_library" %in% colnames(cd)) {
  cd$multiplexed_library <- make.names(as.character(cd$multiplexed_library)) ### library1-library4
}
colData(sce) <- cd
# ============================================================================
# Load NicheNet's ligand-receptor network and ligand-target matrix
# ============================================================================
lr_network_all = readRDS(url(
  "https://zenodo.org/record/10229222/files/lr_network_mouse_allInfo_30112033.rds"
)) %>% 
  mutate(
    ligand = convert_alias_to_symbols(ligand, organism = organism), 
    receptor = convert_alias_to_symbols(receptor, organism = organism))

lr_network_all = lr_network_all  %>% 
  mutate(ligand = make.names(ligand), receptor = make.names(receptor)) 
lr_network = lr_network_all %>% 
  distinct(ligand, receptor)

ligand_target_matrix <- readRDS(ligand_target_matrix_path)

colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% 
  convert_alias_to_symbols(organism = organism) %>% make.names()
rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% 
  convert_alias_to_symbols(organism = organism) %>% make.names()

lr_network = lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
ligand_target_matrix = ligand_target_matrix[, lr_network$ligand %>% unique()]

# ============================================================================
# STEP 3: Define meta parameters 
# ============================================================================

# Define Metadata columns in Seurat object
sample_id   <- make.names("HTO_classification")
group_id    <- make.names("intervention")
celltype_id <- make.names("cluster_annotation")

##Define covariates and batch definition -> No definition in this dataset
covariates = NA
batches <- NA

##Deifne contrasts of interest

contrasts_oi = c("'Injury-Sham','Sham-Injury'")

contrast_tbl <- tibble::tibble(
  contrast = c("Injury-Sham", "Sham-Injury"),
  group      = c("Injury", "Sham")
)

## Define the sender and receiver cell types of interest.
senders_oi   <- unique(sce$cluster_annotation)
receivers_oi <- unique(sce$cluster_annotation)
#sce = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% c(senders_oi, receivers_oi) ] ###skipped as no cells are excluded

# ============================================================================
# STEP 4: Run MultiNicheNet core analysis
# ============================================================================

## Cell type filtering 
min_cells = 10

abundance_info = get_abundance_info(
  sce = sce, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  min_cells = min_cells, 
  senders_oi = senders_oi, receivers_oi = receivers_oi, 
  batches = batches
)


abundance_info$abund_plot_sample
ggsave(filename = file.path(output_dir,"abundance_filtering.jpeg"), dpi = 300)
ggsave(filename = file.path(output_dir,"abundance_filtering.svg"))

###check for sparsity in cells independent and specific for conditions
sample_group_celltype_df = abundance_info$abundance_data %>% 
  filter(n > min_cells) %>% 
  ungroup() %>% 
  distinct(sample_id, group_id) %>% 
  cross_join(
    abundance_info$abundance_data %>% 
      ungroup() %>% 
      distinct(celltype_id)
  ) %>% 
  arrange(sample_id)

abundance_df = sample_group_celltype_df %>% left_join(
  abundance_info$abundance_data %>% ungroup()
)

abundance_df$n[is.na(abundance_df$n)] = 0
abundance_df$keep[is.na(abundance_df$keep)] = FALSE
abundance_df_summarized = abundance_df %>% 
  mutate(keep = as.logical(keep)) %>% 
  group_by(group_id, celltype_id) %>% 
  summarise(samples_present = sum((keep)))

celltypes_absent_one_condition = abundance_df_summarized %>% 
  filter(samples_present == 0) %>% pull(celltype_id) %>% unique() 
# find truly condition-specific cell types by searching for cell types 
# truely absent in at least one condition

celltypes_present_one_condition = abundance_df_summarized %>% 
  filter(samples_present >= 2) %>% pull(celltype_id) %>% unique() 
# require presence in at least 2 samples of one group so 
# it is really present in at least one condition

condition_specific_celltypes = intersect(
  celltypes_absent_one_condition, 
  celltypes_present_one_condition)

total_nr_conditions = SummarizedExperiment::colData(sce)[,group_id] %>% 
  unique() %>% length() 

absent_celltypes = abundance_df_summarized %>% 
  filter(samples_present < 2) %>% 
  group_by(celltype_id) %>% 
  count() %>% 
  filter(n == total_nr_conditions) %>% 
  pull(celltype_id)

print("condition-specific celltypes:")
## [1] "condition-specific celltypes:"
print(condition_specific_celltypes)
## character(0)

print("absent celltypes:")
## [1] "absent celltypes:"
print(absent_celltypes)
## character(0)

###no cells are excluded from the analysis 
analyse_condition_specific_celltypes = FALSE
if(analyse_condition_specific_celltypes == TRUE){
  senders_oi = senders_oi %>% setdiff(absent_celltypes)
  receivers_oi = receivers_oi %>% setdiff(absent_celltypes)
} else {
  senders_oi = senders_oi %>% 
    setdiff(union(absent_celltypes, condition_specific_celltypes))
  receivers_oi = receivers_oi %>% 
    setdiff(union(absent_celltypes, condition_specific_celltypes))
}

sce = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% 
            c(senders_oi, receivers_oi)
]

## Gene filtering 
min_sample_prop = 0.50 #(only consider genes expressed in at least 50% of the samples)
fraction_cutoff = 0.05 #(only consider genes expressed in at least 5% of the cells in a sample)

##calculate gene matrix based on defined cutoffs
frq_list = get_frac_exprs(
  sce = sce, 
  sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id, 
  batches = batches, 
  min_cells = min_cells, 
  fraction_cutoff = fraction_cutoff, min_sample_prop = min_sample_prop)

genes_oi = frq_list$expressed_df %>% 
  filter(expressed == TRUE) %>% pull(gene) %>% unique() 
sce = sce[genes_oi, ]

## Pseudobulking and DEG analysis

###Pseudobulking
abundance_expression_info = process_abundance_expression_info(
  sce = sce, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  min_cells = min_cells, 
  senders_oi = senders_oi, receivers_oi = receivers_oi, 
  lr_network = lr_network, 
  batches = batches, 
  frq_list = frq_list, 
  abundance_info = abundance_info)

DE_info = get_DE_info(
  sce = sce, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  batches = batches, covariates = covariates, 
  contrasts_oi = contrasts_oi, 
  min_cells = min_cells, 
  expressed_df = frq_list$expressed_df)

##inspect DEGs
DE_info$celltype_de$de_output_tidy %>% head()
DE_info$hist_pvals
###distribution looks okay -> no empirical pval estimation needed
empirical_pval = FALSE

#Combine DE information for ligand-senders and receptors-receivers
celltype_de <- DE_info$celltype_de$de_output_tidy
sender_receiver_de = combine_sender_receiver_de(
  sender_de = celltype_de,
  receiver_de = celltype_de,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi,
  lr_network = lr_network
)

##Assess geneset_oi-vs-background ratios for different DE output tresholds prior to the NicheNet ligand activity analysis
logFC_threshold = 0.50
p_val_threshold = 0.05
p_val_adj = FALSE 

geneset_assessment = contrast_tbl$contrast %>% 
  lapply(
    process_geneset_data, 
    celltype_de, logFC_threshold, p_val_adj, p_val_threshold
  ) %>% 
  bind_rows() 
geneset_assessment

table(geneset_assessment$in_range_up) # TRUE16
table(geneset_assessment$in_range_down) # TRUE16

##Perform the ligand activity analysis and ligand-target inference
ligand_activities_targets_DEgenes = suppressMessages(suppressWarnings(
  get_ligand_activities_targets_DEgenes(
    receiver_de = celltype_de,
    receivers_oi = intersect(receivers_oi, celltype_de$cluster_id %>% unique()),
    ligand_target_matrix = ligand_target_matrix,
    logFC_threshold = logFC_threshold,
    p_val_threshold = p_val_threshold,
    p_val_adj = p_val_adj,
    top_n_target = 250,
    verbose = TRUE, 
    n.cores = 10
  )
))


##Prioritization: rank cell-cell communication patterns through multi-criteria prioritization
ligand_activity_down = FALSE ###(decide eacht ime individually!)

sender_receiver_tbl = sender_receiver_de %>% distinct(sender, receiver)

metadata_combined = SummarizedExperiment::colData(sce) %>% tibble::as_tibble()

grouping_tbl = metadata_combined[,c(sample_id, group_id, batches)] %>% 
    tibble::as_tibble() %>% distinct()
  colnames(grouping_tbl) = c("sample","group",batches)

prioritization_tables = suppressMessages(generate_prioritization_tables(
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  contrast_tbl = contrast_tbl,
  sender_receiver_tbl = sender_receiver_tbl,
  grouping_tbl = grouping_tbl,
  scenario = "regular", # all prioritization criteria will be weighted equally
  fraction_cutoff = fraction_cutoff, 
  abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
  abundance_data_sender = abundance_expression_info$abundance_data_sender,
  ligand_activity_down = ligand_activity_down
))

##inspect
prioritization_tables$group_prioritization_tbl %>% head(20)

##Calculate the across-samples expression correlation between ligand-receptor pairs and target genes
lr_target_prior_cor = lr_target_prior_cor_inference(
  receivers_oi = prioritization_tables$group_prioritization_tbl$receiver %>% unique(), 
  abundance_expression_info = abundance_expression_info, 
  celltype_de = celltype_de, 
  grouping_tbl = grouping_tbl, 
  prioritization_tables = prioritization_tables, 
  ligand_target_matrix = ligand_target_matrix, 
  logFC_threshold = logFC_threshold, 
  p_val_threshold = p_val_threshold, 
  p_val_adj = p_val_adj
)

##Save all the output of MultiNicheNet

multinichenet_output = list(
  celltype_info = abundance_expression_info$celltype_info,
  celltype_de = celltype_de,
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de =  sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  prioritization_tables = prioritization_tables,
  grouping_tbl = grouping_tbl,
  lr_target_prior_cor = lr_target_prior_cor
) 

##Visualization of differential cell-cell interactions
  
  ### ChordDiagramms
  prioritized_tbl_oi_all = get_top_n_lr_pairs(
    multinichenet_output$prioritization_tables, 
    top_n = 50, 
    rank_per_group = FALSE
  )
  
  prioritized_tbl_oi = 
    multinichenet_output$prioritization_tables$group_prioritization_tbl %>%
    filter(id %in% prioritized_tbl_oi_all$id) %>%
    distinct(id, sender, receiver, ligand, receptor, group) %>% 
    left_join(prioritized_tbl_oi_all)
  prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] = 0
  
  senders_receivers = union(prioritized_tbl_oi$sender %>% unique(), prioritized_tbl_oi$receiver %>% unique()) %>% sort()
  
  colors_sender = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)
  colors_receiver = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)
  
  circos_list = make_circos_group_comparison(prioritized_tbl_oi, colors_sender, colors_receiver)
  

  ###Specific bubble plots for BCs as receivers // GOI = "Sham"
  group_oi = "Sham"
  
   prioritized_tbl_oi_Injury_50 = get_top_n_lr_pairs(
    multinichenet_output$prioritization_tables, 
    50, 
    groups_oi = group_oi, 
    receivers_oi = "BCs")
  
  plot_oi = make_sample_lr_prod_activity_plots_Omnipath(
    multinichenet_output$prioritization_tables, 
    prioritized_tbl_oi_Injury_50 %>% inner_join(lr_network_all))
  plot_oi
  
  ggsave(paste0(output_dir,"/sample_lr_prod_activity_plot_GOI_Sham_Reciever_BCs.svg"), plot_oi,
         width = 24, height = 8, limitsize = FALSE)
  
  ###Specific bubble plots for BCs as receivers // GOI = "Injury"
  group_oi = "Injury"
  
  prioritized_tbl_oi_Injury_50 = get_top_n_lr_pairs(
    multinichenet_output$prioritization_tables, 
    50, 
    groups_oi = group_oi, 
    receivers_oi = "BCs")
  
  plot_oi = make_sample_lr_prod_activity_plots_Omnipath(
    multinichenet_output$prioritization_tables, 
    prioritized_tbl_oi_Injury_50 %>% inner_join(lr_network_all))
  plot_oi
  
  ggsave(paste0(output_dir,"/sample_lr_prod_activity_plot_GOI_Injury_Reciever_BCs.svg"), plot_oi,
         width = 24, height = 8, limitsize = FALSE)
  
  
  
#####B-cell interaction
# Top 30 ligand-receptor interactions with BCs as receiver (across all senders).

prioritized_tbl_oi_BCs_50 <- get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables,
  top_n = 30,
  receivers_oi = "BCs",
  rank_per_group = T
)

senders_receivers_bcs <- union(
  prioritized_tbl_oi_BCs_50$sender %>% unique(),
  prioritized_tbl_oi_BCs_50$receiver %>% unique()
) %>% sort()

colors_sender_bcs <- RColorBrewer::brewer.pal(n = length(senders_receivers_bcs), name = "Spectral") %>%
  magrittr::set_names(senders_receivers_bcs)
colors_receiver_bcs <- colors_sender_bcs

circos_list_bcs <- make_circos_group_comparison(prioritized_tbl_oi_BCs_50, colors_sender_bcs, colors_receiver_bcs)



