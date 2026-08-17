# scRNA_ImmuneCells_6weeks

scRNA-seq analysis of peripherial CD45+ immune cells in Injury and Sham mice after 6 weeks

## Environment

- R 4.4.0
- Bioconductor 3.19
- Package versions: `renv.lock`

Set the working directory to this folder and restore the library:

```r
# scripts/00_renv_restore.R
setwd("your_wd")
source("renv/activate.R")
renv::restore()
```

Script 10 (SenePy) additionally needs Python packages `senepy`, `anndata` and `scipy` (via `reticulate`).

## Scripts

Run `scripts/01` to `scripts/11` in order. Adapt `setwd` / `input_rds` / `output_rds` at the top of each script.

Raw counts, HTO libraries and processed RDS files can be downloaded from the corresponding Zenodo repository (see `DATA_ACCESS.md`).
