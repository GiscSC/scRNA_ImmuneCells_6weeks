# scRNA-Seq. VPS 6 weeks - Immune Cells
# Script 00: Restore renv environment
# written by Gideon JL Schaefer
# 2026-08-14
# ============================================================================

### define paths
setwd("your_wd")

### activate venv
source("renv/activate.R")

renv::restore()
