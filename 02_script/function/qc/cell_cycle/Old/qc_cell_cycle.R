############################################
## scRNA-seq QC — Cell Cycle (Standalone)   ##
############################################

# Author: Jordan
# Purpose: Cell cycle scoring QC module
# Usage: source("qc_cell_cycle.R")
# Input: Seurat object with UMAP computed
# Output: Updated Seurat object + QC plots

#==============================
# Libraries
#==============================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

#==============================
# Run: Cell cycle scoring
#==============================

compute_cellcycle_qc <- function(
    seu,
    s_genes = cc.genes$s.genes,
    g2m_genes = cc.genes$g2m.genes
) {
  stopifnot(inherits(seu, "Seurat"))
  
  seu <- CellCycleScoring(
    seu,
    s.features = s_genes,
    g2m.features = g2m_genes,
    set.ident = FALSE
  )
  
  return(seu)
}

#==============================
# Visualize: Cell cycle on UMAP
#==============================

plot_cellcycle_umap <- function(
    seu,
    DIR_QC,
    celltype = "Cells",
    reduction = "umap"
) {
  stopifnot(inherits(seu, "Seurat"))
  
  outdir <- file.path(DIR_QC, celltype, "Cell_Cycle")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  p <- DimPlot(seu, reduction = reduction, group.by = "Phase") +
    ggtitle(paste0("Cell cycle phase before QC ", celltype)) +
    facet_wrap(~Phase) +
    NoLegend()
  
  print(p)
  
  ggsave(
    plot = p,
    filename = file.path(
      outdir,
      paste0("UMAP_Cell_Cycle_Phase_before_QC_", celltype, ".png")
    ),
    width = 12,
    height = 6
  )
}

#==============================
# Wrapper (optional)
#==============================

run_cellcycle_qc <- function(
    seu,
    DIR_QC,
    celltype = "Cells"
) {
  seu <- compute_cellcycle_qc(seu)
  plot_cellcycle_umap(seu, DIR_QC = DIR_QC, celltype = celltype)
  return(seu)
}

#==============================
# Example
#==============================
# source("qc_cell_cycle.R")
# seu_CD4 <- run_cellcycle_qc(
#   seu_CD4,
#   DIR_QC = DIR_QC,
#   celltype = "CD4_Tcells"
# )
