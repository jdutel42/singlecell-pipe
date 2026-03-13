##############################################
## scRNA-seq QC — Complete Global Wrapper  ##
##############################################

# Author: Jordan Dutel
# Purpose: Orchestrate the execution of all scRNA-seq QC modules
# Usage: 
#   source("qc_full_qc.R")
#   seu <- run_full_qc(seu, DIR_QC, celltype = "CD4_Tcells", nmads = 3)
#
# Inputs:
# - seu: a preprocessed Seurat object (normalized, PCA and UMAP computed)
# - DIR_QC: root output directory for QC plots
# - celltype: dataset name (e.g. "CD4_Tcells") used for output organization
# - nmads: number of median absolute deviations used to define outliers (default 3)
#
# Outputs:
# - Seurat object updated with QC metadata columns
# - QC plots saved under DIR_QC/<celltype>/ subdirectories

#==============================
# Load QC modules
#==============================

# Adjust these paths according to your project structure
source("qc/qc_mito.R")
source("qc/qc_ribo.R")
source("qc/qc_libsize.R")
source("qc/qc_libdiv.R")
source("qc/qc_libcomplex.R")
source("qc/qc_cell_cycle.R")

#==============================
# Full QC wrapper function
#==============================

run_full_qc <- function(seu, DIR_QC, celltype = "Cells", nmads = 3) {
  stopifnot(inherits(seu, "Seurat"))
  
  message("Starting full QC for dataset: ", celltype)
  
  # Mitochondrial QC
  message(" - Running mitochondrial QC")
  seu <- run_mito_qc(seu, DIR_QC, celltype, nmads)
  
  # Ribosomal QC
  message(" - Running ribosomal QC")
  seu <- run_ribo_qc(seu, DIR_QC, celltype, nmads)
  
  # Library size QC (nCount_RNA)
  message(" - Running library size QC")
  seu <- run_libsize_qc(seu, DIR_QC, celltype, nmads)
  
  # Library diversity QC (nFeature_RNA)
  message(" - Running library diversity QC")
  seu <- run_libdiv_qc(seu, DIR_QC, celltype, nmads)
  
  # Library complexity QC (log10 Genes per UMI)
  message(" - Running library complexity QC")
  seu <- run_libcomplex_qc(seu, DIR_QC, celltype, nmads)
  
  # Cell cycle scoring
  message(" - Running cell cycle scoring")
  seu <- run_cellcycle_qc(seu, DIR_QC, celltype)
  
  message("Full QC completed for dataset: ", celltype)
  
  return(seu)
}
