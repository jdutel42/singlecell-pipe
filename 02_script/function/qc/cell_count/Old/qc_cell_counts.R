# ============================================================
# QC Cell counts - scRNA-seq
# Seurat
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

# ------------------------------------------------------------
# 1. Table of cell counts
# ------------------------------------------------------------
qc_cell_count_table <- function(
    seu,
    group_col
) {
  stopifnot(group_col %in% colnames(seu@meta.data))
  seu@meta.data %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::summarise(cell_count = dplyr::n(), .groups = "drop")
}

# ------------------------------------------------------------
# 2. Barplots of cell counts
# ------------------------------------------------------------
qc_cell_count_barplot <- function(
    seu,
    group_cols,
    outdir = NULL,
    dataset_name = "dataset",
    palette = NULL
) {
  plots <- list()
  i <- 1
  
  if (!is.null(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  }
  
  for (group_col in group_cols) {
    stopifnot(group_col %in% colnames(seu@meta.data))
    
    p <- ggplot(
      seu@meta.data,
      aes(x = .data[[group_col]], fill = .data[[group_col]])
    ) +
      geom_bar() +
      theme_classic() +
      theme(
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")
      ) +
      ggtitle(paste("NCells -", group_col)) +
      NoLegend()
    
    if (!is.null(palette)) {
      p <- p + scale_fill_manual(values = palette, drop = FALSE)
    }
    
    plots[[i]] <- p
    i <- i + 1
    
    if (!is.null(outdir)) {
      ggsave(
        plot = p,
        filename = file.path(
          outdir,
          paste0("Cell_Counts_By_", group_col, "_", dataset_name, "_Before_QC_bis.png")
        ),
        width = 6,
        height = 4
      )
    }
  }
  
  return(plots)
}

# ------------------------------------------------------------
# 3. UMAP facetted by group
# ------------------------------------------------------------
qc_cell_count_umap <- function(
    seu,
    group_col,
    outdir = NULL,
    dataset_name = "dataset",
    reduction = "umap",
    ncol = 6
) {
  stopifnot(group_col %in% colnames(seu@meta.data))
  
  p <- DimPlot(seu, group.by = group_col, reduction = reduction) +
    facet_wrap(stats::as.formula(paste("~", group_col)), ncol = ncol) +
    ggtitle(paste("UMAP by", group_col, "before QC")) +
    NoLegend()
  
  if (!is.null(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    ggsave(
      plot = p,
      filename = file.path(
        outdir,
        paste0("UMAP_Cell_Counts_by_", group_col, "_", dataset_name, "_Before_QC_bis.png")
      ),
      width = 23,
      height = 12
    )
  }
  
  return(p)
}

# ------------------------------------------------------------
# 4️⃣ Fonction principale généralisée pour n'importe quel subset
# ------------------------------------------------------------
qc_cell_counts <- function(
    seu_subset,
    dataset_name = "dataset",
    features_barplot = c("condition", "id", "emulsion"),
    umap_group_col = "emulsion",
    outdir = NULL,
    palette = NULL,
    reduction = "umap",
    umap_ncol = 6
) {
  all_plots <- list()
  i <- 1
  
  # --- Table résumé par emulsion ---
  print(qc_cell_count_table(seu_subset, umap_group_col))
  
  # --- Barplots ---
  barplots <- qc_cell_count_barplot(
    seu = seu_subset,
    group_cols = features_barplot,
    outdir = if (!is.null(outdir)) file.path(outdir, dataset_name, "QC", "Cell_Counts", "Before_QC") else NULL,
    dataset_name = dataset_name,
    palette = palette
  )
  
  for (p in barplots) {
    all_plots[[i]] <- p
    i <- i + 1
  }
  
  # --- UMAP facetted ---
  p_umap <- qc_cell_count_umap(
    seu = seu_subset,
    group_col = umap_group_col,
    outdir = if (!is.null(outdir)) file.path(outdir, dataset_name, "QC", "Cell_Counts", "Before_QC") else NULL,
    dataset_name = dataset_name,
    reduction = reduction,
    ncol = umap_ncol
  )
  
  all_plots[[i]] <- p_umap
  
  return(all_plots)
}