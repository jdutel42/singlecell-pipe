# ============================================================
# QC Ribosomal content - scRNA-seq
# Seurat + scater
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(scater)
  library(dplyr)
  library(ggplot2)
  library(ggridges)
  library(ggbeeswarm)
  library(viridis)
})

# ------------------------------------------------------------
# 1. Compute ribosomal QC + MAD outliers (high + low)
# ------------------------------------------------------------
qc_ribo_compute <- function(
    seu,
    batch_col,
    assay = "RNA",
    ribo_pattern = "^RPL|^RPS",
    nmads = 3,
    meta_colname = "ribo_status"
) {
  
  stopifnot(batch_col %in% colnames(seu@meta.data))
  
  counts <- seu@assays[[assay]]@counts
  
  qc_ribo <- perCellQCMetrics(
    counts,
    subsets = list(
      ribo = grep(ribo_pattern, rownames(counts), value = TRUE)
    )
  )
  
  seu <- AddMetaData(seu, as.data.frame(qc_ribo))
  
  ribo_percent <- seu@meta.data$subsets_ribo_percent
  
  qc_df <- seu@meta.data %>%
    dplyr::select(
      !!batch_col,
      subsets_ribo_percent
    ) %>%
    dplyr::rename(batch = !!batch_col) %>%
    dplyr::mutate(
      high_ribo_outlier = isOutlier(
        ribo_percent,
        nmads = nmads,
        type = "higher"
      ),
      low_ribo_outlier = isOutlier(
        ribo_percent,
        nmads = nmads,
        type = "lower"
      ),
      !!meta_colname := dplyr::case_when(
        high_ribo_outlier ~ "high",
        low_ribo_outlier ~ "low",
        TRUE ~ "normal"
      )
    )
  
  seu <- AddMetaData(
    seu,
    metadata = qc_df[[meta_colname]],
    col.name = meta_colname
  )
  
  return(seu)
}

# ------------------------------------------------------------
# 2. Visualization
# ------------------------------------------------------------
qc_ribo_plots <- function(
    seu,
    batch_col,
    outdir,
    dataset_name = "dataset",
    reduction = "umap",
    meta_ribo = "subsets_ribo_percent",
    meta_status = "ribo_status"
) {
  
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  meta <- seu@meta.data
  
  # ---- Ridge plot
  p1 <- ggplot(
    meta,
    aes(
      x = .data[[meta_ribo]],
      y = .data[[batch_col]],
      fill = after_stat(x)
    )
  ) +
    geom_density_ridges_gradient(
      scale = 3,
      rel_min_height = 0.01
    ) +
    scale_fill_viridis(
      alpha = 0.9,
      option = "C",
      guide = "none"
    ) +
    labs(
      title = "Ribosomal content by batch",
      x = "Percent ribosomal reads",
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.title = element_text(face = "bold")
    )
  
  ggsave(
    plot = p1,
    filename = file.path(
      outdir,
      paste0(
        "Ridgeplot_ribosomal_content_by_",
        batch_col,
        "_",
        dataset_name,
        ".png"
      )
    ),
    width = 8,
    height = 5
  )
  
  # ---- UMAP
  p2 <- FeaturePlot(
    seu,
    features = meta_ribo,
    reduction = reduction,
    order = TRUE,
    cols = c("lightgrey", "red")
  )
  
  ggsave(
    plot = p2,
    filename = file.path(
      outdir,
      paste0(
        "UMAP_ribosomal_content_",
        dataset_name,
        ".png"
      )
    ),
    width = 6,
    height = 5
  )
  
  # ---- Violin
  p3 <- ggplot(
    meta,
    aes(
      x = .data[[batch_col]],
      y = .data[[meta_ribo]],
      color = .data[[meta_status]]
    )
  ) +
    geom_violin(fill = "grey90", color = NA) +
    geom_quasirandom(alpha = 0.5, size = 0.8) +
    scale_y_log10() +
    scale_color_manual(
      values = c(
        "normal" = "green3",
        "high"   = "red3",
        "low"    = "blue3"
      )
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(
        angle = 45,
        vjust = 1,
        hjust = 1
      )
    )
  
  ggsave(
    plot = p3,
    filename = file.path(
      outdir,
      paste0(
        "Violin_ribosomal_percent_by_",
        batch_col,
        "_",
        dataset_name,
        ".png"
      )
    ),
    width = 12,
    height = 6
  )
  
  invisible(list(ridge = p1, umap = p2, violin = p3))
}
