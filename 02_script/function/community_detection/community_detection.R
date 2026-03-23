# ============================================================
# Clustering parameter optimization — WNN multimodal
# ============================================================

library(Seurat)
library(furrr)
library(future)
library(cluster)
library(ggplot2)
library(dplyr)
library(patchwork)
library(clustree)
library(mclust)   # pour ARI

# --- Plan parallèle : multisession (workers isolés, Seurat-safe)
# plan(multisession, workers = min(parallel::detectCores() - 1, 6))
plan(multisession, workers = 2)  # For testing with fewer workers


# ============================================================
# Main function : run clustering + silhouette + ARI bootstrap
# ============================================================

run_combo <- function(params, seu, dims_list, dims_use, base_dir, n_bootstrap = 3) {
  
  
  # ============================================================
  # Parameters
  # ============================================================
  
  k          <- params$k
  resolution <- params$resolution
  
  message(sprintf("[k=%d | res=%.2f] Starting...", k, resolution))
  
  # ============================================================
  # Run WNN + clustering + UMAP
  # ============================================================
  
  # --- 1. WNN + clustering ---
  seu_tmp <- 
  FindMultiModalNeighbors(
    seu,
    reduction.list = list("harmony", "apca"),
    dims.list      = dims_list,
    k.nn           = k,
    verbose        = FALSE
  )
  
  for (res in resolution) {
    seu_tmp <- FindClusters(
      seu_tmp,
      resolution = res,
      graph.name = "wsnn",
      algorithm  = 4,   # Leiden (More stable than Louvain)
      verbose    = FALSE
    )

    RunUMAP(
      seu_tmp,
      reduction = "harmony",
      dims      = dims_use,
      verbose   = FALSE
    )
    
    # ============================================================
    # Calcul metrics
    # ============================================================
    
    clust    <- as.integer(Idents(seu_tmp))
    n_clust  <- length(unique(clust))
    emb      <- Embeddings(seu_tmp, "harmony")[, dims_use]
    
    # --- 2. Approximate silhouette (sampling to avoid dist() on full dataset (very long)) ---
    n_cells  <- nrow(emb)
    n_sample <- min(n_cells, 3000)
    idx      <- sample(n_cells, n_sample)
    
    sil      <- silhouette(clust[idx], dist(emb[idx, ]))
    sil_df   <- data.frame(
      sil     = sil[, "sil_width"],
      cluster = factor(clust[idx])
    )
    
    mean_sil_global <- mean(sil_df$sil)
    prop_neg        <- mean(sil_df$sil < 0)
    
    mean_sil_per_cluster <- sil_df %>%
      group_by(cluster) %>%
      summarise(mean_sil = mean(sil), .groups = "drop")
    
    # --- 3. ARI stability by bootstrap ---
    ari_scores <- vapply(seq_len(n_bootstrap), function(i) {
      idx_boot  <- sample(n_cells, n_cells, replace = TRUE)
      seu_boot  <- seu[, idx_boot]
      seu_boot  <- FindMultiModalNeighbors(seu_boot,
                                           reduction.list = list("harmony", "apca"),
                                           dims.list = dims_list, k.nn = k, verbose = FALSE)
      seu_boot  <- FindClusters(seu_boot, resolution = resolution,
                                graph.name = "wsnn", algorithm = 4, verbose = FALSE)
      mclust::adjustedRandIndex(clust[idx_boot], as.integer(Idents(seu_boot)))
    }, numeric(1))
    
    mean_ari <- mean(ari_scores)
    
    # ============================================================
    # Visualize
    # ============================================================
    
    # Create dir 
    k_dir <- file.path(base_dir, paste0("Neighbor_k", paste0(k)))
    dir.create(k_dir, recursive = TRUE, showWarnings = FALSE)
    
    res_dir <- file.path(k_dir, paste0("Resolution_", sprintf("%.2f", resolution)))
    dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
    # --- 4. Outputs ---
    n_colors   <- max(n_clust, 1)
    fill_colors <- if (n_clust <= length(pal)) pal[1:n_clust] else colors_55[1:n_clust]
    
    # param_name <- sprintf("k%d_res%.2f_dims%d-%d", k, res, min(dims_use), max(dims_use))
    # param_dir  <- file.path(base_dir, param_name)
    # dir.create(param_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Common title for all plots
    title_str <- sprintf("k=%d | res=%.2f | n_clust=%d\nmean_sil=%.3f | prop_neg=%.3f | ARI=%.3f",
                         k, res, n_clust, mean_sil_global, prop_neg, mean_ari)
    
    # ============================================================
    # Silouhette
    # ============================================================
    dir.create(file.path(res_dir, "Silhouette"), recursive = TRUE, showWarnings = FALSE)
    
    # Silhouette classique
    png(file.path(res_dir, "Silhouette", "silhouette_classic.png"),
        width = 2000, height = 1500, res = 300)
    plot(sil, border = NA)
    dev.off()
    
    # 4a. Violin silhouette
    p_violin <- ggplot(sil_df, aes(x = cluster, y = sil, fill = cluster)) +
      geom_violin(trim = FALSE, alpha = 0.8) +
      geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
      scale_fill_manual(values = fill_colors) +
      theme_minimal(base_size = 11) +
      labs(title = title_str, x = "Cluster", y = "Silhouette width") +
      theme(legend.position = "none")
    
    ggsave(file.path(res_dir, "Silhouette", "silhouette_violin.png"),
           p_violin, width = 6, height = 6, dpi = 300)
    
    ggsave(file.path(res_dir, "Silhouette", "silhouette_violin.pdf"),
           p_violin, width = 6, height = 6)
    
    # 4b. Histogram
    p_hist <- ggplot(sil_df, aes(x = sil, fill = cluster)) +
      geom_histogram(bins = 80, alpha = 0.7, position = "identity") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
      scale_fill_manual(values = fill_colors) +
      theme_minimal(base_size = 11) +
      labs(x = "Silhouette width", y = "Cell count", fill = "Cluster") +
      theme(legend.position = "none")
    
    ggsave(file.path(res_dir, "Silhouette", "silhouette_histogram.png"),
           p_hist, width = 6, height = 6, dpi = 300)
    
    ggsave(file.path(res_dir, "Silhouette", "silhouette_histogram.pdf"),
           p_hist, width = 6, height = 6)
    
    # 4c. Mean silhouette per cluster
    p_mean <- ggplot(mean_sil_per_cluster, aes(x = cluster, y = mean_sil, fill = cluster)) +
      geom_col(alpha = 0.85) +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
      scale_fill_manual(values = fill_colors) +
      theme_minimal(base_size = 11) +
      labs(x = "Cluster", y = "Mean silhouette width") +
      theme(legend.position = "none")
    
    ggsave(file.path(res_dir, "Silhouette", "silhouette_mean_per_cluster.png"),
           p_mean, width = 6, height = 6, dpi = 300)
    
    ggsave(file.path(res_dir, "Silhouette", "silhouette_mean_per_cluster.pdf"),
           p_mean, width = 6, height = 6)
    
    # 4d. UMAP silhouette (sur sous-échantillon)
    seu_tmp$silhouette <- NA_real_
    seu_tmp$silhouette[idx] <- sil[, "sil_width"]
    
    p_umap_sil <- FeaturePlot(
      seu_tmp, features = "silhouette",
      cols = c("blue", "white", "red"), order = TRUE, na.value = "grey85"
    ) + ggtitle(title_str)
    
    ggsave(file.path(res_dir, "Silhouette", "UMAP_silhouette.png"),
           p_umap_sil, width = 8, height = 6, dpi = 300)
    
    ggsave(file.path(res_dir, "Silhouette", "UMAP_silhouette.pdf"),
           p_umap_sil, width = 8, height = 6)
    
    # 4e. Combine all silhouette plots
    p_combo <- (p_violin | p_hist) / (p_mean | p_umap_density)
    
    ggsave(file.path(res_dir, "Silhouette", "summary_combo.png"),
           p_combo, width = 14, height = 10, dpi = 300)
    
    ggsave(file.path(res_dir, "Silhouette", "summary_combo.pdf"),
           p_combo, width = 14, height = 10)
    
    # ============================================================
    # Clusters on UMAP + density
    # ============================================================
    dir.create(file.path(res_dir, "UMAP_cluster"), recursive = TRUE, showWarnings = FALSE)
    
    # 4e. UMAP density per cluster
    df_umap <- as.data.frame(Embeddings(seu_tmp, "umap"))
    colnames(df_umap) <- c("UMAP_1", "UMAP_2")
    df_umap$cluster   <- factor(clust)
    df_umap$condition <- seu_tmp$condition
    
    p_umap_density <- ggplot(df_umap, aes(UMAP_1, UMAP_2)) +
      geom_point(aes(color = cluster), size = 0.4, alpha = 0.5) +
      stat_density_2d(bins = 12, contour_var = "density",
                      color = "black", linewidth = 0.2, h = c(0.8, 0.8)) +
      scale_color_manual(values = fill_colors) +
      coord_equal() +
      theme_classic(base_size = 11) +
      labs(title = title_str, x = "UMAP 1", y = "UMAP 2", color = "Cluster") +
      guides(color = guide_legend(override.aes = list(size = 3)))
    
    ggsave(file.path(res_dir, "UMAP_cluster", "UMAP_density.png"),
           p_umap_density, width = 6, height = 6, dpi = 300)
    
    ggsave(file.path(res_dir, "UMAP_cluster", "UMAP_density.pdf"),
           p_umap_density, width = 6, height = 6)
    
    p_umap_density_facet <- ggplot(df_umap, aes(UMAP_1, UMAP_2)) +
      geom_point(aes(color = cluster), size = 0.4, alpha = 0.5) +
      stat_density_2d(bins = 12, contour_var = "density",
                      color = "black", linewidth = 0.2, h = c(0.8, 0.8)) +
      scale_color_manual(values = fill_colors) +
      coord_equal() +
      theme_classic(base_size = 11) +
      labs(title = title_str, x = "UMAP 1", y = "UMAP 2", color = "Cluster") +
      guides(color = guide_legend(override.aes = list(size = 3))) +
      facet_wrap(~ condition) +
      theme(legend.position = "none")
    
    ggsave(file.path(res_dir, "UMAP_cluster", "UMAP_density_facet.png"),
           p_umap_density_facet, width = 12, height = 6, dpi = 300)
    
    ggsave(file.path(res_dir, "UMAP_cluster", "UMAP_density_facet.pdf"),
           p_umap_density_facet, width = 12, height = 6)
    
    # ============================================================
    # Number of cells per cluster
    # ============================================================
    dir.create(file.path(res_dir, "Cluster_sizes"), recursive = TRUE, showWarnings = FALSE)
    
    # Plot number of cells per cluster with percentage labels of total cells
    cell_counts <- seu_tmp@meta.data %>%
      count(seurat_clusters, name = "Cell_Count") %>%
      mutate(
        Percentage = round(100 * Cell_Count / sum(Cell_Count), 6),
        label = paste0(Percentage, "%")
      )
    
    
    ggplot(cell_counts, aes(y = seurat_clusters, x = Cell_Count, fill = seurat_clusters)) +
      geom_col() +
      geom_text(aes(label = label), hjust = -0.1, size = 3) +
      theme_minimal() +
      labs(
        title = "Number (and %) of cells per cluster",
        x = "Number of cells",
        y = "Cluster"
      ) +
      theme(legend.position = "none") +
      scale_fill_manual(values = fill_colors)
    
    # Save the plot
    ggsave(
      filename = file.path(
        res_dir, 
        "Cluster_sizes",
        "Clusters_cell_counts.png"),
      width = 12,
      height = 6,
      dpi = 300)
    
    ggsave(
      filename = file.path(
        res_dir, 
        "Cluster_sizes",
        "Clusters_cell_counts.pdf"),
      width = 12,
      height = 6)
  
    
    
    # ============================================================
    # Batch repartition
    # ============================================================
    dir.create(file.path(res_dir, "Batch_repartition"), recursive = TRUE, showWarnings = FALSE)
    
    # Proportion of cells per cluster for each batch
    cluster_patient <- seu_tmp@meta.data %>%
      group_by(emulsion, seurat_clusters) %>%
      summarise(Cell_Count = n()) %>%
      group_by(emulsion) %>%
      mutate(Percentage = Cell_Count / sum(Cell_Count) * 100)
    
    ggplot(cluster_patient, aes(x = emulsion, y = Percentage, fill = seurat_clusters)) +
      geom_bar(stat = "identity", position = "fill") +
      scale_y_continuous(labels = scales::percent_format()) +
      labs(
        title = "Proportion of cells per cluster for each batch",
        x = "Batch",
        y = "Proportion of cells",
        fill = "Cluster"
      ) + 
      theme(
        axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   vjust = 1)) +
      scale_fill_manual(values = fill_colors)
    
    # Save the plot
    ggsave(
      filename = file.path(
        res_dir,
        "Batch_repartition",
        "Clusters_repartition_per_batch.png"),
      width = 12,
      height = 6,
      dpi = 300)
    
    ggsave(
      filename = file.path(
        res_dir,
        "Batch_repartition",
        "Clusters_repartition_per_batch.pdf"),
      width = 12,
      height = 6)
    
    # ============================================================
    # QC metrics per cluster
    # ============================================================
    dir.create(file.path(res_dir, "QC_metrics"), recursive = TRUE, showWarnings = FALSE)
    
    # All QC metrics (nFeature_RNA, nCount_RNA, percent.mt) by cluster
    VlnPlot(
      seu_tmp,
      features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
      cols = fill_colors,
      group.by = "seurat_clusters",
      pt.size = 0,
      ncol = 3)
    
    ggsave(
      filename = file.path(
        res_dir,
        "QC_metrics",
        "Clusters_QC_metrics.png"),
      width = 12,
      height = 6,
      dpi = 300)
    
    ggsave(
      filename = file.path(
        res_dir,
        "QC_metrics",
        "Clusters_QC_metrics.pdf"),
      width = 12,
      height = 6)
    
    stats_mt <- seu_tmp@meta.data %>%
      group_by(seurat_clusters) %>%
      summarise(
        mean_pct_mt   = mean(percent.mt),
        median_pct_mt = median(percent.mt),
        sd_pct_mt     = sd(percent.mt)
      )
    
    kable(stats_mt, format = "html") %>%
      kable_styling()
    
    write.table(stats_mt, 
                file = file.path(
                  res_dir,
                  "QC_metrics",
                  "cluster_mitochondrial_content_stats.tsv"), 
                row.names = FALSE, 
                sep = "\t")
    
    # Mitochondia
    ggplot(seu_tmp@meta.data,
           aes(
             x = percent.mt,
             y = seurat_clusters,
             fill = seurat_clusters
           )) +
      geom_density_ridges(
        alpha = 0.7,
        scale = 1.2,
        rel_min_height = 0.01,
        color = "white") +
      # theme_ridges() +
      # theme(
      #   legend.position = "none"
      # ) +
      # geom_vline(xintercept = 20, linetype = "dashed", color = "red") +
      labs(
        x = "% mitochondrial genes",
        y = "Clusters") +
      scale_fill_manual(values = fill_colors)
    
    ggsave(
      filename = file.path(
        res_dir,
        "QC_metrics",
        "Clusters_mitochondrial_content_ridges.png"),
      width = 12,
      height = 6,
      dpi = 300)
    
    ggsave(
      filename = file.path(
        res_dir,
        "QC_metrics",
        "Clusters_mitochondrial_content_ridges.pdf"),
      width = 12,
      height = 6)
    
    # Cell cycle phase
    df <- seu_tmp@meta.data %>%
      group_by(seurat_clusters, Phase) %>%
      summarise(n = n()) %>%
      mutate(freq = n / sum(n))
    
    ggplot(df, aes(x = seurat_clusters, y = freq, fill = Phase)) +
      geom_bar(stat = "identity") +
      labs(x = "Cluster", y = "Proportion of cells", fill = "Cell cycle phase") +
      theme_minimal() +
      ggtitle("Distribution of cell cycle phases per cluster") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(
      filename = file.path(
        res_dir,
        "QC_metrics",
        "Clusters_cell_cycle_phase_distribution.png"),
      width = 12,
      height = 6,
      dpi = 300)
    
    ggsave(
      filename = file.path(
        res_dir,
        "QC_metrics",
        "Clusters_cell_cycle_phase_distribution.pdf"),
      width = 12,
      height = 6)
    
    ggplot(df, aes(x = seurat_clusters, y = n, fill = Phase)) +
      geom_bar(stat = "identity") +
      labs(x = "Cluster", y = "Number of cells", fill = "Cell cycle phase") +
      theme_minimal() +
      ggtitle("Cell counts per cell cycle phase and cluster") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(
      filename = file.path(
        res_dir,
        "QC_metrics",
        "Clusters_cell_cycle_phase_counts.png"),
      width = 12,
      height = 6, dpi = 300)
    
    ggsave(
      filename = file.path(
        res_dir,
        "QC_metrics",
        "Clusters_cell_cycle_phase_counts.pdf"),
      width = 12,
      height = 6)
    
    
    mat <- df %>%
      select(seurat_clusters, Phase, freq) %>%
      pivot_wider(names_from = Phase, values_from = freq, values_fill = 0) %>%
      column_to_rownames("seurat_clusters")
    
    pheatmap(as.matrix(mat),
             cluster_rows = TRUE,
             cluster_cols = FALSE,
             display_numbers = TRUE,
             main = "Cell cycle phase proportions per cluster")
    
    ggsave(
      filename = file.path(
        res_dir,
        "QC_metrics",
        "Clusters_cell_cycle_phase_heatmap.png"),
      width = 12,
      height = 6,
      dpi = 300)
    
    ggsave(
      filename = file.path(
        res_dir,
        "QC_metrics",
        "Clusters_cell_cycle_phase_heatmap.pdf"),
      width = 12,
      height = 6)
    
    DimPlot(
      seu_tmp,
      reduction = "umap",
      group.by = "seurat_clusters",
      split.by = "Phase",
      cols = fill_colors,
      pt.size = 0.4,
    ) + ggtitle("UMAP colored by cell cycle phase")
    
    ggsave(
      filename = file.path(
        res_dir,
        "QC_metrics",
        "UMAP_seurat_clusters_split_by_cell_cycle_phase.png"),
      width = 12,
      height = 6,
      dpi = 300)
    
    ggsave(
      filename = file.path(
        res_dir,
        "QC_metrics",
        "UMAP_seurat_clusters_split_by_cell_cycle_phase.pdf"),
      width = 12,
      height = 6)
    
    message(sprintf("[k=%d | res=%.2f] Done. mean_sil=%.3f | ARI=%.3f", k, res, mean_sil_global, mean_ari))
    
  }
  
  # ============================================================
  # Clustree
  # ============================================================
  dir.create(file.path(res_dir, paste0("Clustree_k", k)), recursive = TRUE, showWarnings = FALSE)
  
  # 4f. Clustree classique (all res tested for this object specific k)
  p_clustree <- clustree(seu_tmp, prefix = "wsnn_res.")
  
  ggsave(file.path(res_dir, "Clustree", "clustree.png"),
         p_clustree, width = 8, height = 6, dpi = 300)
  
  ggsave(file.path(res_dir, "Clustree", "clustree.pdf"),
         p_clustree, width = 8, height = 6)
  
  # p_clustree_signature <- clustree(seu_tmp, prefix = "wsnn_res.", node_colour = "mean_sil")
  
  
  # ============================================================
  # Annotation
  # ============================================================
  
  
  
  
  
  
  
  
  
  return(data.frame(
    k              = k,
    resolution     = resolution,
    n_clusters     = n_clust,
    mean_sil       = mean_sil_global,
    prop_negative  = prop_neg,
    mean_ari       = mean_ari,
    sd_ari         = sd(ari_scores)
  ))
}

# ============================================================
# Grid + lancement
# ============================================================

dims_list <- list(1:13, 1:13)
dims_use  <- 1:13

k_values   <- c(20, 30, 40)
res_values <- seq(0.2, 1.0, by = 0.2)
grid       <- expand.grid(k = k_values, resolution = res_values)

base_dir <- file.path(DIR_RESULT, "All_Tcells", "Community_detection", "Tuning_k_res")
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)

all_results <- future_map_dfr(
  split(grid, seq(nrow(grid))),
  run_combo,
  seu      = seu_path,
  dims_list = dims_list,
  dims_use = dims_use,
  base_dir = base_dir,
  .options = furrr_options(seed = TRUE)
)

# ============================================================
# Summary global — heatmap k × resolution
# ============================================================

# Heatmap mean silhouette
p_heatmap_sil <- ggplot(all_results,
                        aes(x = factor(resolution), y = factor(k), fill = mean_sil)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(mean_sil, 3)), size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal(base_size = 11) +
  labs(title = "Mean silhouette — k × resolution",
       x = "Resolution", y = "k", fill = "Mean sil.")

# Heatmap ARI
p_heatmap_ari <- ggplot(all_results,
                        aes(x = factor(resolution), y = factor(k), fill = mean_ari)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(mean_ari, 3)), size = 3) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal(base_size = 11) +
  labs(title = "Mean ARI (bootstrap) — k × resolution",
       x = "Resolution", y = "k", fill = "ARI")

# Heatmap n_clusters
p_heatmap_n <- ggplot(all_results,
                      aes(x = factor(resolution), y = factor(k), fill = n_clusters)) +
  geom_tile(color = "white") +
  geom_text(aes(label = n_clusters), size = 3) +
  scale_fill_viridis_c(option = "plasma") +
  theme_minimal(base_size = 11) +
  labs(title = "N clusters — k × resolution",
       x = "Resolution", y = "k", fill = "N clust.")

p_global <- p_heatmap_sil / p_heatmap_ari / p_heatmap_n

ggsave(file.path(base_dir, "GLOBAL_heatmap_summary.png"),
       p_global, width = 10, height = 14, dpi = 300)

write.csv(all_results,
          file.path(base_dir, "all_results_metrics.csv"),
          row.names = FALSE)