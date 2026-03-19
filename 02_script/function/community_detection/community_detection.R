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
plan(multisession, workers = min(parallel::detectCores() - 1, 6))

# ============================================================
# Main function : run clustering + silhouette + ARI bootstrap
# ============================================================

run_combo <- function(params, seu, dims_list, dims_use, base_dir, n_bootstrap = 3) {
  
  k          <- params$k
  resolution <- params$resolution
  
  message(sprintf("[k=%d | res=%.2f] Starting...", k, resolution))
  
  # --- 1. WNN + clustering ---
  seu_tmp <- 
  FindMultiModalNeighbors(
    seu,
    reduction.list = list("harmony", "apca"),
    dims.list      = dims_list,
    k.nn           = k,
    verbose        = FALSE
  ) %>% 
  FindClusters(
    seu_tmp,
    resolution = resolution,
    graph.name = "wsnn",
    algorithm  = 4,   # Leiden (More stable than Louvain)
    verbose    = FALSE
  ) %>% 
  RunUMAP(
    seu_tmp,
    reduction = "harmony",
    dims      = dims_use,
    verbose   = FALSE
  )
  
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
  
  # --- 4. Outputs ---
  n_colors   <- max(n_clust, 1)
  fill_colors <- if (n_clust <= length(pal)) pal[1:n_clust] else colors_55[1:n_clust]
  
  param_name <- sprintf("k%d_res%.2f_dims%d-%d", k, resolution, min(dims_use), max(dims_use))
  param_dir  <- file.path(base_dir, param_name)
  dir.create(param_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Common title for all plots
  title_str <- sprintf("k=%d | res=%.2f | n_clust=%d\nmean_sil=%.3f | prop_neg=%.3f | ARI=%.3f",
                       k, resolution, n_clust, mean_sil_global, prop_neg, mean_ari)
  
  # 4a. Violin silhouette
  p_violin <- ggplot(sil_df, aes(x = cluster, y = sil, fill = cluster)) +
    geom_violin(trim = FALSE, alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
    scale_fill_manual(values = fill_colors) +
    theme_minimal(base_size = 11) +
    labs(title = title_str, x = "Cluster", y = "Silhouette width") +
    theme(legend.position = "none")
  
  # 4b. Histogram
  p_hist <- ggplot(sil_df, aes(x = sil, fill = cluster)) +
    geom_histogram(bins = 80, alpha = 0.7, position = "identity") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
    scale_fill_manual(values = fill_colors) +
    theme_minimal(base_size = 11) +
    labs(x = "Silhouette width", y = "Cell count", fill = "Cluster") +
    theme(legend.position = "none")
  
  # 4c. Mean silhouette per cluster
  p_mean <- ggplot(mean_sil_per_cluster, aes(x = cluster, y = mean_sil, fill = cluster)) +
    geom_col(alpha = 0.85) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
    scale_fill_manual(values = fill_colors) +
    theme_minimal(base_size = 11) +
    labs(x = "Cluster", y = "Mean silhouette width") +
    theme(legend.position = "none")
  
  # 4d. UMAP silhouette (sur sous-échantillon)
  seu_tmp$silhouette <- NA_real_
  seu_tmp$silhouette[idx] <- sil[, "sil_width"]
  
  p_umap_sil <- FeaturePlot(
    seu_tmp, features = "silhouette",
    cols = c("blue", "white", "red"), order = TRUE, na.value = "grey85"
  ) + ggtitle(title_str)
  
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
  
  # 4f. Clustree (all res tested for this object)
  p_clustree <- clustree(seu_tmp, prefix = "wsnn_res.")
  
  # 4g. Save combined (patchwork)
  p_combo <- (p_violin | p_hist) / (p_mean | p_umap_density)
  ggsave(file.path(param_dir, "summary_combo.png"),
         p_combo, width = 14, height = 10, dpi = 300)
  
  ggsave(file.path(param_dir, "UMAP_silhouette.png"),
         p_umap_sil, width = 8, height = 6, dpi = 300)
  
  # Silhouette classique
  png(file.path(param_dir, "silhouette_classic.png"),
      width = 2000, height = 1500, res = 300)
  plot(sil, border = NA)
  dev.off()
  
  message(sprintf("[k=%d | res=%.2f] Done. mean_sil=%.3f | ARI=%.3f", k, resolution, mean_sil_global, mean_ari))
  
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

base_dir <- file.path(DIR_RESULT, "All_Tcells", "Clustering", "Tuning_k_res")
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)

all_results <- future_map_dfr(
  split(grid, seq(nrow(grid))),
  run_combo,
  seu      = seu_CD4_QC_Moderate,
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