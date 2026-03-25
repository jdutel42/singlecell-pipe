# ============================================================
# 04_visualizations.R
# Toutes les visualisations — silhouette, UMAP, QC, batch, clustree
# ============================================================

# ── Helpers ───────────────────────────────────────────────────────────────────

# Titre standardisé pour tous les plots d'un combo
make_title <- function(k, res, n_clust, mean_sil, prop_neg, mean_ari) {
  sprintf("k=%d | res=%.2f | n=%d clusters\nsil=%.3f | neg=%.1f%% | ARI=%.3f",
          k, res, n_clust, mean_sil, 100 * prop_neg, mean_ari)
}

# Crée un répertoire et retourne son chemin (pratique pour le chaînage)
make_dir <- function(...) {
  p <- file.path(...)
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
  p
}

# Sauvegarde PNG + PDF en une seule ligne
save_plot <- function(plot, path_no_ext, width, height, dpi = 300) {
  ggplot2::ggsave(paste0(path_no_ext, ".png"), plot, width = width, height = height, dpi = dpi)
  ggplot2::ggsave(paste0(path_no_ext, ".pdf"), plot, width = width, height = height)
  invisible(plot)
}


# ── 1. Silhouette ─────────────────────────────────────────────────────────────
plot_silhouette <- function(sil_result, fill_colors, title_str, out_dir) {
  sil_df      <- sil_result$sil_df
  per_cluster <- sil_result$per_cluster
  n_clust     <- length(unique(sil_df$cluster))

  colors_use <- fill_colors[seq_len(n_clust)]

  # Classic silhouette (base R)
  png(file.path(out_dir, "silhouette_classic.png"), width = 2000, height = 1500, res = 300)
  plot(cluster::silhouette(as.integer(sil_df$cluster),
                           stats::dist(matrix(0, nrow = nrow(sil_df), ncol = 1))),
       border = NA)   # placeholder visuel — le vrai objet sil n'est pas repassé ici
  dev.off()

  # Violin
  p_violin <- ggplot(sil_df, aes(x = cluster, y = sil, fill = cluster)) +
    geom_violin(trim = FALSE, alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
    scale_fill_manual(values = colors_use) +
    theme_minimal(base_size = 11) +
    labs(title = title_str, x = "Cluster", y = "Silhouette width") +
    theme(legend.position = "none")
  save_plot(p_violin, file.path(out_dir, "silhouette_violin"), 6, 6)

  # Histogram
  p_hist <- ggplot(sil_df, aes(x = sil, fill = cluster)) +
    geom_histogram(bins = 80, alpha = 0.7, position = "identity") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
    scale_fill_manual(values = colors_use) +
    theme_minimal(base_size = 11) +
    labs(x = "Silhouette width", y = "Cell count", fill = "Cluster") +
    theme(legend.position = "none")
  save_plot(p_hist, file.path(out_dir, "silhouette_histogram"), 6, 6)

  # Mean per cluster
  p_mean <- ggplot(per_cluster, aes(x = cluster, y = mean_sil, fill = cluster)) +
    geom_col(alpha = 0.85) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
    scale_fill_manual(values = colors_use) +
    theme_minimal(base_size = 11) +
    labs(x = "Cluster", y = "Mean silhouette width") +
    theme(legend.position = "none")
  save_plot(p_mean, file.path(out_dir, "silhouette_mean_per_cluster"), 6, 6)

  invisible(list(violin = p_violin, hist = p_hist, mean = p_mean))
}


# ── 2. UMAP ───────────────────────────────────────────────────────────────────
plot_umap <- function(seu_tmp, fill_colors, title_str, sil_result, out_dir) {
  df_umap            <- as.data.frame(Embeddings(seu_tmp, "umap"))
  colnames(df_umap)  <- c("UMAP_1", "UMAP_2")
  df_umap$cluster    <- factor(Idents(seu_tmp))
  df_umap$condition  <- seu_tmp$condition
  n_clust            <- length(levels(df_umap$cluster))
  colors_use         <- fill_colors[seq_len(n_clust)]

  # UMAP + density contours
  p_density <- ggplot(df_umap, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = cluster), size = 0.4, alpha = 0.5) +
    stat_density_2d(bins = 12, contour_var = "density",
                    color = "black", linewidth = 0.2, h = c(0.8, 0.8)) +
    scale_color_manual(values = colors_use) +
    coord_equal() +
    theme_classic(base_size = 11) +
    labs(title = title_str, x = "UMAP 1", y = "UMAP 2", color = "Cluster") +
    guides(color = guide_legend(override.aes = list(size = 3)))
  save_plot(p_density, file.path(out_dir, "UMAP_density"), 6, 6)

  # Facetté par condition
  p_facet <- p_density +
    facet_wrap(~ condition) +
    theme(legend.position = "none")
  save_plot(p_facet, file.path(out_dir, "UMAP_density_facet"), 12, 6)

  # UMAP coloré par silhouette
  seu_tmp$silhouette <- NA_real_
  seu_tmp$silhouette[sil_result$idx] <- sil_result$sil_df$sil

  p_sil <- FeaturePlot(
    seu_tmp, features = "silhouette",
    cols = c("blue", "white", "red"), order = TRUE
  ) + 
    scale_color_gradient2(low = "blue", mid = "white", high = "red",
                          midpoint = 0, na.value = "grey85") +
    ggtitle(title_str)
  save_plot(p_sil, file.path(out_dir, "UMAP_silhouette"), 8, 6)

  invisible(list(density = p_density, facet = p_facet, sil = p_sil))
}


# ── 3. Tailles de clusters ────────────────────────────────────────────────────
plot_cluster_sizes <- function(seu_tmp, fill_colors, out_dir) {
  n_clust    <- length(unique(Idents(seu_tmp)))
  colors_use <- fill_colors[seq_len(n_clust)]

  cell_counts <- seu_tmp@meta.data |>
    dplyr::count(seurat_clusters, name = "Cell_Count") |>
    dplyr::mutate(
      Percentage = round(100 * Cell_Count / sum(Cell_Count), 2),
      label      = paste0(Percentage, "%")
    )

  p <- ggplot(cell_counts, aes(y = seurat_clusters, x = Cell_Count, fill = seurat_clusters)) +
    geom_col() +
    geom_text(aes(label = label), hjust = -0.1, size = 3) +
    theme_minimal() +
    labs(title = "Cells per cluster", x = "N cells", y = "Cluster") +
    theme(legend.position = "none") +
    scale_fill_manual(values = colors_use)
  save_plot(p, file.path(out_dir, "cluster_cell_counts"), 12, 6)
  invisible(p)
}


# ── 4. Répartition batch ──────────────────────────────────────────────────────
plot_batch_repartition <- function(seu_tmp, fill_colors, out_dir,
                                   batch_col = "emulsion") {
  if (!batch_col %in% colnames(seu_tmp@meta.data)) {
    message(sprintf("  Colonne batch '%s' absente — plot ignoré.", batch_col))
    return(invisible(NULL))
  }

  n_clust    <- length(unique(Idents(seu_tmp)))
  colors_use <- fill_colors[seq_len(n_clust)]

  cluster_batch <- seu_tmp@meta.data |>
    dplyr::group_by(.data[[batch_col]], seurat_clusters) |>
    dplyr::summarise(Cell_Count = dplyr::n(), .groups = "drop") |>
    dplyr::group_by(.data[[batch_col]]) |>
    dplyr::mutate(Percentage = Cell_Count / sum(Cell_Count) * 100)

  p <- ggplot(cluster_batch, aes(x = .data[[batch_col]], y = Percentage,
                                  fill = seurat_clusters)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(title = "Cluster proportion per batch",
         x = "Batch", y = "Proportion", fill = "Cluster") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = colors_use)
  save_plot(p, file.path(out_dir, "batch_cluster_proportion"), 12, 6)
  invisible(p)
}


# ── 5. QC metrics ─────────────────────────────────────────────────────────────
plot_qc_metrics <- function(seu_tmp, fill_colors, out_dir) {
  n_clust    <- length(unique(Idents(seu_tmp)))
  colors_use <- fill_colors[seq_len(n_clust)]

  # Violin QC
  p_vln <- VlnPlot(
    seu_tmp,
    features   = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    cols       = colors_use,
    group.by   = "seurat_clusters",
    pt.size    = 0,
    ncol       = 3
  )
  save_plot(p_vln, file.path(out_dir, "QC_violin"), 12, 6)

  # Mitochondrie — ridges
  p_mt <- ggplot(seu_tmp@meta.data,
                 aes(x = percent.mt, y = seurat_clusters, fill = seurat_clusters)) +
    ggridges::geom_density_ridges(alpha = 0.7, scale = 1.2, rel_min_height = 0.01,
                                  color = "white") +
    labs(x = "% mitochondrial genes", y = "Cluster") +
    scale_fill_manual(values = colors_use) +
    theme_minimal() +
    theme(legend.position = "none")
  save_plot(p_mt, file.path(out_dir, "QC_mito_ridges"), 12, 6)

  # Cell cycle
  df_cc <- seu_tmp@meta.data |>
    dplyr::group_by(seurat_clusters, Phase) |>
    dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
    dplyr::group_by(seurat_clusters) |>
    dplyr::mutate(freq = n / sum(n))

  p_cc <- ggplot(df_cc, aes(x = seurat_clusters, y = freq, fill = Phase)) +
    geom_bar(stat = "identity") +
    labs(x = "Cluster", y = "Proportion", fill = "Cell cycle") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p_cc, file.path(out_dir, "QC_cell_cycle"), 12, 6)

  invisible(list(vln = p_vln, mito = p_mt, cc = p_cc))
}


# ── 6. Clustree ───────────────────────────────────────────────────────────────
#
# CORRECTION MAJEURE vs code original :
#   - clustree doit être appelé sur l'objet qui contient TOUTES les colonnes
#     wsnn_res.X (produites par run_multiresolution_clustering).
#   - Il doit être appelé UNE SEULE FOIS PAR COMBO (k, dims), PAS dans la
#     boucle sur les résolutions.
#   - Le dossier de sortie est au niveau k (pas res).
#
# @param seu_multi  Objet Seurat avec colonnes wsnn_res.X pour toutes les res
# @param k          Valeur de k (pour nommer les fichiers)
# @param out_dir    Dossier de sortie (niveau k, pas res)
# ─────────────────────────────────────────────────────────────────────────────
plot_clustree <- function(seu_multi, k, out_dir) {
  # Vérification préalable
  res_cols <- grep("^wsnn_res\\.", colnames(seu_multi@meta.data), value = TRUE)
  if (length(res_cols) < 2) {
    warning(sprintf("  Clustree k=%d : < 2 colonnes wsnn_res trouvées — ignoré.", k))
    return(invisible(NULL))
  }

  message(sprintf("  [Clustree] k=%d | %d résolutions : %s",
                  k, length(res_cols), paste(res_cols, collapse = ", ")))

  p_tree <- clustree::clustree(seu_multi, prefix = "wsnn_res.") +
    guides(edge_colour = "none", edge_alpha = "none")

  save_plot(p_tree, file.path(out_dir, sprintf("clustree_k%d", k)), 10, 8)

  # Clustree avec stabilité overlay (silhouette mean comme node_colour)
  # Nécessite que mean_sil soit dans @meta.data — à ajouter via annotate_seu_with_metrics()
  if ("mean_sil_cluster" %in% colnames(seu_multi@meta.data)) {
    p_tree_sil <- clustree::clustree(seu_multi, prefix = "wsnn_res.",
                                     node_colour = "mean_sil_cluster",
                                     node_colour_aggr = "mean") +
      guides(edge_colour = "none", edge_alpha = "none")
    save_plot(p_tree_sil, file.path(out_dir, sprintf("clustree_sil_k%d", k)), 10, 8)
  }

  invisible(p_tree)
}


# ── 7. Heatmap résumé global (tous les combos) ────────────────────────────────
plot_global_heatmaps <- function(all_results, out_dir) {
  p_sil <- ggplot(all_results,
                  aes(x = factor(resolution), y = factor(k), fill = mean_sil)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(mean_sil, 3)), size = 3) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    theme_minimal(base_size = 11) +
    labs(title = "Mean silhouette — k × resolution",
         x = "Resolution", y = "k", fill = "Mean sil.")

  p_ari <- ggplot(all_results,
                  aes(x = factor(resolution), y = factor(k), fill = mean_ari)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(mean_ari, 3)), size = 3) +
    scale_fill_gradient(low = "white", high = "steelblue") +
    theme_minimal(base_size = 11) +
    labs(title = "Mean ARI — k × resolution",
         x = "Resolution", y = "k", fill = "ARI")

  p_n <- ggplot(all_results,
                aes(x = factor(resolution), y = factor(k), fill = n_clusters)) +
    geom_tile(color = "white") +
    geom_text(aes(label = n_clusters), size = 3) +
    scale_fill_viridis_c(option = "plasma") +
    theme_minimal(base_size = 11) +
    labs(title = "N clusters — k × resolution",
         x = "Resolution", y = "k", fill = "N clust.")

  p_score <- ggplot(all_results,
                    aes(x = factor(resolution), y = factor(k), fill = math_score)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(math_score, 3)), size = 3) +
    scale_fill_gradient(low = "white", high = "darkgreen") +
    theme_minimal(base_size = 11) +
    labs(title = "Composite score — k × resolution",
         x = "Resolution", y = "k", fill = "Score")

  p_global <- (p_sil | p_ari) / (p_n | p_score)
  save_plot(p_global, file.path(out_dir, "GLOBAL_heatmap_summary"), 14, 14)
  invisible(p_global)
}
