# ============================================================
# 02_metrics_math.R
# Métriques mathématiques : silhouette (approximée) + ARI bootstrap
# ============================================================

# ── 1. Silhouette approximée ──────────────────────────────────────────────────
#
# dist() sur l'ensemble complet est O(n²) en mémoire → infaisable au-delà de
# ~5 000 cellules. On sous-échantillonne de manière stratifiée (proportions
# de cluster préservées) pour que les petits clusters ne soient pas ignorés.
#
# @param emb       Matrice (n_cells × n_dims) de l'embedding (ex: harmony)
# @param clusters  Vecteur entier de labels de clusters (même longueur que nrow(emb))
# @param n_sample  Nombre max de cellules à échantillonner
# @param seed      Graine pour reproductibilité
# @return          Liste : $sil_df (data.frame avec sil et cluster),
#                          $mean_global, $prop_negative,
#                          $per_cluster (data.frame mean_sil par cluster),
#                          $idx (indices échantillonnés — pour annoter l'objet Seurat)
# ─────────────────────────────────────────────────────────────────────────────
compute_silhouette <- function(emb, clusters, n_sample = 3000, seed = 42) {
  set.seed(seed)
  n_cells <- nrow(emb)

  # Échantillonnage stratifié par cluster
  cluster_factor <- factor(clusters)
  idx <- unlist(lapply(levels(cluster_factor), function(cl) {
    cell_idx <- which(clusters == as.integer(cl))
    n_draw   <- max(1, round(n_sample * length(cell_idx) / n_cells))
    sample(cell_idx, min(length(cell_idx), n_draw))
  }))
  idx <- sort(unique(idx))

  sil <- cluster::silhouette(clusters[idx], dist(emb[idx, ]))

  sil_df <- data.frame(
    sil     = sil[, "sil_width"],
    cluster = factor(clusters[idx])
  )

  per_cluster <- sil_df |>
    dplyr::group_by(cluster) |>
    dplyr::summarise(
      mean_sil  = mean(sil),
      sd_sil    = sd(sil),
      prop_neg  = mean(sil < 0),
      n_cells   = dplyr::n(),
      .groups   = "drop"
    )

  list(
    sil_df       = sil_df,
    mean_global  = mean(sil_df$sil),
    prop_negative = mean(sil_df$sil < 0),
    per_cluster  = per_cluster,
    idx          = idx          # pour annoter seu_tmp$silhouette
  )
}


# ── 2. ARI bootstrap ─────────────────────────────────────────────────────────
#
# Principe : on rééchantillonne les cellules avec remise, on re-clusterise
# sur le sous-échantillon, on compare les labels avec l'ARI.
#
# CORRECTION vs code original :
#   - Le bootstrap utilisait `resolution` (vecteur entier) au lieu de `res`
#     (résolution courante) → les clusters bootstrap ne correspondaient pas.
#   - Le WNN est recalculé sur le sous-objet bootstrap : c'est correct,
#     mais coûteux. On conserve ce comportement car c'est le seul moyen
#     d'évaluer la stabilité du graph + clustering ensemble.
#   - On mappe les labels bootstrap → labels originaux via la cellule de
#     référence (idx_boot contient des doublons : on prend l'union unique).
#
# @param seu          Objet Seurat original (avec réductions harmony et apca)
# @param dims_list    list(dims_harmony, dims_apca)
# @param k            Nombre de voisins
# @param res          Résolution SCALAIRE
# @param ref_clusters Vecteur de labels de référence (entiers, n_cells)
# @param n_bootstrap  Nombre d'itérations bootstrap
# @param seed         Graine de base (incrémentée à chaque itération)
# @return             Liste : $scores (vecteur ARI), $mean, $sd
# ─────────────────────────────────────────────────────────────────────────────
compute_ari_bootstrap <- function(seu, dims_list, k, res, ref_clusters,
                                  n_bootstrap = 5, seed = 42) {
  n_cells <- ncol(seu)

  ari_scores <- vapply(seq_len(n_bootstrap), function(i) {
    set.seed(seed + i)
    idx_boot <- sample(n_cells, n_cells, replace = TRUE)
    seu_boot <- seu[, idx_boot]

    # Recalcul WNN + clustering sur le sous-objet bootstrap
    tryCatch({
      seu_boot <- FindMultiModalNeighbors(
        seu_boot,
        reduction.list = list("harmony", "apca"),
        dims.list      = dims_list,
        k.nn           = k,
        verbose        = FALSE
      )
      seu_boot <- FindClusters(
        seu_boot,
        resolution = res,          # ← résolution scalaire (bug corrigé)
        graph.name = "wsnn",
        algorithm  = 4,
        verbose    = FALSE
      )

      # Comparaison : labels originaux[idx_boot] vs labels bootstrap
      mclust::adjustedRandIndex(
        ref_clusters[idx_boot],
        as.integer(Idents(seu_boot))
      )
    }, error = function(e) {
      warning(sprintf("Bootstrap iter %d échouée : %s", i, conditionMessage(e)))
      NA_real_
    })
  }, numeric(1))

  ari_clean <- ari_scores[!is.na(ari_scores)]

  list(
    scores = ari_scores,
    mean   = if (length(ari_clean) > 0) mean(ari_clean) else NA_real_,
    sd     = if (length(ari_clean) > 0) sd(ari_clean)   else NA_real_
  )
}


# ── 3. Score composite mathématique ──────────────────────────────────────────
#
# Combine silhouette globale et ARI en un score unique normalisé [0, 1].
# Paramètres de pondération ajustables selon la priorité (stabilité vs séparation).
#
# @param mean_sil   Silhouette moyenne globale ∈ [-1, 1]
# @param mean_ari   ARI moyen bootstrap ∈ [0, 1]
# @param prop_neg   Proportion de cellules avec silhouette < 0 ∈ [0, 1]
# @param w_sil      Poids silhouette (défaut 0.4)
# @param w_ari      Poids ARI (défaut 0.4)
# @param w_neg      Poids pénalité prop_negative (défaut 0.2)
# @return           Score composite ∈ [0, 1]
# ─────────────────────────────────────────────────────────────────────────────
compute_math_score <- function(mean_sil, mean_ari, prop_neg,
                               w_sil = 0.4, w_ari = 0.4, w_neg = 0.2) {
  # Normalise silhouette de [-1,1] vers [0,1]
  sil_norm <- (mean_sil + 1) / 2
  # ARI est déjà dans [0,1]
  # Pénalité : 1 - prop_neg (moins de cellules mal assignées = mieux)
  score <- w_sil * sil_norm + w_ari * mean_ari + w_neg * (1 - prop_neg)
  round(score, 4)
}
