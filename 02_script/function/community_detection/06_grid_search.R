# ============================================================
# 06_grid_search.R
# Orchestrateur principal du grid search
#
# Architecture du pipeline :
#
#   Pour chaque combo (k, dims) :
#     1. build_wnn_graph()           — 1 fois, coûteux
#     2. run_multiresolution_clustering() — toutes les res d'un coup
#     3. compute_umap()              — 1 fois par combo
#     4. plot_clustree()             — 1 fois par combo (toutes res visibles)
#     Pour chaque résolution :
#       5. compute_silhouette()      — métriques math
#       6. compute_ari_bootstrap()   — métriques math
#       7. plots UMAP, QC, batch, silhouette
#       8. Retourne une ligne de résultats
#   Fin boucle combos
#   9. Scoring composite + sélection des top candidats
#   10. Sur top-K candidats : compute_marker_quality() + compute_module_scores()
#       + compute_adt_rna_coherence() → bio_score
#   11. Classement final + export CSV
# ============================================================

# ── Paramètres du grid search ────────────────────────────────────────────────
GRID_PARAMS <- list(
  k_values   = c(20, 30, 40),
  res_values = seq(0.2, 1.2, by = 0.2),
  dims_configs = list(
    "13dims" = list(dims_list = list(1:13, 1:13), dims_use = 1:13),
    "20dims" = list(dims_list = list(1:20, 1:20), dims_use = 1:20)
  )
)

N_BOOTSTRAP   <- 5      # ARI bootstrap iterations (augmenter à 10-20 en prod)
N_TOP_COMBOS  <- 5      # Nombre de combos sur lesquels calculer les métriques bio
SIL_N_SAMPLE  <- 3000   # Cellules max pour silhouette approx.


# ── run_one_combo : traite une paire (k, dims_config) ────────────────────────
#
# Retourne un data.frame avec une ligne par résolution testée.
# Note : une seule passe WNN + clustering multi-résolution pour tout le combo.
# ─────────────────────────────────────────────────────────────────────────────
run_one_combo <- function(k, dims_config_name, seu, grid_params, base_dir,
                          n_bootstrap = N_BOOTSTRAP) {
  dc       <- grid_params$dims_configs[[dims_config_name]]
  dims_list <- dc$dims_list
  dims_use  <- dc$dims_use
  res_values <- grid_params$res_values

  message(sprintf("\n══ Combo k=%d | %s ══", k, dims_config_name))

  # ── Étape 1-3 : WNN + clustering multi-res + UMAP (1 seule fois) ───────────
  obj <- prepare_clustering_object(seu, dims_list, dims_use, k, res_values)
  seu_multi <- obj$seu

  # ── Étape 4 : Clustree (APRÈS avoir toutes les résolutions dans l'objet) ───
  k_dir <- make_dir(base_dir, sprintf("k%d_%s", k, dims_config_name))
  clustree_dir <- make_dir(k_dir, "Clustree")
  plot_clustree(seu_multi, k = k, out_dir = clustree_dir)

  # ── Étape 5-8 : Métriques + visualisations par résolution ──────────────────
  emb <- Embeddings(seu_multi, "harmony")[, dims_use]

  results_list <- lapply(res_values, function(res) {
    message(sprintf("  [res=%.2f]", res))

    res_col <- paste0("wsnn_res.", res)
    if (!res_col %in% colnames(seu_multi@meta.data)) {
      warning(sprintf("  Colonne '%s' manquante — résolution ignorée.", res_col))
      return(NULL)
    }

    # Fixer les identités sur la résolution courante
    Idents(seu_multi) <- res_col
    seu_multi$seurat_clusters <- Idents(seu_multi)  # pour les fonctions qui lisent seurat_clusters

    clust   <- as.integer(Idents(seu_multi))
    n_clust <- length(unique(clust))

    # Silhouette
    sil_res <- compute_silhouette(emb, clust, n_sample = SIL_N_SAMPLE)

    # ARI bootstrap
    ari_res <- compute_ari_bootstrap(
      seu        = seu,            # objet ORIGINAL (sans graph recalculé)
      dims_list  = dims_list,
      k          = k,
      res        = res,            # ← scalaire (bug corrigé)
      ref_clusters = clust,
      n_bootstrap  = n_bootstrap
    )

    # Score mathématique composite
    math_score <- compute_math_score(
      mean_sil  = sil_res$mean_global,
      mean_ari  = ari_res$mean,
      prop_neg  = sil_res$prop_negative
    )

    title_str <- make_title(k, res, n_clust, sil_res$mean_global,
                            sil_res$prop_negative, ari_res$mean)
    fill_cols <- pick_colors(n_clust)

    # Dossier de sortie pour cette résolution
    res_dir <- make_dir(k_dir, sprintf("res_%.2f", res))

    # Visualisations
    plot_silhouette(sil_res, fill_cols, title_str, make_dir(res_dir, "Silhouette"))
    plot_umap(seu_multi, fill_cols, title_str, sil_res, make_dir(res_dir, "UMAP"))
    plot_cluster_sizes(seu_multi, fill_cols, make_dir(res_dir, "Cluster_sizes"))
    plot_batch_repartition(seu_multi, fill_cols, make_dir(res_dir, "Batch"))
    plot_qc_metrics(seu_multi, fill_cols, make_dir(res_dir, "QC"))

    message(sprintf("  [k=%d | res=%.2f | %s] sil=%.3f | ARI=%.3f | score=%.3f",
                    k, res, dims_config_name,
                    sil_res$mean_global, ari_res$mean, math_score))

    data.frame(
      k              = k,
      dims_config    = dims_config_name,
      resolution     = res,
      n_clusters     = n_clust,
      mean_sil       = sil_res$mean_global,
      prop_negative  = sil_res$prop_negative,
      mean_ari       = ari_res$mean,
      sd_ari         = ari_res$sd,
      math_score     = math_score,
      bio_score      = NA_real_,   # rempli lors de l'étape bio
      final_score    = NA_real_,
      stringsAsFactors = FALSE
    )
  })

  dplyr::bind_rows(Filter(Negate(is.null), results_list))
}


# ── Étape biologique sur les top candidats ────────────────────────────────────
#
# Appelée UNE SEULE FOIS après le grid search complet.
# On recalcule le WNN uniquement pour les top combos.
# ─────────────────────────────────────────────────────────────────────────────
run_bio_evaluation <- function(top_combos, seu_path, grid_params, base_dir) {
  message("\n══ Évaluation biologique des top candidats ══")
  

  # Force le mode séquentiel — FindAllMarkers + FindMultiModalNeighbors
  # ont leur propre parallélisation interne qui entre en conflit avec future
  plan(sequential)
  
  results_bio <- lapply(seq_len(nrow(top_combos)), function(i) {
    row        <- top_combos[i, ]
    dc         <- grid_params$dims_configs[[row$dims_config]]
    dims_list  <- dc$dims_list
    dims_use   <- dc$dims_use
    k          <- row$k
    res        <- row$resolution

    message(sprintf("  [Bio] k=%d | %s | res=%.2f", k, row$dims_config, res))

    seu <- qread(seu_path)  # Charger une fois pour tous les combos bio
    
    # Recalcul WNN + clustering pour ce combo précis
    seu_tmp <- prepare_clustering_object(seu, dims_list, dims_use, k, res)$seu
    Idents(seu_tmp) <- paste0("wsnn_res.", res)
    seu_tmp$seurat_clusters <- Idents(seu_tmp)
    
    # Libère seu immédiatement — on n'en a plus besoin
    rm(seu); gc()

    # Module scores
    mod_res <- compute_module_scores(seu_tmp, BIO_SIGNATURES)
    seu_tmp <- mod_res$seu

    # ADT/RNA cohérence
    adt_res <- compute_adt_rna_coherence(seu_tmp, ADT_RNA_DICT)

    # Marker quality (coûteux — uniquement sur top candidats)
    DefaultAssay(seu_tmp) <- "SCT"
    mk_res <- compute_marker_quality(seu_tmp)

    bio_score <- compute_bio_score(
      marker_score  = mk_res$score,
      module_score  = mod_res$score,
      adt_rna_score = adt_res$score
    )

    # Plots additionnels bio
    bio_dir <- make_dir(base_dir,
                        sprintf("k%d_%s", k, row$dims_config),
                        sprintf("res_%.2f", res),
                        "Bio")
    if (!is.null(mk_res$markers)) {
      plot_volcano_per_cluster(mk_res$markers, make_dir(bio_dir, "Volcano"))
      plot_marker_heatmap(seu_tmp, mk_res$markers, out_dir = make_dir(bio_dir, "Heatmap"))
    }
    plot_module_score_umaps(seu_tmp, mod_res$score_col_names, make_dir(bio_dir, "ModuleScores"))
    
    # Libère seu immédiatement — on n'en a plus besoin
    rm(seu); gc()

    data.frame(
      k             = k,
      dims_config   = row$dims_config,
      resolution    = res,
      bio_score     = bio_score,
      marker_score  = mk_res$score,
      module_score  = mod_res$score,
      adt_rna_score = adt_res$score
    )
  })
  
  # Remet la parallélisation pour d'éventuels appels suivants
  plan(multisession, workers = 2)

  dplyr::bind_rows(results_bio)
}


# ── Pipeline principal ────────────────────────────────────────────────────────
run_grid_search <- function(seu_path, grid_params = GRID_PARAMS, base_dir,
                            n_bootstrap = N_BOOTSTRAP,
                            n_top_combos = N_TOP_COMBOS,
                            w_math = 0.5, w_bio = 0.5) {

  dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Augmente la limite pour les gros objets Seurat
  options(future.globals.maxSize = 2000 * 1024^2)  # 2 GB par worker
  plan(multisession, workers = 2)                   # 2 workers max pour éviter les problèmes de mémoire

  # ── Phase 1 : Grid search mathématique ──────────────────────────────────────
  # Parallélisation sur les combos (k × dims_config) — pas sur les résolutions
  # car FindMultiModalNeighbors n'est pas thread-safe en mémoire partagée.
  combos <- expand.grid(
    k               = grid_params$k_values,
    dims_config     = names(grid_params$dims_configs),
    stringsAsFactors = FALSE
  )

  all_results <- furrr::future_map_dfr(
    seq_len(nrow(combos)),
    function(i) {
      seu_local <- qread(seu_path)   # chaque worker charge son propre exemplaire
      run_one_combo(
        k                = combos$k[i],
        dims_config_name = combos$dims_config[i],
        seu              = seu_local,
        grid_params      = grid_params,
        base_dir         = base_dir,
        n_bootstrap      = n_bootstrap
      )
    },
    .options = furrr::furrr_options(seed = TRUE)
  )

  # ── Phase 2 : Scoring composite + sélection top candidats ───────────────────
  all_results <- all_results |>
    dplyr::arrange(dplyr::desc(math_score))

  top_combos <- head(all_results, n_top_combos)

  # ── Phase 3 : Évaluation biologique sur top candidats ───────────────────────
  # NB : séquentielle (FindAllMarkers ne se parallélise pas proprement)
  bio_results <- run_bio_evaluation(top_combos, seu_path, grid_params, base_dir)

  # ── Phase 4 : Score final combiné ───────────────────────────────────────────
  all_results <- all_results |>
    dplyr::left_join(bio_results |> dplyr::select(k, dims_config, resolution,
                                                   bio_score, marker_score,
                                                   module_score, adt_rna_score),
                     by = c("k", "dims_config", "resolution")) |>
    dplyr::mutate(
      bio_score   = dplyr::coalesce(bio_score.y, bio_score.x),
      final_score = dplyr::case_when(
        !is.na(bio_score) ~ w_math * math_score + w_bio * bio_score,
        TRUE              ~ math_score
      )
    ) |>
    dplyr::select(-dplyr::ends_with(".x"), -dplyr::ends_with(".y")) |>
    dplyr::arrange(dplyr::desc(final_score))

  # ── Phase 5 : Outputs globaux ───────────────────────────────────────────────
  plot_global_heatmaps(all_results, base_dir)

  write.csv(all_results,
            file.path(base_dir, "grid_search_results.csv"),
            row.names = FALSE)

  message("\n══ Grid search terminé ══")
  message(sprintf("  Meilleur combo : k=%d | %s | res=%.2f | final_score=%.4f",
                  all_results$k[1], all_results$dims_config[1],
                  all_results$resolution[1], all_results$final_score[1]))

  all_results
}


# ── Point d'entrée ─────────────────────────────────────────────────────────
# Décommenter pour lancer :
#
# source("00_config.R")
# # seu <- readRDS("path/to/your/seurat_object.rds")
#
# results <- run_grid_search(
#   seu         = seu,
#   grid_params = GRID_PARAMS,
#   base_dir    = file.path(DIR_RESULT, "Community_detection", "grid_search"),
#   n_bootstrap = 5,
#   n_top_combos = 5
# )
#
# # Meilleur combo :
# best <- results[1, ]
# print(best)
