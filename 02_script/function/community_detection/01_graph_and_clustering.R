# ============================================================
# 01_graph_and_clustering.R
# Étape 1 : Construction du graph WNN (coûteuse, 1 seule fois par combo k/dims)
# Étape 2 : Clustering multi-résolution (stocke wsnn_res.X pour clustree)
# Étape 3 : UMAP (1 seul calcul par combo k/dims)
# ============================================================


# ── 1. Construction du graph WNN ──────────────────────────────────────────────
#
# POURQUOI séparer cette étape ?
#   - FindMultiModalNeighbors est l'opération la plus coûteuse.
#   - Pour un combo (k, dims) fixe, on ne doit la calculer qu'UNE FOIS,
#     puis tester toutes les résolutions sur le même graph.
#   - Cette fonction retourne un objet Seurat enrichi avec les slots
#     "wnn.nn", "wnn.graph", "wsnn".
#
# @param seu        Objet Seurat avec les réductions "harmony" et "apca"
# @param dims_list  list(dims_harmony, dims_apca), ex. list(1:13, 1:13)
# @param k          Nombre de voisins WNN
# @return           seu avec graph WNN calculé
# ─────────────────────────────────────────────────────────────────────────────
build_wnn_graph <- function(seu, dims_list, k) {
  message(sprintf("  [WNN] k=%d | dims_harmony=%d | dims_apca=%d",
                  k, max(dims_list[[1]]), max(dims_list[[2]])))

  FindMultiModalNeighbors(
    seu,
    reduction.list = list("harmony", "apca"),
    dims.list      = dims_list,
    k.nn           = k,
    verbose        = FALSE
  )
}


# ── 2. Clustering multi-résolution (pour clustree) ───────────────────────────
#
# POURQUOI une boucle externe et non resolution = res_values ?
#   FindClusters avec un vecteur de résolutions ajoute bien les colonnes
#   wsnn_res.X dans @meta.data — MAIS l'objet retourné a Idents() fixé
#   à la DERNIÈRE résolution, ce qui peut induire en erreur.
#   La boucle ci-dessous est explicite et garantit que chaque colonne
#   wsnn_res.X est correctement nommée et présente avant clustree.
#
# @param seu_wnn     Objet avec graph WNN déjà calculé
# @param res_values  Vecteur de résolutions, ex. seq(0.2, 1.2, 0.2)
# @return            seu avec colonnes wsnn_res.0.2, wsnn_res.0.4 … dans @meta.data
# ─────────────────────────────────────────────────────────────────────────────
run_multiresolution_clustering <- function(seu_wnn, res_values) {
  # FindClusters accepte nativement un vecteur — il ajoute toutes les colonnes
  # wsnn_res.X en une seule passe, ce qui est plus efficace que la boucle.
  # On l'utilise donc directement, mais on vérifie ensuite la présence des colonnes.
  seu_wnn <- FindClusters(
    seu_wnn,
    resolution = res_values,
    graph.name = "wsnn",
    algorithm  = 4,      # Leiden — plus stable que Louvain pour scRNA-seq
    verbose    = FALSE
  )

  # Vérification : toutes les colonnes wsnn_res.X doivent exister
  expected_cols <- paste0("wsnn_res.", res_values)
  missing       <- setdiff(expected_cols, colnames(seu_wnn@meta.data))
  if (length(missing) > 0) {
    warning("Colonnes manquantes après FindClusters : ", paste(missing, collapse = ", "))
  }

  seu_wnn
}


# ── 3. UMAP (1 calcul par combo k/dims) ──────────────────────────────────────
#
# UMAP est calculé sur la réduction harmony (ou wnn.umap si disponible).
# Il doit être calculé UNE FOIS par combo et réutilisé pour toutes les
# résolutions — pas recalculé à chaque résolution.
#
# @param seu_clustered  Objet Seurat avec clustering multi-résolution
# @param dims_use       Dimensions harmony utilisées pour UMAP
# @param reduction      "harmony" ou "wnn" (pour RunUMAP basé sur WNN graph)
# @return               seu avec réduction "umap" calculée
# ─────────────────────────────────────────────────────────────────────────────
compute_umap <- function(seu_clustered, dims_use, reduction = "harmony") {
  if (reduction == "wnn") {
    # UMAP basé sur le graph WNN directement (recommandé pour données multimodales)
    RunUMAP(
      seu_clustered,
      nn.name = "weighted.nn",
      reduction.name = "umap",
      verbose = FALSE
    )
  } else {
    RunUMAP(
      seu_clustered,
      reduction      = reduction,
      dims           = dims_use,
      reduction.name = "umap",
      verbose        = FALSE
    )
  }
}


# ── 4. Wrapper : graph + clustering multi-résolution + UMAP ──────────────────
#
# Fonction principale appelée une fois par combo (k, dims).
# Retourne un objet Seurat prêt pour :
#   - clustree (toutes les colonnes wsnn_res.X présentes)
#   - le calcul de métriques par résolution (silhouette, ARI, bio)
#
# @param seu        Objet Seurat de référence (non modifié)
# @param dims_list  list(dims_harmony, dims_apca)
# @param dims_use   Dimensions harmony pour UMAP et silhouette
# @param k          Nombre de voisins
# @param res_values Vecteur de résolutions à tester
# @return           Liste : $seu (objet enrichi), $k, $dims_use, $res_values
# ─────────────────────────────────────────────────────────────────────────────
prepare_clustering_object <- function(seu, dims_list, dims_use, k, res_values) {
  message(sprintf("[k=%d | dims=1:%d] Construction WNN + clustering multi-résolution...",
                  k, max(dims_use)))

  seu_out <- seu |>
    build_wnn_graph(dims_list = dims_list, k = k) |>
    run_multiresolution_clustering(res_values = res_values) |>
    compute_umap(dims_use = dims_use, reduction = "wnn")

  list(
    seu       = seu_out,
    k         = k,
    dims_use  = dims_use,
    res_values = res_values
  )
}
