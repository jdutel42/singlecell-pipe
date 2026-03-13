library(qs)
library(Seurat)
library(ggplot2)
library(patchwork)
library(clustree)
library(harmony)
library(patchwork)




seu <- qread("/mnt/crct-share/crct13/PPA/PPA_work/project/01_HCvsMM/01A_HCvsMM_10x/results/Tcells_analysis/seurat_objects/QCed_seurat_Tcells_clean.qs")

# RNA
seu <- NormalizeData(seu, assay = "RNA")
seu <- FindVariableFeatures(seu, assay = "RNA")
seu <- ScaleData(seu, assay = "RNA")
seu <- RunPCA(seu)

ElbowPlot(seu, ndims = 50)

#batch correction with Harmony (RNA)
seu <- RunHarmony(seu, 
                  group.by.vars = "id", 
                  plot_convergence = TRUE)


# ADT
seu <- NormalizeData(seu, normalization.method = "CLR", margin = 2, assay = "ADT")
DefaultAssay(seu) <- "ADT"
VariableFeatures(seu) <- rownames(seu[["ADT"]])
seu <- ScaleData(seu, assay = "ADT")
seu <- RunPCA(seu, reduction.name = "apca")


ElbowPlot(seu, ndims = 30, reduction = "apca")

# Neihbors WNN (commun RNA + ADT)
seu <- FindMultiModalNeighbors(
  seu,
  reduction.list = list("harmony", "apca"),
  dims.list = list(1:16, 1:13)
)

# Clustering (commun RNA + ADT on WNN graph)
# Explore different resolutions
resolutions <- seq(0.01, 0.1, by = 0.01)
for (res in resolutions) {
  seu <- FindClusters(seu, graph.name = "wsnn", resolution = res, verbose = FALSE)
}
clustree(seu, prefix = "wsnn_res.")
# Choose resolution 0.06
seu$seurat_clusters <- seu$wsnn_res.0.06
# reorder cluster numbers
seu$seurat_clusters <- factor(seu$seurat_clusters, levels = c(
  "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"
))

# UMAP
seu <- RunUMAP(seu,
               nn.name = "weighted.nn",
               min.dist = 0.1,
               spread = 1.0,
               n.neighbors = 50,
               seed.use = 42)

# Visualization
# Dossier pour sauvegarder
output_dir <- "../03_result/UMAP_metadata_harmony"
dir.create(output_dir, showWarnings = FALSE)
metadata <- c("condition", "source", "id", "predicted.celltype.l1", "predicted.celltype.l2", "predicted.celltype.l3", "seurat_clusters")
for (meta in metadata) {
  p <- DimPlot(seu, reduction = "umap", group.by = meta) + ggtitle(meta)
  # Nom de fichier
  file_name <- paste0(output_dir, "/UMAP_metadata_", meta, ".png")
  
  ggsave(file_name, plot = p, width = 6, height = 5)
}


############

# Paramètres à tester
min.dists <- c(0.04, 0.06, 0.08, 0.1, 0.12, 0.14)
spreads <- c(0.5, 1.0, 1.5)
n.neighbors <- c(60)
# Dossier pour sauvegarder
output_dir <- "../03_result/UMAP_tests"
dir.create(output_dir, showWarnings = FALSE)

for (md in min.dists) {
  for (sp in spreads) {
    for (nn in n.neighbors) {
      seu <- RunUMAP(seu,
                     nn.name = "weighted.nn",
                     min.dist = md,
                     spread = sp,
                     n.neighbors = 60,
                     seed.use = 42,
                     reduction.name = paste0("umap_md", md, "_sp", sp, "_nn", nn))
      
      p <- DimPlot(seu,
                   reduction = paste0("umap_md", md, "_sp", sp, "_nn", nn),
                   group.by = "seurat_clusters") +
        ggtitle(paste("min.dist =", md, "spread =", sp, "n.neighbors =", nn))
      
      # Nom de fichier
      file_name <- paste0(output_dir, "/UMAP_md", md, "_sp", sp, "_nn", nn, ".png")
      
      ggsave(file_name, plot = p, width = 6, height = 5)
    }
  }
}
seu$active.umap <- seu$umap_md0.1_sp1.5_nn60
#Reorder cluster numbers numerically to 15 by numbers of cells inside
seu$seurat_clusters <- factor(seu$seurat_clusters, levels = names(sort(table(seu$seurat_clusters), decreasing = TRUE)))
seu$seurat_clusters <- factor(seu$seurat_clusters, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"))


############

# Visualize expression of ADT vs RNA for CD4
# DefaultAssay(seu) <- "ADT"
# rownames(seu[["ADT"]])
# p1 <- FeaturePlot(seu, "CD4-TotalSeqB", cols = c("lightgrey", "darkgreen")) + ggtitle("CD4 protein")
# DefaultAssay(seu) <- "RNA"
# grep("CD8", rownames(seu[["RNA"]]), value = TRUE)
# p2 <- FeaturePlot(seu, "CD4") + ggtitle("CD4 RNA")
# 
# # place plots side-by-side
# p1 | p2

# Visualize expression of multiple ADT markers
FeaturePlot(seu,
            features = c(
              "CD4-TotalSeqB",
              "CD8a-TotalSeqB",
              "CD45RA-TotalSeqB",
              "CD197-CCR7-TotalSeqB",
              "CD62L-TotalSeqB",
              "CD25-TotalSeqB",
              "TIGIT-VSTM3-TotalSeqB",
              "CD226-DNAM-1-TotalSeqB"),
            ncol = 3)

DefaultAssay(seu) <- "RNA"
FeaturePlot(seu, features = c("CD4", "CD8A", "FOXP3"), ncol = 3, reduction = "umap")
DefaultAssay(seu) <- "ADT"
FeaturePlot(seu, features = c("CD4-TotalSeqB", "CD8a-TotalSeqB"), ncol = 2, reduction = "umap")
DimPlot(seu, group.by = "condition", reduction = "umap")
DimPlot(seu, group.by = "predicted.celltype.l2", reduction = "umap")

dict <- list(
  "CD16-TotalSeqB" = "FCGR3A",
  "CD335-NKp46-TotalSeqB" = "NCR1",
  "CD197-CCR7-TotalSeqB" = "CCR7",
  "CD4-TotalSeqB" = "CD4",
  "CD8a-TotalSeqB" = "CD8A",
  "CD45RA-TotalSeqB" = "PTPRC",
  "CD62L-TotalSeqB" = "SELL",
  "CD127-IL-7Ra-TotalSeqB" = "IL7R",
  "CD25-TotalSeqB" = "IL2RA",
  "TIGIT-VSTM3-TotalSeqB" = "TIGIT",
  "CD226-DNAM-1-TotalSeqB" = "CD226",
  "CD38-TotalSeqB" = "CD38"
)


dir.create("../03_result/ADT_features_harmony", showWarnings = FALSE)

for (prot in names(dict)) {
  gene <- dict[[prot]]
  
  # ADT plot
  DefaultAssay(seu) <- "ADT"
  p1 <- FeaturePlot(seu, features = prot, reduction = "umap") + ggtitle(paste0(prot, " (ADT)"))
  
  # RNA plot
  DefaultAssay(seu) <- "RNA"
  # Vérifie que le gène existe dans RNA
  if (gene %in% rownames(seu[["RNA"]])) {
    p2 <- FeaturePlot(seu, features = gene, reduction = "umap") + ggtitle(paste0(gene, " (RNA)"))
  } else {
    p2 <- ggplot() + 
      ggtitle(paste0(gene, " (RNA) - Not found")) + 
      theme_void()
  }
  
  # Combine plots côte à côte
  p <- p1 | p2
  
  # Sauvegarde
  ggsave(filename = paste0("../03_result/ADT_features_harmony/", prot, "_featureplot.png"),
         plot = p, width = 10, height = 4)
}



# Cluster comparison

## Marqueur diff RNA
DefaultAssay(seu) <- "RNA"
seu <- FindVariableFeatures(
  seu,
  selection.method = "vst",
  nfeatures = 3000
)
seu <- ScaleData(
  seu,
  features = VariableFeatures(seu)
)

markers_rna <- FindAllMarkers(
  seu,
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.25
)

dir.create("../03_result/cluster_markers", showWarnings = FALSE)
write.csv(markers_rna, file = "../03_result/cluster_markers/cluster_markers_RNA.csv")

features_use <- unique(top10$gene)

seu <- ScaleData(
  seu,
  features = features_use
)

DoHeatmap(
  seu,
  features = features_use,
  group.by = "seurat_clusters"
)


## Marqueur diff ADT
DefaultAssay(seu) <- "ADT"
markers_adt <- FindAllMarkers(
  seu,
  only.pos = TRUE,
  test.use = "wilcox"
)
write.csv(markers_adt, file = "../03_result/cluster_markers/cluster_markers_ADT.csv")

#  DotPlot ADT markers
DotPlot(seu, 
        features = rownames(seu[["ADT"]]),
        group.by = "seurat_clusters") + RotatedAxis()

# VioloinPlot ADT markers
VlnPlot(
  seu,
  features = c("CD4-TotalSeqB", "CD8a-TotalSeqB"),
  assay = "ADT",
  group.by = "seurat_clusters",
  pt.size = 0
)
VlnPlot(
  seu,
  features = c("CD38-TotalSeqB", "CD226-DNAM-1-TotalSeqB", "TIGIT-VSTM3-TotalSeqB"),
  assay = "ADT",
  group.by = "seurat_clusters",
  pt.size = 0
)

# Compare clusters per condition
prop <- prop.table(table(seu$seurat_clusters, seu$condition), margin = 2)
prop

library(reshape2)
df <- melt(prop)

ggplot(df, aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(stat = "identity") +
  ylab("Condition") +
  xlab("Proportion") +
  ggtitle("Proportion of clusters per condition")


    Idents(seu) <- "condition"
  markers_c2 <- FindMarkers(
    seu,
    ident.1 = "MM",
    ident.2 = "HC",
    subset.ident = "2"
  )
  

