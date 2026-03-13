############
# PACKAGES #
############

# BiocManager::install("DropletUtils")
library(DropletUtils)
# BiocManager::install("SingleCellExperiment")
library(SingleCellExperiment)
library(scater)
library(scran)
library(ggplot2)
library(qs)
library(dplyr)
library(Matrix)
library(patchwork)
# BiocManager::install("SingleR")
library(SingleR)
# BiocManager::install("celldex")
# install.packages("celldex")
library(celldex)

#########
# PATHS #
#########

setwd("~/Documents/Project/Thomas/20260116_scRNAseq_HD_MM_reanalysis/")

QCed_seurat_path <- "/mnt/crct-share/crct13/PPA/PPA_work/project/01_HCvsMM/01A_HCvsMM_10x/results/seurat_objects/QCed_seurat.qs"
sce_obj_path <- "02_processed_data/HD_vs_MM_scRNAseq_reanalysis_SCE.qs"
markers_path <- "03_result/HD_vs_MM_scRNAseq_reanalysis_markers_Louvain_res03.qs"

#################################################################################
#################################################################################
#################################################################################

# 1. Load data

# Read Seurat object
QCed_seurat <- qread(QCed_seurat_path)
sce_obj <- qread(sce_obj_path)

# As I want to reperform analysis and I don't really know what was done before, I will start from the raw count matrix (and metadata)
# Extract infos from Seurat object
# counts <- QCed_seurat@assays$RNA@counts
# meta <- QCed_seurat@meta.data
# 
# writeMM(counts, 
#         "01_raw_data/counts.mtx")
# 
# write.table(rownames(counts),
#             file = "01_raw_data/genes.tsv",
#             quote = FALSE,
#             row.names = FALSE,
#             col.names = FALSE)
# 
# write.table(colnames(counts),
#             file = "01_raw_data/barcodes.tsv",
#             quote = FALSE,
#             row.names = FALSE,
#             col.names = FALSE)
# 
# write.csv(meta,
#           "01_raw_data/cell_metadata.csv",
#           row.names = TRUE)

# Read matrix file
counts <- readMM("01_raw_data/counts.mtx")
# Read gene names
genes <- read.table("01_raw_data/genes.tsv",
                    header = FALSE,
                    stringsAsFactors = FALSE)
# Read cell barcodes
barcodes <- read.table("01_raw_data/barcodes.tsv",
                       header = FALSE,
                       stringsAsFactors = FALSE)
# Read cell metadata
meta <- read.csv("01_raw_data/cell_metadata.csv",
                 row.names = 1,
                 stringsAsFactors = FALSE)
# Create SingleCellExperiment object
sce_obj <- SingleCellExperiment(
  assays = list(counts = counts),
  colData = meta)
# Assign gene names to rowData
rowData(sce_obj)$Symbol <- genes$V1
rownames(sce_obj) <- genes$V1
# Assign cell barcodes to colnames
colnames(sce_obj) <- barcodes$V1
# Assign metadata to colData
colData(sce_obj) <- DataFrame(meta)

# Verify the SCE object
sce_obj

# Delete some colData columns that are not needed / or don't understand what they are
colnames(colData(sce_obj))
cols_to_remove <- c("nCount_ADT", "nFeature_ADT", "RNA_snn_res.0.8", 
                    "seurat_clusters", "ident", "scDblFinder.cluster",
                    "scDblFinder.class", "scDblFinder.score", "scDblFinder.difficulty",
                    "scDblFinder.cxds_score", "scDblFinder.mostLikelyOrigin", "scDblFinder.originAmbiguous",
                    "wsnn_res.0.5", "nCount_SCT", "nFeature_SCT", "SCT.weight", "ADT.weight")
colData(sce_obj) <- colData(sce_obj)[ , !colnames(colData(sce_obj)) %in% cols_to_remove ]
# Add a new column in metadata that combine condition column and source column
sce_obj$sample <- paste0(sce_obj$condition, "_", sce_obj$source)

meta <- colData(sce_obj)

sce_obj












# 2. Initial Exploration

# Dimensions: genes x cells
dim(sce_obj)

# Preview of rowData content (gene-level metadata)
head(rowData(sce_obj))

# Preview of colData content (cell-level metadata)
head(colData(sce_obj))

# Available assays (expression matrices)
assayNames(sce_obj)

# Access to the raw count matrix
counts <- assay(sce_obj, "counts")

# Preview of the count matrix (subset)
counts[1:5, 1:5]



# 3. Quality Control

## 3.1. Pre-QC

# Mitochondrial genes identification
mito_genes <- grepl("^MT-", rowData(sce_obj)$Symbol)
# Calculate QC metrics (with mitochondrial genes)
sce_obj <- addPerCellQCMetrics(
  sce_obj,
  subsets = list(Mito = mito_genes))
# View QC metrics
head(colData(sce_obj))

## 3.2. Visualization of pre-QC metrics

# Histogram of total counts per cell
total_plot <-
  hist(sce_obj$sum, breaks = 50,
       main = "Total Counts/UMI per Cell (Library Size)",
       xlab = "Total Counts",
       ylab = "Number of Cells")
# Histogram of detected genes per cell
detected_plot <-
  hist(sce_obj$detected, breaks = 50,
       main = "Detected Genes per Cell",
       xlab = "Number of Detected Genes",
       ylab = "Number of Cells")
# Boxplot of mitochondrial gene percentage
mito_plot <-
  boxplot(sce_obj$subsets_Mito_percent,
          main = "Mitochondrial Gene Percentage",
          ylab = "Percentage (%)")
# Combine plots
par(mfrow = c(3, 1))

# According 10X, tresholds for QC filtering are often:
# - Library size between 1000-20000 UMIs
# - Detected genes between 500-6000 genes
# - Mitochondrial content < 20%

# QC seems ok




# 4. Normalization
# Compute size factors using deconvolution method
sce_obj <- computeSumFactors(sce_obj)
# Normalize counts
sce_obj <- logNormCounts(sce_obj)
# New assay "logcounts" added to the SCE object





# 5. Feature Selection
# Identify highly variable genes (take time to run)
# dec <- modelGeneVar(sce_obj,
#                     assay.type = "logcounts",
#                     block = sce_obj$sample,
#                     span = 0.3)
# saveRDS(dec, file = "modelGeneVar_blockSample_span0.3.rds")
dec <- readRDS("modelGeneVar_blockSample_span0.3.rds")

# View top highly variable genes
head(dec[order(dec$bio, decreasing = TRUE), ])
# Select top 2000 highly variable genes
top_hvgs <- getTopHVGs(dec,
                       n = 2000)








# 5. Downstream Analysis

# 5.1. PCA
# PCA on highly variable genes
set.seed(123)  # For reproducibility
sce_obj <- runPCA(sce_obj,
                  subset_row = top_hvgs,
                  ncomponents = 30
                  )
# Visualize PCA
p1 <- plotPCA(sce_obj, colour_by = "sample")
p2 <- plotPCA(sce_obj, colour_by = "id")
p3 <- plotPCA(sce_obj, colour_by = "source")
p4 <- plotPCA(sce_obj, colour_by = "condition")

(p1 | p2) / 
(p3 | p4)


# 5.2. UMAP
# UMAP on PCA results
set.seed(123)  # For reproducibility
sce_obj <- runUMAP(sce_obj,
                   dimred = "PCA",
                   ncomponents = 2,
                   n_neighbors = 30,
                   min_dist = 0.3)
# Visualize UMAP
p5 <- plotUMAP(sce_obj, colour_by = "sample")
p6 <- plotUMAP(sce_obj, colour_by = "id")
p7 <- plotUMAP(sce_obj, colour_by = "source")
p8 <- plotUMAP(sce_obj, colour_by = "condition")

(p5 | p6) / 
(p7 | p8)



# 5.3 Neigboor graph
snn <- buildSNNGraph(sce_obj,
                     use.dimred = "PCA")



# 5.4. Louvain Clustering
for (res in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) {
  clusters <- igraph::cluster_louvain(snn,
                                      resolution = res)
  cluster_colname <- paste0("cluster_res", gsub("\\.", "", as.character(res)))
  sce_obj[[cluster_colname]] <- factor(clusters$membership)
  print(paste0("Clusters computed for resolution = ", res))
}

# clusters <- igraph::cluster_louvain(snn,
#                                     resolution = 0.1)
# sce_obj$cluster_res01 <- factor(clusters$membership)



# Visualize UMAP with clusters
for (res in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) {
  cluster_colname <- paste0("cluster_res", gsub("\\.", "", as.character(res)))
  p <- plotUMAP(sce_obj, colour_by = cluster_colname) +
    ggtitle(paste0("Louvain Clustering (res = ", res, ")"))
  # download plot
  ggsave(filename = paste0("03_result/Clusters_UMAP/UMAP_Louvain_Clustering_res_", res, ".png"),
         plot = p,
         width = 6,
         height = 5)
  print(paste0("UMAP plot saved for resolution = ", res))
}
# res = 0.3 seems to be ok for start, if we need to be more specific we can change it later


# (p9 <- plotUMAP(sce_obj, colour_by = "cluster"))
# (p10 <- plotUMAP(sce_obj, colour_by = "sample"))
# (p9 | p10)



# 5.5. Clusters markers
markers <- findMarkers(sce_obj,
                       groups = sce_obj$cluster_res03)

markers <- qread(markers_path)


# 5.6. Annotation of clusters
# View markers for cluster 1
markers_cluster1 <- markers[["1"]]
head(markers_cluster1)



# Load reference dataset
ref <- HumanPrimaryCellAtlasData()

# Automatic annotation using SingleR
pred <- SingleR(
  test = sce_obj,
  ref = ref,
  labels = ref$label.main,
  # clusters = sce_obj$cluster_res03
)

sce_obj$auto_label <- pred$labels[sce_obj$cluster_res03]
# Visualize UMAP with automatic labels
p11 <- plotUMAP(sce_obj, colour_by = "auto_label") +
  ggtitle("Automatic Cell Type Annotation (SingleR)")
p11

# Verify if my predictions macth with previous annotation in col predicted.celltype.l3
# cell_types <- c(
#   "gdT_2", "CD4 Naive", "CD4 TCM_1", "NK_2",
#   "CD8 TEM_5", "CD8 TEM_4", "CD4 TEM_1", "CD8 TEM_2",
#   "CD8 TEM_1", "CD4 TCM_3", "NK_1", "CD4 CTL",
#   "CD4 TEM_2", "Treg Memory", "gdT_4", "CD8 Naive",
#   "NK_CD56bright", "CD4 TEM_3", "CD8 Naive_2", "Treg Naive",
#   "NK_4", "CD8 TEM_3", "CD8 TCM_1", "gdT_1",
#   "CD8 TEM_6", "NK_3", "CD4 TCM_2", "MAIT",
#   "CD8 TCM_3", "CD8 TCM_2", "dnT_2", "NK Proliferating",
#   "gdT_3", "CD4 Proliferating", "dnT_1", "CD8 Proliferating"
# )
# 
# correspondance <- data.frame(
#   fine_cell_type = cell_types,
#   broad_cell_type = dplyr::case_when(
#     grepl("^NK", cell_types)            ~ "NK_cells",
#     grepl("CD4|CD8|Treg|gdT|dnT|MAIT", cell_types) ~ "T_cells",
#     TRUE                                ~ NA_character_
#   ),
#   stringsAsFactors = FALSE
# )
# 
# correspondance

# See hom many of my previous annotations match with automatic annotations
table(sce_obj$predicted.celltype.l3,
      sce_obj$auto_label)

# 6. Save object & Environmental Session Info

# Save sce_obj
qsave(sce_obj,
      "02_processed_data/HD_vs_MM_scRNAseq_reanalysis_SCE.qs")

# Save markers
qsave(markers,
      "03_result/HD_vs_MM_scRNAseq_reanalysis_markers_Louvain_res03.qs")

# Save predictions
qsave(pred,
      "03_result/HD_vs_MM_scRNAseq_reanalysis_predictions_cell_type_cluster_res03.qs")

# Session info
sessioninfo::session_info()






