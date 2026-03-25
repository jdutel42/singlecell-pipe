# ============================================================
# 00_config.R — Paramètres globaux et palettes
# ============================================================
# Charger ce fichier en premier dans tous les scripts.

library(Seurat)
library(furrr)
library(future)
library(cluster)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(patchwork)
library(clustree)
library(mclust)
library(ComplexHeatmap)
library(ggrepel)
library(ggridges)
library(scales)

# ── Parallélisation ──────────────────────────────────────────
# Leiden (igraph) n'est pas thread-safe dans les workers R.
# multisession (processus isolés) est le seul mode sûr avec Seurat.
plan(multisession, workers = min(parallel::detectCores() - 1, 6))

# ── Palettes ─────────────────────────────────────────────────
PAL_55 <- c(
  "#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00",
  "#A65628","#F781BF","#999999","#66C2A5","#FC8D62",
  "#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494",
  "#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E",
  "#E6AB02","#A6761D","#666666","#8DD3C7","#BEBADA",
  "#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5",
  "#D9D9D9","#BC80BD","#CCEBC5","#FFED6F","#1F78B4",
  "#33A02C","#FB9A99","#E31A1C","#FDBF6F","#CAB2D6",
  "#6A3D9A","#FFFF99","#B15928","#8E0152","#DE77AE",
  "#7FBC41","#4EB265","#4393C3","#2166AC","#B2182B",
  "#F4A582","#92C5DE","#FDDBC7","#D1E5F0","#542788"
)

# Choisit automatiquement assez de couleurs pour n clusters
pick_colors <- function(n) {
  if (n <= length(PAL_55)) PAL_55[seq_len(n)]
  else grDevices::colorRampPalette(PAL_55)(n)
}

# ── Signatures biologiques ────────────────────────────────────
# Adapter selon le type cellulaire ciblé.
BIO_SIGNATURES <- list(
  CD4_naive       = c("CCR7", "SELL", "TCF7", "LEF1", "IL7R"),
  CD4_Tcm         = c("CCR7", "IL7R", "CD27", "CD28"),
  CD4_Tem         = c("GZMK", "CCL5", "NKG7", "CXCR3"),
  CD4_Treg        = c("FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF18"),
  CD4_Tfh         = c("CXCR5", "BCL6", "PDCD1", "ICOS", "IL21"),
  CD4_Th17        = c("RORC", "IL17A", "IL17F", "CCR6", "RORA"),
  CD4_Th1         = c("IFNG", "TBX21", "CXCR3", "IL12RB2"),
  CD8_naive       = c("CCR7", "SELL", "TCF7", "LEF1"),
  CD8_effector    = c("GZMB", "PRF1", "IFNG", "NKG7", "GNLY"),
  CD8_exhausted   = c("PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4"),
  Cycling         = c("MKI67", "TOP2A", "STMN1", "PCNA", "CDK1")
)

# ── Dictionnaire ADT → RNA (CITE-seq) ────────────────────────
ADT_RNA_DICT <- list(
  "CD4-TotalSeqB"              = "CD4",
  "CD8a-TotalSeqB"             = "CD8A",
  "CD197-CCR7-TotalSeqB"       = "CCR7",
  "CD45RA-TotalSeqB"           = "PTPRC",
  "CD62L-TotalSeqB"            = "SELL",
  "CD127-IL-7Ra-TotalSeqB"     = "IL7R",
  "CD25-TotalSeqB"             = "IL2RA",
  "TIGIT-VSTM3-TotalSeqB"      = "TIGIT",
  "CD226-DNAM-1-TotalSeqB"     = "CD226",
  "CD38-TotalSeqB"             = "CD38",
  "CD16-TotalSeqB"             = "FCGR3A",
  "CD335-NKp46-TotalSeqB"      = "NCR1"
)
