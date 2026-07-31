#!/usr/bin/env Rscript
# =============================================================================
# All Patients: Broad Clustering, Harmony (no asinh), k20
# Converted from 05_broad_clustering_harmony_0100_k20.Rmd for qsub / Rscript use
# Tazein Shah
# =============================================================================
# Clustering only: loads the harmonized SPE, runs Louvain on the Harmony
# embeddings, and saves the clustered SPE. No plotting, annotation, or
# findMarkers here, those stay in the Rmd for interactive review.

set.seed(42)

suppressPackageStartupMessages({
  library(imcRtools)
  library(SingleCellExperiment)
  library(scran)
  library(bluster)
})

# -----------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------
spe_path <- "/restricted/projectnb/camplab/projects/20220504_Suzuki_LN/TumorIMC/analysis/batch_correction/spe_combined_harmony_r5.rds"

cluster_path <- "/restricted/projectnb/camplab/projects/20220504_Suzuki_LN/TumorIMC/analysis/final_cell_annotation/h5_cluster.rds"

# -----------------------------------------------------------------------
# Load Harmonized Data
# -----------------------------------------------------------------------
spe <- readRDS(spe_path)

cat("Cells:", ncol(spe), "\n")
cat("Reduced dims available:", paste(reducedDimNames(spe), collapse = ", "), "\n")

# -----------------------------------------------------------------------
# Louvain Clustering on Harmony Embeddings
# -----------------------------------------------------------------------
if (file.exists(cluster_path)) {
  cat("Clustering file already exists, skipping clustering:", cluster_path, "\n")
} else {
  set.seed(220225)
  
  harmony_mat <- reducedDim(spe, "Harmony")
  cat("Harmony embedding dimensions:", dim(harmony_mat), "\n")
  
  g     <- buildSNNGraph(t(harmony_mat), k = 20, type = "rank")
  clust <- igraph::cluster_louvain(g)
  spe$broad_clusters <- factor(clust$membership)
  
  cat("Number of broad clusters:", nlevels(spe$broad_clusters), "\n")
  print(table(spe$broad_clusters))
  
  dir.create(dirname(cluster_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(spe, cluster_path)
  cat("Saved clustered intermediate SPE to:", cluster_path, "\n")
}

cat("\nDone.\n")