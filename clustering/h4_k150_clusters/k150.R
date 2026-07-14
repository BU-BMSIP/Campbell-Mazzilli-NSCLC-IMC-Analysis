#!/usr/bin/env Rscript
# Plain R script version for qsub batch testing
# Tazein Shah, 07/13/2026

library(imcRtools)
library(scran)
library(bluster)
library(BiocParallel)

cat("Starting k = 150 clustering script\n")
flush.console()

## 1. Load harmonized data
spe_path <- "/restricted/projectnb/camplab/projects/20220504_Suzuki_LN/TumorIMC/analysis/batch_correction/spe_combined_harmony_r4.rds"

if (!file.exists(spe_path)) {
  stop("File not found at: ", spe_path, "\nCheck that /restricted is mounted on this compute node.")
}

spe <- readRDS(spe_path)
cat("Loaded SPE with", ncol(spe), "cells\n")
flush.console()

## 2. Harmony embedding
harmony_mat <- reducedDim(spe, "Harmony")
cat("Harmony embedding dimensions:", dim(harmony_mat), "\n")
flush.console()

## 3. k value
k <- 150

## 4. Output path
final_path <- "/restricted/projectnb/camplab/projects/20220504_Suzuki_LN/TumorIMC/analysis/clustering/h4_k150_clusters/harmony_k150_clusters.rds"

## 5. Build graph and cluster
cat("\n--- k =", k, "---\n")
flush.console()

cat("Building SNN graph...\n")
flush.console()
g <- buildSNNGraph(t(harmony_mat), k = k, type = "rank",
                   BPPARAM = MulticoreParam(workers = 2))

cat("Running Louvain community detection...\n")
flush.console()
clust <- igraph::cluster_louvain(g)

spe$broad_clusters <- factor(clust$membership)
n_clusters <- nlevels(spe$broad_clusters)

cat("k =", k, "produced", n_clusters, "clusters\n")
flush.console()

## 6. Save final rds
saveRDS(spe, final_path)
cat("Saved final rds to:", final_path, "\n")

cat("\nDone.\n")