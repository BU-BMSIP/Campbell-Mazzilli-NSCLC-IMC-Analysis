#!/usr/bin/env Rscript
# All Patients: Broad Clustering k Sweep, Harmony (no asinh)
# Plain R script version for qsub batch testing
# Tazein Shah, 06/29/2026

set.seed(42)

## 1. Load Libraries

library(imcRtools)
library(tidyverse)
library(dittoSeq)
library(viridis)
library(scater)
library(scran)
library(bluster)
library(patchwork)
library(RColorBrewer)
library(BiocParallel)

## 2. Load Harmonized Data

spe_path <- "/restricted/projectnb/camplab/projects/20220504_Suzuki_LN/TumorIMC/analysis/batch_correction/spe_combined_harmony_r3.rds"

spe <- readRDS(spe_path)

cat("Cells:", ncol(spe), "\n")
cat("Reduced dims available:", paste(reducedDimNames(spe), collapse = ", "), "\n")
cat("Assays available:", paste(assayNames(spe), collapse = ", "), "\n")

## 3. Marker Sets

core_markers <- c("SMA", "CD4", "CD8a", "CD19", "CD15",
                  "CD14", "CD86", "LAMP3", "CD163", "EpCAM")
core_markers <- core_markers[core_markers %in% rownames(spe)]

full_markers <- c(
  "SMA", "CD4", "CD8a", "CD19", "CD15", "CD14", "CD86", "LAMP3",
  "Ki67", "MKI67", "CD138", "PD1", "PDL1", "FoxP3", "CD68", "CD163",
  "TIM3", "CD141", "CD11c", "CD56", "CD117", "EpCAM", "CD38",
  "CD45RO", "LAG3", "CD1c", "GATA3", "TTF1orCD16"
)
full_markers <- full_markers[full_markers %in% rownames(spe)]

cat("Core markers:", paste(core_markers, collapse = ", "), "\n")
cat("Full panel markers found:", length(full_markers), "\n")

## 4. Preprocessing: 99th Percentile Cap + Zscore

cap_99 <- function(mat) {
  t(apply(mat, 1, function(x) {
    q99 <- quantile(x, 0.99)
    pmin(x, q99)
  }))
}

zscore_rows <- function(mat) {
  t(scale(t(mat)))
}

mat_0100 <- counts(spe) |> cap_99() |> zscore_rows()
assay(spe, "scaled_0100") <- mat_0100
cat("Assays in SPE:", paste(assayNames(spe), collapse = ", "), "\n")

## 5. UMAP (reuse existing)

if ("UMAP_harmony" %in% reducedDimNames(spe)) {
  cat("Reusing existing UMAP_harmony.\n")
} else {
  set.seed(220225)
  spe <- runUMAP(spe, dimred = "Harmony", name = "UMAP_harmony")
  cat("UMAP computed.\n")
}

## 6. Subsample for Testing
# Set use_subsample to TRUE to test k values quickly on a smaller set of
# cells before committing to a full run. Set to FALSE for the full dataset.

use_subsample <- TRUE
cells_per_roi_n <- 5000

if (use_subsample) {
  set.seed(220225)
  cells_per_roi <- split(seq_len(ncol(spe)), spe$sample_id)
  sampled_cells <- unlist(lapply(cells_per_roi, function(idx) {
    sample(idx, min(cells_per_roi_n, length(idx)))
  }))
  spe <- spe[, sampled_cells]
  cat("Subsampled to", ncol(spe), "cells (", cells_per_roi_n, "per ROI )\n")
} else {
  cat("Running on full dataset,", ncol(spe), "cells\n")
}

## 7. k Values to Test
# Edit this vector to test different k values without touching the rest of the script

k_values <- c(200)

## 8. Parallel Backend

n_workers <- 8
bp <- MulticoreParam(workers = n_workers)
cat("Using MulticoreParam with", n_workers, "workers\n")

## 9. Run k Sweep

harmony_mat <- reducedDim(spe, "Harmony")

results_path <- "/restricted/projectnb/camplab/projects/20220504_Suzuki_LN/TumorIMC/analysis/clustering/ksweep_cluster_counts.csv"

k_results <- list()

for (k_val in k_values) {
  col_name <- paste0("clusters_k", k_val)
  cat("\n=== k =", k_val, "===\n")
  flush.console()
  
  set.seed(220225)
  g <- buildSNNGraph(t(harmony_mat), k = k_val, type = "rank", BPPARAM = bp)
  clust <- igraph::cluster_louvain(g)
  
  spe[[col_name]] <- factor(clust$membership)
  n_clusters <- nlevels(spe[[col_name]])
  
  cat("k =", k_val, "produced", n_clusters, "clusters\n")
  print(table(spe[[col_name]]))
  flush.console()
  
  k_results[[col_name]] <- n_clusters
  
  # append result immediately so progress is saved even if the job gets
  # killed partway through a later k value
  row <- data.frame(k = k_val, n_clusters = n_clusters,
                    subsample = use_subsample, timestamp = Sys.time())
  write.table(row, results_path, sep = ",", row.names = FALSE,
              col.names = !file.exists(results_path), append = file.exists(results_path))
  cat("Logged result to:", results_path, "\n")
  flush.console()
}

## 10. Summary of Cluster Counts per k

cat("\n=== Summary ===\n")
summary_df <- data.frame(
  k = k_values,
  n_clusters = unlist(k_results)
)
print(summary_df)

## 11. Save Intermediate
# Saved with all k columns so the winning k can be picked later without rerunning the sweep
# If use_subsample is TRUE, this intermediate only contains the subsampled
# cells, meant for quick k testing, not for final annotation.

intermediate_path <- if (use_subsample) {
  "/restricted/projectnb/camplab/projects/20220504_Suzuki_LN/TumorIMC/analysis/clustering/harmony_0100_ksweep_subsample_intermediate.rds"
} else {
  "/restricted/projectnb/camplab/projects/20220504_Suzuki_LN/TumorIMC/analysis/clustering/harmony_0100_ksweep_intermediate.rds"
}

saveRDS(spe, intermediate_path)
cat("\nSaved intermediate to:", intermediate_path, "\n")

## Notes
# This script only runs the sweep and saves the intermediate RDS.
# Dot plots, findMarkers, and annotation stay in the Rmd for interactive
# review, since those steps depend on visually inspecting output and
# picking a winning k, which does not fit a noninteractive batch job.
# Once a winning k is found on the subsample, set use_subsample to FALSE
# and rerun on the full dataset to get the actual clustering used for
# annotation. Then open the Rmd, load the resulting intermediate RDS
# instead of rerunning buildSNNGraph, and continue from chunk 8 onward.
# Cluster counts per k are appended to ksweep_cluster_counts.csv as each
# k value finishes, so results are visible even if the job is killed
# partway through a later k value.

cat("\nDone.\n")