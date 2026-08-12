# BMSIP 2026: Immune changes in tumor tissue and lymph nodes associated with aggressive non-small cell lung cancer (NSCLC)

**PIs:** Dr. Josh Campbell & Dr. Sarah Mazzilli
**Interns:** Vaidehi Gupta & Tazein Shah

## Project Overview

Building on prior work from the Campbell lab characterizing immune
dysregulation in regional lymph nodes of NSCLC patients (Xi et al.,
2026), this project analyzes matched IMC data from primary tumor tissue
and adjacent normal regions across 11 patients. The goal is to identify
immune cell populations and spatial cellular niches within the tumor
microenvironment, compare findings across disease stages (Benign/IA2/IA3
vs. IB–IIB), and correlate tumor immune features with previously
characterized lymph node findings.

## Analysis Folder Organization

### batch_correction

5 batch correction runs of the tumor IMC data, plus overcorrection
checks. Latest (run 5):

- **Batch correction (run 5)**: combines all 12 steinbock outputs
  (P10–P19 + both Patient20 acquisitions, 5/19 added back this round),
  tags `disease_state`, resolves a `sample_id` collision risk between
  the two Patient20 acquisitions, drops Patient17's faulty ROI 1/ROI
  4, filters by area (20–500 px), then runs PCA/UMAP before and after
  Harmony (batch var = `sample_id`). Saves UMAPs to
  `figures/batch_correction_umaps/run5/` and
  `spe_combined_harmony_r5.rds`.
- **Overcorrection check (run 5)**: three diagnostics —
  - marker expression before/after Harmony (CD19, CD38, CD8a,
    FoxP3 + clean controls), including a spillover check on CD19
    against same-element isotope neighbors
  - per-cluster patient enrichment heatmaps (>2x with >50 cells
    flagged), with a Patient20 acquisition comparison and a
    Patient17/Cluster 19 deep dive
  - LISI (batch-LISI should rise, cluster-LISI hold/improve) on a
    20% patient-stratified subsample

### celltype_prop_analysis

sccomp cell type composition analysis on the finalized (r3) annotation,
plus poster figures:

- **sccomp (Disease State)**: `~ DiseaseState + (1 | Patient)`, two
  contrasts (Lower & Higher vs. Benign; Higher vs. Lower). Forest
  plots per contrast (FDR < 0.05, CI excludes 0) and a stacked bar
  chart of per-sample proportions by disease state. Saves
  `sccomp_celltype_vs_benign.rds`,
  `sccomp_celltype_higher_v_lower.rds`.
- **Poster figures (Higher vs. Lower)**: reuses that sccomp result for
  a curated diverging bar chart (+ a "full" version), and small-
  multiples boxplots by patient. Includes a display-only rename (Tumor
  Epithelial TTF1++ → "Lung Epithelial (TTF1+)", per PI feedback).

### clustering

k-value sweeps (k = 20, 150, 200) for picking the clustering resolution,
an uncorrected-data baseline, and cluster diagnostics for the B cell,
unknown, T cell, and epithelial populations.

- **k-sweep batch scripts** (qsub): Harmony embedding → SNN graph →
  Louvain, with an optional per-ROI subsample for quick testing.
  Outputs land in per-k subfolders (`k20_clusters/`, `k150_clusters/`,
  `k200_clusters/`, `h4_k150_clusters/`).
- **Uncorrected baseline clustering**: Louvain on `PCA_uncorrected`
  (34 patient-driven clusters) vs. Harmony-corrected clusters — shows
  how much Harmony reorganized the data.
- **Unknown cluster diagnostics**: is k150's "Unknown" cluster real or
  a resolution artifact? Cross-checks k20/k150/k200, plus Ki67,
  composition, and area checks on the stable core.
- **B cell cluster comparison**: cross-resolution check on B
  cell/plasma clusters — confirms two stable subtypes
  (CD38/FoxP3-leaning vs. CD138/CD15-leaning); k20 is the coarser
  outlier.
- **Epithelial cluster diagnostics**: cluster 10 likely mislabeled
  ("Transitional Myeloid" → actually epithelial, 93–94% land in Tumor
  Epithelial at k150/k200); cluster 18 small/patient-enriched;
  clusters 13/16 screened as possible duplicates.
- **T cell cluster diagnostics**: cluster 11 ("Treg") is actually
  B/plasma-lineage, not T cell; cluster 20 ("Exhausted T Cell") too
  small/weak to annotate confidently; lighter checks on clusters 8 and
  12.

### final_cell_annotation

Final cell annotations, moved here after `clustering`. Contains one
final broad cluster run, one subclustering folder (bcells, cluster15,
dc_mregdc, epithelial, myeloid, tcells), and three hierarchical
annotation runs (`r1`–`r3`), each with an `annotation_inspection`
subfolder of marker violin plots.

- **Broad clustering (Harmony V5)**: 20 broad clusters from
  `h5_cluster.rds` — preprocessing, UMAPs, violin/dot plots,
  `findMarkers()`, and ambiguous-cluster diagnostics (e.g. cluster 17
  flagged as low-quality/poorly-stained, not a real population). Ends
  with the manual annotation lookup table; saves
  `h5_cluster_annot.rds`.
- **Violin plots by cell type category**: groups the ~38 fine-grained
  r3 labels into 5 broad categories, then plots per-marker violins
  within each.
- **Cell type distribution across patients/ROIs**: composition bar
  charts by patient and ROI, a broad-vs-fine nesting consistency
  check, and flags for annotations dominated by one patient (>25%) or
  ROI (>7%).
- **T cell subcluster relabeling (r2 → r3)**: recodes 7 T cell label
  strings per PI feedback (naive/memory/exhaustion naming), no cell
  IDs changed. Saves `r3_annotation/spe_final_celltype_r3.rds`.

### figures

Batch correction UMAPs and cell segmentation images for all patients.

### images

Outlier distribution ridgeplots for all markers.

### patient20_QC

QC notebook (`Patient20_QC.Rmd`) comparing Patient20's two acquisitions
(4/25, 5 ROIs vs. 5/19, 6 ROIs): raw/filtered pixel distributions,
spatial marker maps, segmentation QC (cell counts, area filtering, SNR),
mask overlays, and expression heatmaps. **Conclusion:** both pass QC and
are kept; 5/19 chosen as the reference for cohort-wide batch effect
assessment.

### QC_scripts

Notebooks comparing raw / filtered / post-segmentation marker
intensities:

- **Cross-patient comparison**: ridgeplots across P12/14/17/18/19/20
  for all four processing stages.
- **Per-patient pixel-level QC**: histogram grid, sparsity summary
  table, raw-vs-filtered density plots (parameterized by patient).
- **Per-patient cell-level distributions**: same idea at the
  post-segmentation level, matched against the raw pixel files.

### segmentation_scripts

Cell segmentation mask visualizations for the full cohort:

- **`cell_segmentation.Rmd`**: per-patient report with marker overlays
  matched to Xi et al. (2026) cell types (SMA, CD4, CD8a, CD19, CD15,
  CD14+CD86, CD14+CD163, LAMP3).
- **`render_cell_seg.R`**: batch-renders that report across all
  patients to `cell_seg_reports/`.
- **Four Marker Overview**: DNA191/EpCAM/CD68/SMA overlays for every
  patient, saved to `figures/cellseg/`.

### spatial_analysis

Spatial/niche analysis of the annotated tumor IMC data. Three runs
(`run_1`–`run_3`); `run_3` is current, with three notebooks that hand
off in sequence:

1. **`spe_spatial_graph_r3.Rmd`**: builds the cell-cell spatial graph
   (`type = "expansion"`, `threshold = 20`) from the r2 annotation.
   Saves `spe_spatial_graph_r3.rds`.
2. **Niche generation + differential abundance**: k-means into 50
   niches, builds Stage/N-status/disease-state covariates, fits
   `sccomp` for both Stage and N-status, includes patient-driven-
   cluster risk checks. Saves to `niche_k50_thresh20/`; Stage is the
   primary model.
3. **`niche_annotation_composition_summary.Rmd`**: composition
   summaries for Stage-significant niches (dot plot, heatmap,
   patient-dominance check, top-cell-type tables, bar charts), and
   auto-drafted `#Niche N: ...` comments to
   `niche_annotation_drafts.txt`.

### visualizations

Cytomapper/spatial overlays of annotated cell types, plus poster
figures:

- **Patient18 ROI 1 overlays**: 13-category display grouping,
  all-cells + per-cell-type highlight plots.
- **Cytomapper cell type overlays**: Patient10 (Benign) vs. Patient12
  (Higher-stage) — mask boundaries (`plotCells()`) and centroid points
  (`plotSpatial()`) for a curated cell type shortlist.
- Poster figures: patient metadata summary
  (`data_information_summary.Rmd`), niche composition, and other
  summary plots.

## References

Xi et al. (2026). *medRxiv*.
<https://doi.org/10.64898/2026.01.12.25343268v1>

Windhager et al. (2023). *Nature Protocols*.
<https://doi.org/10.1038/s41596-023-00881-0>