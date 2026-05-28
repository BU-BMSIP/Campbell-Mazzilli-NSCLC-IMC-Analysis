# BMSIP 2026 Internship
**Project:** Immune changes in tumor tissue and lymph nodes associated with 
aggressive non-small cell lung cancer (NSCLC)  
**PIs:** Dr. Josh Campbell & Dr. Sarah Mazzilli  
**Interns:** Vaidehi Gupta & Tazein Shah  
**Start Date:** May 11, 2026  

## Project Overview
Building on prior work from the Campbell lab characterizing 
immune dysregulation in regional lymph nodes of NSCLC patients (Xi et al., 2026), 
this project analyzes matched IMC data from primary tumor tissue and adjacent 
normal regions across 11 patients. The goal is to identify immune cell populations 
and spatial cellular niches within the tumor microenvironment, compare findings 
across disease stages (IA vs. IB–IIIA), and correlate tumor immune features with 
previously characterized lymph node findings.

## Timeline (subject to change)
**Week 1-2 — Data Exploration**
- Review pre-print biological background and IMC methodology
- Visualize raw data using MCDViewer (subset of ROIs)
- Generate histograms and density plots, check for artifacts

**Week 3-5 — Preprocessing & QC**
- Run steinbock pipeline on tumor samples (preprocessing → segmentation → feature extraction)
- Quality control report (segmentation, image-level, cell-level)
- Import data into R using imcRtools

**Week 6-7 — Cell Phenotyping**
- Batch correction if needed
- Cell phenotyping using Celda (split between Vaidehi & Tazein)
- Cell type composition analysis

**Week 8-10 — Spatial Analysis & Deliverables**
- Spatial niche analysis using imcRtools
- Correlation with lymph node findings (Xi et al., 2026)
- Write full project report
- Prepare poster for August 2026 Poster Presentation 

## Methods
Imaging Mass Cytometry (IMC) · steinbock · DeepCell/Mesmer · 
imcRtools · Celda · Single-cell analysis · Spatial niche identification · R/Python

## References
- Xi et al. (2026) *medRxiv* https://doi.org/10.64898/2026.01.12.25343268v1
- Windhager et al. (2023) *Nature Protocols* https://doi.org/10.1038/s41596-023-00881-0
- Greenwald et al. (2022) *Nature Biotechnology* https://doi.org/10.1038/s41587-021-01094-0
