# Amrute_FAP_CAD_BiTE_2025

**Author:** Junedh M. Amrute  
**Project:** Targeting Modulated Vascular Smooth Muscle Cells in Atherosclerosis via FAP-Directed Immunotherapy
**Repository:** [jamrute/Amrute_FAP_CAD_BiTE_2025](https://github.com/jamrute/Amrute_FAP_CAD_BiTE_2025)  
**Date:** 2025  

---

## 🧠 Overview
This repository contains R scripts, analyses, and pipelines for "Targeting Modulated Vascular Smooth Muscle Cells in Atherosclerosis via FAP-Directed Immunotherapy"

---

## ⚙️ Requirements
- **R version:** ≥ 4.3  
- **Operating system:** macOS or Linux  
- **Hardware:** ≥ 64 GB RAM recommended for single-cell datasets  

---

## 📦 Dependencies
Install all required packages using:

```r
install.packages(c(
  "Seurat", "SeuratDisk", "sctransform", "harmony", "ArchR",
  "dplyr", "tidyverse", "Matrix", "data.table", "ggplot2",
  "ggpubr", "patchwork", "RColorBrewer", "ggsci", "scales",
  "pheatmap", "viridis", "biomaRt", "org.Hs.eg.db", "org.Mm.eg.db",
  "clusterProfiler", "ReactomePA", "enrichplot", "DOSE", "nichenetr"
))
remotes::install_github("GreenleafLab/ArchR")
remotes::install_github("saeyslab/nichenetr")

## Usage

git clone https://github.com/jamrute/Amrute_FAP_CAD_BiTE_2025.git
cd Amrute_FAP_CAD_BiTE_2025

## Structure

- CITE-seq Analysis:
  - QC: mergeFiles.Rmd (merge cell ranger output of the CITE-seq raw data and QC), run_scrublet_human_coronary.ipynb (run scrublet to remove doublets)
  - Integrate/annotate: global_annotate.Rmd (integrate data and use WNN clustering and perform DE expression to annotate populations)
  - GWAS enrichment: gwas_enrichment.Rmd
  - SMC analysis: postProcess_SMC.Rmd (analysis of SMC cell states: Fig. 2 of manuscript), stroma_palantir.ipynb (pseudotime of SMC states)
  - Myeloid analysis:

- Mouse-Human Integration:
  - Athero_regression_mouse_mapping_E-MTAB-12019.Rmd
  - Elencar_mouse_to_human_athero.Rmd
  - mouse_mapping_Qutermous.Rmd

- Spatial:
  - Visium FFPE: carotid_NCVR_integration.Rmd, coronary_tangram.ipynb, visium.preprocess.tangram.py, visium_FFPE_spatial_seurat.Rmd
  - Xenium: xenium.Rmd, xenium_morphology.Rmd

- BiTE therapy scRNAseq
  - murine_athero_BiTE_global.Rmd
  - murine_athero_BiTE_TCR.Rmd

## Contact
Junedh Amrute, jamrute@wustl.edu
