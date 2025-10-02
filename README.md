# SNH_AIH_spatial_transcriptomics

# Visium Pre-processing the Data and Figure Generation

This repository provides an **easy-to-use pipeline for pre-processing Visium data and generating figures**. It is designed for liver disease research (e.g., autoimmune hepatitis (AIH), seronegative (SN) liver disease, and donor (D) samples) and can produce UMAPs, barplots, boxplots, heatmaps and dotplots.

---

## 📂 Repository Structure

The repository files are structured by the below:

```
├── .gitignore                 
├── 1.QC/r                                      # Script for inital QC
├── 2a.integration.light.r                      # Script for light integration
├── 2b.integration.strict.r                     # Script for strict integration
├── 3a.cell_types_public.r                      # Script to label cell types
├── 3b.spatial_deconvolution_strict.r           # Script to strictly deconvolute cell types
├── 3c.spatial_deconvolution_light.r            # Script to lightly deconvolute cell types
├── 3d.pseudobulk_linneages_deconvolution.r     # Script to pseudobulk deconvoluted cell types
├── 3f_spatial_deconvolution_lineagges_light.r  # Script to pseudobulk deconvoluted cell types
├── 4a_region_barcodes.r                        # Script to generate barcodes from regions of interest (parenchyma vs non-parenchyma)
├── 4b_region_barcodes_pseudo.r                 # Script to pseudobulk chosen regions
├── 5a.plots.r                                  # Script to generate all figures
├── LICENSE                                     # License information for this repository
└── README.md/                                  # This document
```

---

## ⚙️ Requirements

- **R ≥ 4.2**
- Required packages:

```r
install.packages(c("dplyr", "ggplot2", "pheatmap", "patchwork", "ggrepel"))
BiocManager::install(c("SingleR", "CellChat", "Seurat"))
```

•	Input data: a Seurat object with:
	•	Annotated cell types in meta.data$final_cell_types
	•	Expression layers including data and counts

---

## 📜 License

This repository is provided for academic and research purposes only. Please cite appropriately if used in publications.
