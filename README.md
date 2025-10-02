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

### R Version
- R >= 4.2.0

### Required R Packages
Make sure the following packages are installed before running the scripts:

- `Seurat` (>= 5.0)
- `ggplot2`
- `dplyr`
- `tidyr`
- `patchwork`
- `Matrix`
- `stringr`
- `cowplot`
- `ggrepel`
- `gridExtra`
- `hdf5r`
- `spatstat` (if using spatial statistics)
- `SeuratDisk` (if working with `.h5Seurat` or converting from AnnData)

You can install missing packages with:
```r
install.packages(c("ggplot2", "dplyr", "tidyr", "patchwork", 
                   "Matrix", "stringr", "cowplot", "ggrepel", "gridExtra"))
```

And from Bioconductor if needed:

```r
install.packages("Seurat")
install.packages("SeuratDisk")
install.packages("spatstat.geom")
install.packages("spatstat.core")
```

### Input Data

The scripts assume 10x Genomics Space Ranger outputs in the following format:

data/
  sample1/
    filtered_feature_bc_matrix.h5
    spatial/
      tissue_hires_image.png
      scalefactors_json.json
  sample2/
    ...

---

## 📜 License

This repository is provided for academic and research purposes only. Please cite appropriately if used in publications.
