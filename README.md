# SNH_AIH_spatial_transcriptomics
Workflows and code for spatial and meta-transcriptomic analysis of seronegative hepatitis and autoimmune hepatitis samples.


# Spatial Transcriptomics Figure Generation

This repository provides an **easy-to-use pipeline for pre-processing CosMx data and generating figures**. It is designed for liver disease research (e.g., autoimmune hepatitis (AIH), seronegative (SN) liver disease, and donor (D) samples) and can produce UMAPs, barplots, boxplots, heatmaps, dotplots, and CellChat-based cell-cell interaction visualizations.  

The scripts are structured for **reproducibility**, with figures saved automatically to a dedicated folder. This repository is ideal for immunology researchers looking to visualize cell types, gene expression patterns, and intercellular communication in spatial transcriptomics data.

---

## 📂 Repository Structure

The repository files are structured by the below:

```
├── .gitignore                 # 
├── 2.1.combine_all.r          # Script to combine all samples
├── 2.2.cell_label.r           # Script to label cell types
├── 2.3.singleR.r              # Script to use singleR to identify cell types
├── 2.4.immune.r               # Script to split cells into immune and non-immune to help with cell labelling
├── 2.5.polygon.r              # Script to 
├── 2.6.de.r                   # Script to do differentials
├── combined_plots.r           # Script to generate all figures
├── LICENSE                    # License information for this repository
└── README.md/                 # This document
```

---

## Data structure 

```
├── data/                          
│   ├── seurat_obj_processed.rds   # Preprocessed Seurat object
├── scripts/                       # Script for workflow
│   └── 2.1.combine_all.r          # Script to combine all samples
│   └── 2.2.cell_label.r           # Script to label cell types
│   └── 2.3.singleR.r              # Script to use singleR to identify cell types
│   └── 2.4.immune.r               # Script to split cells into immune and non-immune to help with cell labelling
│   └── 2.5.polygon.r              # Script to 
│   └── 2.6.de.r                   # Script to do differentials
│   └── combined_plots.r           # Script to generate all figures
└── figures/                       # Output folder for all generated plots
```

- **data/**: Store your preprocessed Seurat or CosMx object here. Must contain a `final_cell_types` column in `meta.data` for annotated cell types.
- **figures/**: All output figures (UMAP, barplot, boxplot, heatmap, dotplot, CellChat networks) will be automatically saved here.
- **scripts/**: Contains the combined plotting script. Designed to run as a single R script for all figure types.

---

## ⚙️ Requirements

- **R ≥ 4.2**
- Required packages:

```r
install.packages(c("dplyr", "ggplot2", "pheatmap", "patchwork", "ggrepel"))
BiocManager::install(c("SingleR", "CellChat", "Seurat"))
```

•	Input data: a Seurat/CosMx object with:
	•	Annotated cell types in meta.data$final_cell_types
	•	Expression layers including data and counts if using CosMx data

---

🚀 Usage
	1.	Clone the repository:

```
git clone https://github.com/AmberBozward/SNH_AIH_spatial_transcriptomics.git
cd spatial-figures
```

	2.	Place your preprocessed Seurat/CosMx object in data/:

```r
data/seurat_obj_processed.rds
```

	3.	Run the combined plotting script in R:

```r
source("scripts/combined_plots.R")
```

	4.	Check the figures/ folder for all generated outputs:
	•	UMAPs
	•	Barplots
	•	Boxplots
	•	Heatmaps
	•	Dotplots
	•	CellChat network plots

---

📝 Customization
	•	Modify genes for plots: Edit the genes_to_plot vector in the script for boxplots and dotplots.
	•	Change cell types: Ensure final_cell_types contains the desired annotations for coloring and subsetting.

---

📜 License

This repository is provided for academic and research purposes only. Please cite appropriately if used in publications.
