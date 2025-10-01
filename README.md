# SNH_AIH_spatial_transcriptomics

This repository contains the analysis code and workflows associated with the manuscript:

“Spatial transcriptomics links hepatocyte-macrophage interactions to viral signatures in seronegative hepatitis”
Submitted to Nature Communications, 2025

The repository includes workflows for CosMx single-cell spatial transcriptomics, Visium whole-transcriptome profiling, and metatranscriptomic analysis of explanted liver tissue, along with scripts for data processing, integration, and figure generation.

⸻

**Repository Structure**

SNH_spatial_transcriptomics/

	•	CosMx/ – Single-cell spatial transcriptomics workflows
    •	preprocessing/ – Quality control, segmentation, filtering
	  •	analysis/ – Clustering, cell-cell interactions, pathway analysis
	  •	visualization/ – UMAPs, heatmaps, spatial plots
	•	Visium/ – Whole-transcriptome spatial workflows
	  •	preprocessing/ – QC, alignment, gene filtering
	  •	integration/ – Integration with CosMx and histology
	  •	visualization/ – Spatial maps, differential expression plots
	•	Metatranscriptomics/ – Microbial and viral transcriptomics workflows
	  •	preprocessing/ – Host read removal, QC
	  •	analysis/ – HERV/viral gene quantification
	  •	visualization/ – Heatmaps, volcano plots
	•	Supplementary/ – Miscellaneous scripts, helper functions
	•	environment.yml – Conda environment for reproducibility
	•	LICENSE – MIT License
	•	.gitignore – Ignored files
	•	README.md – This document
  
⸻

**Contributors**

Dr Amber Bozward — CosMx and Visium workflows, integration, figure generation

Dr Mahboobeh Behruznia - Metatranscriptomics sequencing analysis

John Cole - CosMx and Visium workflows, integration, figure generation

Chiranjit Das - Figure generation

Kylie Savoye - CosMx neighbourhood analysis, figure generation

Professor Ye Oo -
  
⸻

**License**

This repository is licensed under the MIT License (see LICENSE file).

⸻

**Citation**

If you use this code, please cite:

Bozward et al., “Spatial transcriptomics links hepatocyte-macrophage interactions to viral signatures in seronegative hepatitis”, Nature Communications (2025).
