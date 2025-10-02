library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(devtools)
library(hdf5r)
library(glmGamPoi)
library(SingleR)
library(celldex)
library(presto)
library(EnhancedVolcano)
library(presto)
library(harmony)
library(SingleR)
library(celldex)


#### Parameters ####


# set memory available
options(future.globals.maxSize = 100000 * 1024^2)


#### Load data ####


# folders
setwd("C:/Users/JohnJ/Desktop/analysis/Ye_Oo/Amber_B_CosMx_2024/")
de_folder = "analysis/combined_all/singleR/"

# load sc object
sc.data = readRDS("analysis/combined_all/all_integrated_named.RDS")



#### SingleR ####


# Download avaialable human primary cell atlas markers
hpca.se = celldex::HumanPrimaryCellAtlasData()

# Extract expression information
exp.ma = sc.data@assays$RNA$data

# Get predictions
predictions_broad = SingleR(test = exp.ma, ref = hpca.se,labels = hpca.se$label.main)
predictions_fine = SingleR(test = exp.ma, ref = hpca.se,labels = hpca.se$label.fine)

# add predictions to metadata
sc.data = AddMetaData(sc.data, predictions_broad$labels, col.name = "predictions_broad")


# UMAP
p1 = DimPlot(sc.data, reduction = "integrated_umap", group.by = "sample_name", pt.size = 1)
p2 = DimPlot(sc.data, reduction = "integrated_umap", group.by = "sample_group", pt.size = 1)
p3 = DimPlot(sc.data, reduction = "integrated_umap", group.by = "predictions_broad", pt.size = 1, label = TRUE)
p4 = DimPlot(sc.data, reduction = "integrated_umap", group.by = "named_clusters_cell_type", pt.size = 1, label = TRUE)

ggp = p1 + p2 + p3 + p4
png(paste0(de_folder, "/RNA_named_UMAP.png"), height = 4000, width = 4000)
print(ggp)
dev.off()



#### Save RDS ####

saveRDS(sc.data, file="analysis/combined_all/all_integrated_named.RDS")

























