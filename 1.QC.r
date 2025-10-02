#### Load Libraries ####

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(devtools)
library(hdf5r)

#### Guides and Manuals ####

# https://satijalab.org/seurat/reference/read10x_image
# https://satijalab.org/seurat/reference/read10x_image
# https://satijalab.org/seurat/reference/load10x_spatial
# https://yu-tong-wang.github.io/talk/sc_st_data_analysis_R.html


#### https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_07_spatial.html
#### https://satijalab.org/seurat/articles/spatial_vignette.html
#### https://satijalab.org/seurat/articles/integration_large_datasets.html

#### Function to load data ####

slice_counter = 0 
load_visium = function(image_path, image_filename, data_path, data_filename, slice_name)
  
{
  
  # loads the low res image
  sample_image = Read10X_Image(image_path, image.name=image_filename)
  
  # loads the spatial data
  sample_data = Load10X_Spatial(data_path, filename=data_filename, image = sample_image, slice = paste0("slice",slice_counter))
  slice_counter = slice_counter + 1
  
  # adds identity
  sample_data = SetIdent(sample_data, value = slice_name)
  sample_data[["orig.ident"]]$orig.ident = slice_name
  
  # get percent mitochondria
  sample_data = PercentageFeatureSet(sample_data, "^MT-", col.name = "percent_mito")
  sample_data = PercentageFeatureSet(sample_data, "HBB", col.name = "percent_hb")
  
  return(sample_data)
  
}



#### Load data ####

## WD
setwd("C:/Users/JohnJ/Desktop/analysis/Amber_Bozward_Visium_2023/data/")


## load data

data_AH_1607 = load_visium("autoimmune\\1607\\spatial", "tissue_hires_image.png", "autoimmune\\1607\\", "filtered_feature_bc_matrix.h5", "AH_1607")
data_AH_2709 = load_visium("autoimmune\\2709\\spatial", "tissue_hires_image.png", "autoimmune\\2709\\", "filtered_feature_bc_matrix.h5", "AH_2709")
data_AH_7538 = load_visium("autoimmune\\7538\\spatial", "tissue_hires_image.png", "autoimmune\\7538\\", "filtered_feature_bc_matrix.h5", "AH_7538")
data_AH_6079 = load_visium("autoimmune\\6079\\spatial", "tissue_hires_image.png", "autoimmune\\6079\\", "filtered_feature_bc_matrix.h5", "AH_6079")
data_AH_7014 = load_visium("autoimmune\\7014\\spatial", "tissue_hires_image.png", "autoimmune\\7014\\", "filtered_feature_bc_matrix.h5", "AH_7014")

data_SN_606 = load_visium("seronegative\\606\\spatial", "tissue_hires_image.png", "seronegative\\606\\", "filtered_feature_bc_matrix.h5", "SN_606")
data_SN_827 = load_visium("seronegative\\827\\spatial", "tissue_hires_image.png", "seronegative\\827\\", "filtered_feature_bc_matrix.h5", "SN_827")
data_SN_1042 = load_visium("seronegative\\1042\\spatial", "tissue_hires_image.png", "seronegative\\1042\\", "filtered_feature_bc_matrix.h5", "SN_1042")
data_SN_2739 = load_visium("seronegative\\2739\\spatial", "tissue_hires_image.png", "seronegative\\2739\\", "filtered_feature_bc_matrix.h5", "SN_2739")
data_SN_3303 = load_visium("seronegative\\3303\\spatial", "tissue_hires_image.png", "seronegative\\3303\\", "filtered_feature_bc_matrix.h5", "SN_3303")
data_SN_4393 = load_visium("seronegative\\4393\\spatial", "tissue_hires_image.png", "seronegative\\4393\\", "filtered_feature_bc_matrix.h5", "SN_4393")


## Upscale Images
data_AH_1607@images$slice0@scale.factors$lowres = data_AH_1607@images$slice0@scale.factors$hires
data_AH_2709@images$slice0@scale.factors$lowres = data_AH_2709@images$slice0@scale.factors$hires
data_AH_7538@images$slice0@scale.factors$lowres = data_AH_7538@images$slice0@scale.factors$hires
data_AH_6079@images$slice0@scale.factors$lowres = data_AH_6079@images$slice0@scale.factors$hires
data_AH_7014@images$slice0@scale.factors$lowres = data_AH_7014@images$slice0@scale.factors$hires

data_SN_606@images$slice0@scale.factors$lowres = data_SN_606@images$slice0@scale.factors$hires
data_SN_827@images$slice0@scale.factors$lowres = data_SN_827@images$slice0@scale.factors$hires
data_SN_1042@images$slice0@scale.factors$lowres = data_SN_1042@images$slice0@scale.factors$hires
data_SN_2739@images$slice0@scale.factors$lowres = data_SN_2739@images$slice0@scale.factors$hires
data_SN_3303@images$slice0@scale.factors$lowres = data_SN_3303@images$slice0@scale.factors$hires
data_SN_4393@images$slice0@scale.factors$lowres = data_SN_4393@images$slice0@scale.factors$hires






#### merge data ####

## merge
liver.merge = merge(data_AH_1607, c(data_AH_2709, data_AH_7538, data_AH_6079,data_AH_7014,data_SN_606,data_SN_827,data_SN_1042,data_SN_2739,data_SN_3303,data_SN_4393), merge.data = TRUE, add.cell.ids = c("AH_1607", "AH_2709", "AH_7538", "AH_6079","AH_7014","SN_606","SN_827","SN_1042","SN_2739","SN_3303","SN_4393"))


#### QC and Filter ####

## WD
setwd("C:/Users/JohnJ/Desktop/analysis/Amber_Bozward_Visium_2023/analysis/")

## QC plots
plot1 = VlnPlot(liver.merge, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito", "percent_hb"), pt.size = 0.1, ncol = 4) + NoLegend()
plot2 = SpatialFeaturePlot(liver.merge, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito","percent_hb"), pt.size.factor = 6, stroke = NA)

## Filter spots
liver.merge = liver.merge[, liver.merge$nFeature_Spatial > 500 & liver.merge$percent_mito < 25 & liver.merge$percent_hb < 20]

## QC plots
plot3 = SpatialFeaturePlot(liver.merge, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito","percent_hb"), pt.size.factor = 6, stroke = NA)


## Most expressed genes
C = liver.merge@assays$Spatial@counts
C@x = C@x/rep.int(colSums(C), diff(C@p))
most_expressed <- order(Matrix::rowSums(C), decreasing = T)[20:1]

png("QC/QC_top_genes.png", height = 500, width = 500)
print(boxplot(as.matrix(t(as.matrix(C[most_expressed, ]))), cex = 0.1, las = 1, xlab = "% total count per cell",col = (scales::hue_pal())(20)[20:1], horizontal = TRUE))
dev.off()

## remove contaminant genes
dim(liver.merge)

# Filter Mitocondrial
liver.merge = liver.merge[!grepl("^MT-", rownames(liver.merge)), ]

# Filter Heamoglobin genes
liver.merge = liver.merge[!grepl("HBB", rownames(liver.merge)), ]
liver.merge = liver.merge[!grepl("HBA1", rownames(liver.merge)), ]
liver.merge = liver.merge[!grepl("HBA2", rownames(liver.merge)), ]


dim(liver.merge)


## save QC plots

png("1a_QC/QC_violin_plots.png", height = 1500, width = 2500)
print(plot1)
dev.off()

png("1a_QC/QC_spatial_all_voxels.png", height = 5000, width = 5000)
print(plot2)
dev.off()

png("1a_QC/QC_spatial_filtered_voxels.png", height = 5000, width = 5000)
print(plot3)
dev.off()




#### Determine optimal components ####

## Jack Straw
#liver.merge = JackStraw(liver.merge, num.replicate = 100)
#liver.merge = ScoreJackStraw(liver.merge, dims = 1:20)
#JackStrawPlot(liver.merge, dims = 1:15)




#### Transform ####

## Transform
liver.merge = SCTransform(liver.merge, assay = "Spatial", verbose = TRUE, method = "poisson")

## QC plots
plot5 = VlnPlot(liver.merge, features = c("ACTB", "ALB", "SERPINA1", "FCN1", "LYZ", "S100A8"), pt.size = 0.1, ncol = 6) + NoLegend()
plot6 = SpatialFeaturePlot(liver.merge, features = c("ACTB", "ALB", "SERPINA1", "FCN1", "LYZ", "S100A8"), pt.size.factor = 6, stroke = NA)

png("1a_QC/QC_marker_genes.png", height = 2500, width = 2500)
print(plot5)
dev.off()

png("1a_QC/QC_spatial_marker_genes.png", height = 5000, width = 5000)
print(plot6)
dev.off()



#### Reduce and cluster ####


# RUn PCA
liver.merge = RunPCA(liver.merge, assay = "SCT", verbose = FALSE)
# compute K nearest neighbors (KNN)
liver.merge = FindNeighbors(liver.merge, reduction = "pca", dims = 1:30)
# Leiden algorithm for community detection
liver.merge = FindClusters(liver.merge, verbose = FALSE)
# PCA result is the default UMAP input, use dimensions 1:30 as input features
liver.merge = RunUMAP(liver.merge, reduction = "pca", dims = 1:30)


# Plot PCA heatmap
png("1a_QC/pca_heatmap.png", height = 1000, width = 1000)
DimHeatmap(liver.merge, dims = 1:15, cells = 2500, balanced = TRUE) + theme(text = element_text(size = 20))
dev.off()


# Plot UMAP
plot = DimPlot(liver.merge, reduction = "umap", group.by = c("ident", "orig.ident"))
png("1a_QC/cluster_umap.png", height = 1250, width = 2500)
print(plot)
dev.off()


# Plot spatial clusters
plot = SpatialDimPlot(liver.merge, pt.size.factor = 6, stroke = NA)
png("1a_QC/spatial_clusters.png", height = 1000, width = 7500)
print(plot)
dev.off()


# Plot clusters heatmap
liver.merge = PrepSCTFindMarkers(liver.merge)
liver.markers = FindAllMarkers(liver.merge, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
liver.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

png("1a_QC/clusters_heatmap.png", height = 1000, width = 1000)
DoHeatmap(liver.merge, features = top10$gene) + NoLegend() + theme(text = element_text(size = 20))
dev.off()



# save
saveRDS(liver.merge, file = "liver_final.rds")





