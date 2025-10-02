#### Load Libraries ####

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(devtools)
library(hdf5r)

#### Guides and Manuals ####

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
  
  # get percent mitochondria and blood
  sample_data = PercentageFeatureSet(sample_data, "^MT-", col.name = "percent_mito")
  sample_data = PercentageFeatureSet(sample_data, "HBB", col.name = "percent_hb")
  
  ## Filter spots
  sample_data = sample_data[, sample_data$nFeature_Spatial > 500 & sample_data$percent_mito < 25 & sample_data$percent_hb < 20]

  # Should we filter albumin???
  
  # Filter Mitocondrial
  sample_data = sample_data[!grepl("^MT-", rownames(sample_data)), ]
  
  # Filter Heamoglobin genes
  sample_data = sample_data[!grepl("HBB", rownames(sample_data)), ]
  sample_data = sample_data[!grepl("HBA1", rownames(sample_data)), ]
  sample_data = sample_data[!grepl("HBA2", rownames(sample_data)), ]
  
  return(sample_data)
  
}



#### Load data ####

## WD
setwd("C:/Users/admin/Desktop/analysis/Ye_Oo/Amber_Bozward_visium_2023/liver_data/")


## load data

data_AH_2709 = load_visium("autoimmune\\2709\\spatial", "tissue_lowres_image.png", "autoimmune\\2709\\", "filtered_feature_bc_matrix.h5", "AH_2709")
data_AH_6079 = load_visium("autoimmune\\6079\\spatial", "tissue_lowres_image.png", "autoimmune\\6079\\", "filtered_feature_bc_matrix.h5", "AH_6079")
data_AH_7014 = load_visium("autoimmune\\7014\\spatial", "tissue_lowres_image.png", "autoimmune\\7014\\", "filtered_feature_bc_matrix.h5", "AH_7014")

data_SN_606 = load_visium("seroneg\\606\\spatial", "tissue_lowres_image.png", "seroneg\\606\\", "filtered_feature_bc_matrix.h5", "SN_606")
data_SN_827 = load_visium("seroneg\\827\\spatial", "tissue_lowres_image.png", "seroneg\\827\\", "filtered_feature_bc_matrix.h5", "SN_827")
data_SN_1042 = load_visium("seroneg\\1042\\spatial", "tissue_lowres_image.png", "seroneg\\1042\\", "filtered_feature_bc_matrix.h5", "SN_1042")
data_SN_2739 = load_visium("seroneg\\2739\\spatial", "tissue_lowres_image.png", "seroneg\\2739\\", "filtered_feature_bc_matrix.h5", "SN_2739")
data_SN_3303 = load_visium("seroneg\\3303\\spatial", "tissue_lowres_image.png", "seroneg\\3303\\", "filtered_feature_bc_matrix.h5", "SN_3303")
data_SN_4393 = load_visium("seroneg\\4393\\spatial", "tissue_lowres_image.png", "seroneg\\4393\\", "filtered_feature_bc_matrix.h5", "SN_4393")



#### integrate data ####

# create a list of the original data that we loaded to start with
st.list = list(AH_2709 = data_AH_2709, AH_6079 = data_AH_6079,AH_7014 = data_AH_7014,SN_606 = data_SN_606,SN_827 = data_SN_827,SN_1042 = data_SN_1042,SN_2739 = data_SN_2739,SN_3303 = data_SN_3303,SN_4393 = data_SN_4393)

# run SCT on all datasets
st.list = lapply(st.list, SCTransform, assay = "Spatial", method = "poisson")

# tidy up 
rm(data_AH_2709, data_AH_6079,data_AH_7014,data_SN_606,data_SN_827,data_SN_1042,data_SN_2739,data_SN_3303,data_SN_4393)
gc()


# need to set maxSize for PrepSCTIntegration to work
options(future.globals.maxSize = 4000 * 1024^2)  # set allowed size to 2K MiB
st.features = SelectIntegrationFeatures(st.list, nfeatures = 3000, verbose = FALSE)
st.list = PrepSCTIntegration(object.list = st.list, anchor.features = st.features,verbose = FALSE)

# get anchors  - data_AH_2709, data_SN_3303
int.anchors = FindIntegrationAnchors(object.list = st.list, normalization.method = "SCT",verbose = FALSE, anchor.features = st.features,  reference = c(1, 8))

# tidy up
rm(st.list)
gc()

# integrate
liver.integrated = IntegrateData(anchorset = int.anchors, normalization.method = "SCT",verbose = FALSE)


# tidy up 
rm(int.anchors)
gc()


#### Reduce and cluster Integrated Data ####

## WD
setwd("C:/Users/admin/Desktop/analysis/Ye_Oo/Amber_Bozward_visium_2023/analysis/")


liver.integrated = ScaleData(liver.integrated, verbose = FALSE)
liver.integrated = RunPCA(liver.integrated, verbose = FALSE)
liver.integrated = FindNeighbors(liver.integrated, reduction = "pca", dims = 1:30)
liver.integrated = FindClusters(liver.integrated, verbose = FALSE)
liver.integrated = RunUMAP(liver.integrated,  reduction = "pca", dims = 1:30)

# Plot PCA heatmap
png("DR/strict_integrated_pca_heatmap.png", height = 1500, width = 1500)
DimHeatmap(liver.integrated, dims = 1:15, cells = 2500, balanced = TRUE) + theme(text = element_text(size = 20))
dev.off()


# Plot UMAP
plot = DimPlot(liver.integrated, reduction = "umap", group.by = c("ident", "orig.ident"))
png("DR/strict_integrated_cluster_umap.png", height = 1250, width = 2500)
print(plot)
dev.off()


# Plot spatial clusters
plot = SpatialDimPlot(liver.integrated, image.alpha = 0, stroke = NA)
png("DR/strict_integrated_spatial_clusters.png", height = 1000, width = 7500)
print(plot)
dev.off()


# Plot clusters heatmap
liver.integrated = PrepSCTFindMarkers(liver.integrated)
liver.markers = FindAllMarkers(liver.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
liver.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

png("DR/strict_integrated_clusters_heatmap.png", height = 2000, width = 2000)
DoHeatmap(liver.integrated, features = top10$gene) + NoLegend() + theme(text = element_text(size = 20))
dev.off()




#### Save integrated table for analysis ####

saveRDS(liver.integrated, file = "liver_integrated_strict.rds")






#### Hepatocyte and Kuphfer markers ####

plot = SpatialFeaturePlot(liver.integrated, features = c("ALB","FGA", "TF", "C7", "MYL9", "CD74"), image.alpha = 0, stroke = NA)
png("spatial_features/strict_integrated_2Cluster_markers_spatial.png", height = 3000, width = 2500)
print(plot)
dev.off()


plot = SpatialFeaturePlot(liver.integrated, features = c("APOA2","CYP2A6","IGKC","CFTR","MYL9","CYP2E1","CXCL8","COL3A1","FTL","EPO","CCL21","MALAT1"), image.alpha = 0, stroke = NA)
png("spatial_features/strict_integrated_12Cluster_markers_spatial.png", height = 6000, width = 2500)
print(plot)
dev.off()





