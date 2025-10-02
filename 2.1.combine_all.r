

#### Load Libraries ####

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


#### Variables ####

# folders

wd = "C:/Users/JohnJ/Desktop/analysis/Ye_Oo/Amber_B_CosMx_2024/"
qc_folder = "analysis/combined_all/QC/"
markers_folder = "analysis/combined_all/markers/"


# set memory available
options(future.globals.maxSize = 1000000 * 1024^2)  # set allowed size


#### Load data ####

# set wd
setwd("C:/Users/JohnJ/Desktop/analysis/Ye_Oo/Amber_B_CosMx_2024/")


# RDS objects
sample_paths = c("C:/Users/JohnJ/Desktop/analysis/Ye_Oo/Amber_B_CosMx_2024/analysis/AIH_1607/filtered_named.RDS",
                 "C:/Users/JohnJ/Desktop/analysis/Ye_Oo/Amber_B_CosMx_2024/analysis/AIH_6079/filtered_named.RDS",
                 "C:/Users/JohnJ/Desktop/analysis/Ye_Oo/Amber_B_CosMx_2024/analysis/SN_606/filtered_named.RDS",
                 "C:/Users/JohnJ/Desktop/analysis/Ye_Oo/Amber_B_CosMx_2024/analysis/SN_1042/filtered_named.RDS", 
                 "C:/Users/JohnJ/Desktop/analysis/Ye_Oo/Amber_B_CosMx_2024/analysis/SN_2739/filtered_named.RDS",
                 "C:/Users/JohnJ/Desktop/analysis/Ye_Oo/Amber_B_CosMx_2024/analysis/D_7678/filtered_named.RDS",
                 "C:/Users/JohnJ/Desktop/analysis/Ye_Oo/Amber_B_CosMx_2024/analysis/D_6414/filtered_named.RDS",
                 "C:/Users/JohnJ/Desktop/analysis/Ye_Oo/Amber_B_CosMx_2024/analysis/D_6446/filtered_named.RDS")


# sample names
sample_names = c("AIH_1607","AIH_6079","SN_606","SN_1042","SN_2739","D_7678","D_6414","D_6446")

# groups
sample_groups = c("AIH","AIH","SN","SN","SN","D","D","D")

# empty list
sc.data = list()

# loop through samples and add to seurat list
for (index in 1:length(sample_paths))
{
  sample_path = sample_paths[index]
  sample_name = sample_names[index]
  sample_group = sample_groups[index]
  
  sc.data.sample = readRDS(sample_path)
  sc.data.sample = AddMetaData(sc.data.sample, metadat = sample_name, col.name = "sample_name")
  sc.data.sample = AddMetaData(sc.data.sample, metadat = sample_group, col.name = "sample_group")
  
  sc.data[[sample_name]] = sc.data.sample
}

# merge
sc.data = merge(x = sc.data[[1]], y = sc.data[-1])



#### QC data ####


# set to RNA
DefaultAssay(sc.data) = 'RNA'

# active ident
sc.data = SetIdent(sc.data, value = "orig.ident")

# QC violins
ggp = VlnPlot(sc.data, features = c("nCount_RNA", "nFeature_RNA"), pt.size = 0, ncol = 4) + NoLegend()
png(paste0(qc_folder, "/QC_RNA_violin.png"), height = 750, width = 1250)
print(ggp)
dev.off()

# QC scatter
ggp = FeatureScatter(sc.data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
png(paste0(qc_folder, "/QC_RNA_scatter.png"), height = 750, width = 1250)
print(ggp)
dev.off()


#### Variable Features ####


# Perform log-normalization
sc.data = NormalizeData(sc.data)
sc.data = ScaleData(sc.data)

# find variable features
sc.data = FindVariableFeatures(sc.data, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top50 = head(VariableFeatures(sc.data), 50)

# plot variable features with and without labels
p1 = VariableFeaturePlot(sc.data)
p2 = LabelPoints(plot = p1, points = top50, repel = TRUE)
ggp = p1 + p2
png(paste0(qc_folder, "/QC_RNA_variable_features.png"), height = 1000, width = 1800)
print(ggp)
dev.off()


#### Normalise data & PCA ####


# set to RNA
DefaultAssay(sc.data) = 'RNA'

# Run the PCA 
sc.data = RunPCA(object = sc.data, assay ="RNA", reduction.name='pca_unint')

# PCA Elbow plot
ggp = ElbowPlot(sc.data, reduction = "pca_unint")
png(paste0(qc_folder, "/QC_unfiltered_RNA_PCA_elbow.png"), height = 500, width = 500)
print(ggp)
dev.off()

# Dimensions to use
max_dims = 15


# PCA scatter plot
ggp = DimPlot(sc.data, reduction = "pca_unint")
png(paste0(qc_folder, "/QC_unfiltered_RNA_PCA_scatter.png"), height = 1000, width = 1000)
print(ggp)
dev.off()

# PCA loadings
ggp = VizDimLoadings(sc.data, dims = 1:max_dims, reduction = "pca_unint")
png(paste0(qc_folder, "/QC_unfiltered_RNA_PCA_loadings.png"), height = 1000, width = 1800)
print(ggp)
dev.off()

# PCA heatmap
png(paste0(qc_folder, "/QC_unfiltered_RNA_PCA_heatmaps.png"), height = 1800, width = 1000)
print(DimHeatmap(sc.data, dims = 1:max_dims, reduction = "pca_unint"))
dev.off()



#### Unintegrated UMAP ####


## Run UMAP
sc.data = FindNeighbors(object = sc.data, dims = 1:max_dims, reduction = "pca_unint")
sc.data = FindClusters(object = sc.data, cluster.name = "unintegrated_clusters")
sc.data = RunUMAP(object = sc.data, dims = 1:max_dims, reduction = "pca_unint", reduction.name = "umap_unint")

# active ident
sc.data = SetIdent(sc.data, value = "unintegrated_clusters")

# standard UMAPs
p1 = DimPlot(sc.data, reduction = "umap_unint", group.by = "sample_name", pt.size = 2.5, label = TRUE)
p2 = DimPlot(sc.data, reduction = "umap_unint", group.by = "sample_group", pt.size = 2.5, label = TRUE)
p3 = DimPlot(sc.data, reduction = "umap_unint", group.by = "unintegrated_clusters", pt.size = 2.5, label = TRUE)
p4 = DimPlot(sc.data, reduction = "umap_unint", group.by = "named_clusters_cell_type", pt.size = 2.5, label = TRUE)

ggp = p1 + p2 + p3 + p4
png(paste0(qc_folder, "/QC_RNA_unintegrated_UMAP.png"), height = 4000, width = 4000)
print(ggp)
dev.off()



#### Integrate Data ####

## set to RNA
DefaultAssay(sc.data) = 'RNA'

## Integrate - using harmony
sc.data = RunHarmony(sc.data,
                     group.by.vars = c("sample_name"),
                     theta = c(2),
                     reduction = "pca_unint", assay.use = "RNA", reduction.save = "integrated_pca")




#### Integrated UMAP ####


# re-cluster after integration
sc.data = FindNeighbors(sc.data, reduction = "integrated_pca", dims = 1:max_dims)
sc.data = FindClusters(sc.data, resolution = 0.75, verbose =TRUE, cluster.name = "integrated_clusters" )
sc.data = RunUMAP(sc.data, spread = 6, min.dist = 0.25, dims = 1:max_dims, reduction = "integrated_pca", reduction.name ="integrated_umap")

# Integrated UMAP
p1 = DimPlot(sc.data, reduction = "integrated_umap", group.by = "sample_name", pt.size = 1)
p2 = DimPlot(sc.data, reduction = "integrated_umap", group.by = "sample_group", pt.size = 1)
p3 = DimPlot(sc.data, reduction = "integrated_umap", group.by = "integrated_clusters", pt.size = 1, label = TRUE)
p4 = DimPlot(sc.data, reduction = "integrated_umap", group.by = "named_clusters_cell_type", pt.size = 1, label = TRUE)

ggp = p1 + p2 + p3 + p4
png(paste0(qc_folder, "/QC_RNA_integrated_UMAP.png"), height = 4000, width = 4000)
print(ggp)
dev.off()



#### Integrated Spatial ####

# plots
sc.data.sub = subset(sc.data, orig.ident == "AIH_1607")
png(paste0(qc_folder, "/AIH_1607_segmentation.png"), width = 4000, height = 4000)
ImageDimPlot(sc.data.sub)
dev.off()

sc.data.sub = subset(sc.data, orig.ident == "AIH_6079")
png(paste0(qc_folder, "/AIH_6079_segmentation.png"), width = 4000, height = 4000)
ImageDimPlot(sc.data.sub)
dev.off()

sc.data.sub = subset(sc.data, orig.ident == "SN_606")
png(paste0(qc_folder, "/SN_606_segmentation.png"), width = 4000, height = 4000)
ImageDimPlot(sc.data.sub)
dev.off()

sc.data.sub = subset(sc.data, orig.ident == "SN_1042")
png(paste0(qc_folder, "/SN_1042_segmentation.png"), width = 4000, height = 4000)
ImageDimPlot(sc.data.sub)
dev.off()

sc.data.sub = subset(sc.data, orig.ident == "SN_2739")
png(paste0(qc_folder, "/SN_2739_segmentation.png"), width = 4000, height = 4000)
ImageDimPlot(sc.data.sub)
dev.off()

sc.data.sub = subset(sc.data, orig.ident == "D_7678")
png(paste0(qc_folder, "/D_7678_segmentation.png"), width = 4000, height = 4000)
ImageDimPlot(sc.data.sub)
dev.off()

sc.data.sub = subset(sc.data, orig.ident == "D_6414")
png(paste0(qc_folder, "/D_6414_segmentation.png"), width = 4000, height = 4000)
ImageDimPlot(sc.data.sub)
dev.off()

sc.data.sub = subset(sc.data, orig.ident == "D_6446")
png(paste0(qc_folder, "/D_6446_segmentation.png"), width = 4000, height = 4000)
ImageDimPlot(sc.data.sub)
dev.off()


#### Cluster markers ####

## set to RNA
DefaultAssay(sc.data) = 'RNA'

# get markers
sc.data = JoinLayers(sc.data)
markers_RNA = FindAllMarkers(sc.data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers_RNA %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC) -> top15_RNA

# save markers
write.table(markers_RNA, paste0(markers_folder, "/RNA_cluster_markers.csv"), quote=FALSE, sep = "\t")
write.table(top15_RNA, paste0(markers_folder, "/RNA_cluster_markers_top15.csv"), quote=FALSE, sep = "\t")

# Plot clusters heatmap
png(paste0(markers_folder, "/RNA_cluster_markers_heatmap.png"), height = 4000, width = 2000)
DoHeatmap(sc.data, features = top15_RNA$gene) + NoLegend() + theme(text = element_text(size = 20))
dev.off()



#### Marker Feature Plots ####

# get the clusters
clusters = as.vector(unique(sc.data@meta.data$integrated_clusters))

# loop through the clusters
for (cluster_index in 1:length(clusters))
{
  # current cluster
  cluster_name = clusters[cluster_index]
  print(cluster_name)
  
  # get markers
  top15_RNA_cluster = subset(top15_RNA, cluster == cluster_name)$gene
  
  # test for enough RNA markers
  if (length(top15_RNA_cluster) > 0)
  {
    # make feature plot RNA
    ggp = FeaturePlot(sc.data, features = top15_RNA_cluster, reduction = "integrated_umap", pt.size = 3)
    png(paste0(markers_folder, "/marker_UMAP/cluster_", cluster_name, "_markers.png"), height = 4000, width = 4000)
    print(ggp)
    dev.off()
    
  }
}




#### Save RDS ####

saveRDS(sc.data, file="analysis/combined_all/all_integrated.RDS")


