#### Load Libraries ####

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(devtools)
library(hdf5r)

#### Guides and Manuals ####


# https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_07_spatial.html
# https://satijalab.org/seurat/articles/spatial_vignette.html
# https://www.livercellatlas.org/download.php
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04574-5


####--------------------------- PART 1. MERGED -------------------------------------####


#### Load cirrhotic single cell dataset ####

## wd
setwd("C:/Users/JohnJ/Desktop/analysis/Amber_Bozward_Visium_2023/analysis/public/GSE158723/")


## load
WF_158 = Read10X(data.dir = "WF_158/")
WF_158 = CreateSeuratObject(counts = WF_158, project = "WF_158", min.cells = 3, min.features = 200)
NPF1_156 = Read10X(data.dir = "NPF1_156/")
NPF1_156 = CreateSeuratObject(counts = NPF1_156, project = "NPF1_156", min.cells = 3, min.features = 200)
NPF1_158 = Read10X(data.dir = "NPF1_158/")
NPF1_158 = CreateSeuratObject(counts = NPF1_158, project = "NPF1_158", min.cells = 3, min.features = 200)
NPF2_156 = Read10X(data.dir = "NPF2_156/")
NPF2_156 = CreateSeuratObject(counts = NPF2_156, project = "NPF2_156", min.cells = 3, min.features = 200)
NPF2_158 = Read10X(data.dir = "NPF2_158/")
NPF2_158 = CreateSeuratObject(counts = NPF2_158, project = "NPF2_158", min.cells = 3, min.features = 200)

PF_158 = Read10X(data.dir = "PF_158/")
PF_158 = CreateSeuratObject(counts = PF_158, project = "PF_158", min.cells = 3, min.features = 200)
PF1_156 = Read10X(data.dir = "PF1_156/")
PF1_156 = CreateSeuratObject(counts = PF1_156, project = "PF1_156", min.cells = 3, min.features = 200)
PF2_156 = Read10X(data.dir = "PF2_156/")
PF2_156 = CreateSeuratObject(counts = PF2_156, project = "PF2_156", min.cells = 3, min.features = 200)




## merge
sc.merge = merge(WF_158,c(NPF1_156,NPF1_158,NPF2_156,NPF2_158,PF_158,PF1_156,PF2_156), merge.data = TRUE, add.cell.ids = c("WF_158","NPF1_156","NPF1_158","NPF2_156","NPF2_158","PF_158","PF1_156","PF2_156"))

# tidy
rm(WF_158,NPF1_156,NPF1_158,NPF2_156,NPF2_158,PF_158,PF1_156,PF2_156)
gc()


#### QC ####

## Percent mt
sc.merge[["percent.mt"]] = PercentageFeatureSet(sc.merge, pattern = "^MT-")

## QC violin
plot = VlnPlot(sc.merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
png("QC_violin_plots.png", height = 750, width = 1500)
print(plot)
dev.off()

## Feature scatters
plot = FeatureScatter(sc.merge, feature1 = "nCount_RNA", feature2 = "percent.mt")
png("QC_scatter_count_vs_mt.png", height = 750, width = 1500)
print(plot)
dev.off()

plot = FeatureScatter(sc.merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
png("QC_scatter_count_vs_rna.png", height = 750, width = 1500)
print(plot)
dev.off()

## Filter cells
sc.merge = subset(sc.merge, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 25)


#### Transform ####

sc.merge = SCTransform(sc.merge, assay = "RNA", verbose = TRUE, method = "poisson")


#### Variable Features ####


sc.merge = FindVariableFeatures(sc.merge, selection.method = "SCT", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 = head(VariableFeatures(sc.merge), 10)

# plot variable features with and without labels
plot = VariableFeaturePlot(sc.merge)
png("QC_varibale_features.png", height = 750, width = 1500)
print(plot)
dev.off()

plot = LabelPoints(plot = plot, points = top10, repel = TRUE)
png("QC_varibale_features_labelled.png", height = 750, width = 1500)
print(plot)
dev.off()


#### Dimensionality Reduction - Merged ####

## PCA
sc.merge = RunPCA(sc.merge, features = VariableFeatures(object = sc.merge))

plot = DimPlot(sc.merge, reduction = "pca")
png("DR_merged_PCA_scatter.png", height = 1000, width = 1000)
print(plot)
dev.off()

png("DR_merged_PCA_heatmap.png", height = 2500, width = 1500)
DimHeatmap(sc.merge, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()


## Elbow

plot = ElbowPlot(sc.merge)
png("DR_merged_elbow.png", height = 1000, width = 1000)
print(plot)
dev.off()


## Cluster

sc.merge = FindNeighbors(sc.merge, dims = 1:12)
sc.merge = FindClusters(sc.merge, resolution = 0.5)

## UMAP

sc.merge = RunUMAP(sc.merge, dims = 1:12)

# clusters umap
plot = DimPlot(sc.merge, reduction = "umap", group.by = c("ident", "orig.ident"))
png("DR_merged_UMAP_cluster.png", height = 1000, width = 1000)
print(plot)
dev.off()



####--------------------------- PART 2. INTEGRATION -------------------------------------####


## Tidy previous analysis
rm(sc.merge)
gc()



#### Load Data ####


## function to load and filter data
load_sc = function(sc_path, sc_name)
{
  sc_data = Read10X(data.dir = sc_path)
  sc_data = CreateSeuratObject(counts = sc_data, project = sc_name, min.cells = 3, min.features = 200)
  sc_data[["percent.mt"]] = PercentageFeatureSet(sc_data, pattern = "^MT-")
  sc_data = subset(sc_data, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 25)
  return(sc_data)
  
}

## Load and filter data
WF_158 = load_sc("WF_158/","WF_158")
NPF1_156 = load_sc("NPF1_156/","NPF1_156")
NPF1_158 = load_sc("NPF1_158/","NPF1_158")
NPF2_156 = load_sc("NPF2_156/","NPF2_156")
NPF2_158 = load_sc("NPF2_158/","NPF2_158")
PF_158 = load_sc("PF_158/","PF_158")
PF1_156 = load_sc("PF1_156/","PF1_156")
PF2_156 = load_sc("PF2_156/","PF2_156")


#### Integrate ####

# create list
sc.list = list(WF_158 = WF_158, NPF1_156 = NPF1_156, NPF1_158 = NPF1_158, NPF2_156 = NPF2_156, NPF2_158 = NPF2_158, PF_158 = PF_158, PF1_156 = PF1_156, PF2_156 = PF2_156)

# run SCT on all datasets
sc.list = lapply(sc.list, SCTransform, assay = "RNA", method = "poisson")

# Tidy up
rm(WF_158,NPF1_156,NPF1_158,NPF2_156,NPF2_158,PF_158,PF1_156,PF2_156)
gc()


# need to set maxSize for PrepSCTIntegration to work
options(future.globals.maxSize = 4000 * 1024^2)  # set allowed size to 2K MiB
st.features = SelectIntegrationFeatures(sc.list, nfeatures = 2000, verbose = FALSE)
sc.list = PrepSCTIntegration(object.list = sc.list, anchor.features = st.features,verbose = FALSE)

# get anchors  - pf2_156 (8) and pf_158 (6)
int.anchors = FindIntegrationAnchors(object.list = sc.list, normalization.method = "SCT",verbose = FALSE, anchor.features = st.features,  reference = c(6, 8))

# tidy up
rm(sc.list)
gc()

# integrate
sc.integrated = IntegrateData(anchorset = int.anchors, normalization.method = "SCT",verbose = FALSE)

# tidy up 
rm(int.anchors)
gc()


#### Dimensionality Reduction - Integrated ####


sc.integrated = ScaleData(sc.integrated, verbose = FALSE)
sc.integrated = RunPCA(sc.integrated, verbose = FALSE)
sc.integrated = FindNeighbors(sc.integrated, reduction = "pca", dims = 1:30)
sc.integrated = FindClusters(sc.integrated, verbose = FALSE)
sc.integrated = RunUMAP(sc.integrated,  reduction = "pca", dims = 1:30)

# Plot PCA heatmap
png("DR_integrated_pca_heatmap.png", height = 1500, width = 1500)
DimHeatmap(sc.integrated, dims = 1:15, cells = 2500, balanced = TRUE) + theme(text = element_text(size = 20))
dev.off()


# Plot UMAP
plot = DimPlot(sc.integrated, reduction = "umap", group.by = c("ident", "orig.ident"))
png("DR_integrated_cluster_umap.png", height = 1250, width = 2500)
print(plot)
dev.off()


# Plot clusters heatmap
sc.integrated = PrepSCTFindMarkers(sc.integrated)
liver.markers = FindAllMarkers(sc.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
liver.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

png("DR_integrated_clusters_heatmap.png", height = 2000, width = 2000)
DoHeatmap(sc.integrated, features = top10$gene) + NoLegend() + theme(text = element_text(size = 20))
dev.off()



#### Linneages plots ####


# LP
plot = FeaturePlot(sc.integrated, features = c("IGKC","JCHAIN", "IGHA2","IGHM","KLRB1"), cols = c("purple","black","yellow"))
png("Linneages_LP.png", height = 1250, width = 2500)
print(plot)
dev.off()

# MP
plot = FeaturePlot(sc.integrated, features = c("HLA-DPA1","CD74", "HLA-DRA", "C1QA", "C1QB"), cols = c("purple","black","yellow"))
png("Linneages_MP.png", height = 1250, width = 2500)
print(plot)
dev.off()

# HSC
plot = FeaturePlot(sc.integrated, features = c("CXCL14","IGFBP3", "IGFBP7", "DCN", "COLEC11"), cols = c("purple","black","yellow"))
png("Linneages_HSC.png", height = 1250, width = 2500)
print(plot)
dev.off()

# EC
plot = FeaturePlot(sc.integrated, features = c("FCN2","FCN3", "AKAP12", "CRHBP", "DNASE1L3"), cols = c("purple","black","yellow"))
png("Linneages_EC.png", height = 1250, width = 2500)
print(plot)
dev.off()

# CC
plot = FeaturePlot(sc.integrated, features = c("TACSTD2","KRT19","CXCL1","KRT7","CXCL8"), cols = c("purple","black","yellow"))
png("Linneages_CC.png", height = 1250, width = 2500)
print(plot)
dev.off()

# HEP
plot = FeaturePlot(sc.integrated, features = c("TF","ORM2","HP","SAA1","ORM1"), cols = c("purple","black","yellow"))
png("Linneages_HEP.png", height = 1250, width = 2500)
print(plot)
dev.off()


####--------------------------- PART 3. CELL TYPES -------------------------------------####


#### Get Cell Types ####

## get umap table for ggplot
umap_tx = data.frame(sc.integrated@reductions$umap@cell.embeddings)
umap_tx$ident = sc.integrated@active.ident


# plot per cluster umaps

ggp = ggplot(umap_tx, aes(x=UMAP_1, y=UMAP_2, colour = ident)) + geom_point() + theme(legend.position = "none") + xlim(-15,15) + ylim(-18,8)
png("cell_types/UMAP.png", height = 1250, width = 1250)
print(ggp)
dev.off()

for (cluster in 0:max(as.numeric(umap_tx$ident)))
{
  ggp = ggplot(subset(umap_tx, ident == cluster), aes(x=UMAP_1, y=UMAP_2)) + geom_point() + xlim(-15,15) + ylim(-18,8)
  png(paste0("cell_types/UMAP_",cluster,".png"), height = 1250, width = 1250)
  print(ggp)
  dev.off()
}



## manually match
new.cluster.ids = c("Hep.1", "Hep.2", "Hep.3", "Hep.4", "Hep.5", "Hep.6", "Hep.7", "MP.1", "MP.2", "Hep.8", "Hep.9", "CC.1", "EC.1", "Hep.10", "LP.1", "LP.2", "EC.2", "MP.3", "MP.4", "HSC.1", "HSC.2", "LP.3", "LP.4", "LP.5", "EC.3")

# rename
names(new.cluster.ids) = levels(sc.integrated)
sc.integrated = RenameIdents(sc.integrated, new.cluster.ids)

# plot
plot = DimPlot(sc.integrated, reduction = "umap", label = TRUE, pt.size = 0.5)
png("DR_integrated_cell_type_umap.png", height = 1250, width = 1250)
print(plot)
dev.off()


# create new column for linneage
umap_tx = data.frame(sc.integrated@reductions$umap@cell.embeddings)
umap_tx$ident = sc.integrated@active.ident
umap_tx$linneage = as.character(umap_tx$ident)
for (row in 1:nrow(umap_tx))
{
  umap_tx[row,"linneage"] = unlist(strsplit(as.character(umap_tx[row,"ident"]),split="\\."))[1]
}

# plot
ggp = ggplot(umap_tx, aes(x=UMAP_1, y=UMAP_2, colour = linneage)) + geom_point()
png("DR_integrated_linneage_umap.png", height = 1250, width = 1250)
print(ggp)
dev.off()


#### Save integrated table ####

#saveRDS(sc.integrated, file = "single_cell_integrated.rds")



#### Get markers for named clusters ####

# load 
sc.integrated = readRDS(file = "single_cell_integrated.rds")

# Plot clusters heatmap
sc.integrated = PrepSCTFindMarkers(sc.integrated)
liver.markers = FindAllMarkers(sc.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
liver.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

png("DR_integrated_clusters_named_heatmap.png", height = 2000, width = 2000)
DoHeatmap(sc.integrated, features = top10$gene) + NoLegend() + theme(text = element_text(size = 20))
dev.off()


# get table of top 25 markers
liver.markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC) -> top25
write.table(top25, file = "top_25_markers.csv", sep = "\t", quote = FALSE)
write.table(liver.markers, file = "all_markers.csv", sep = "\t", quote = FALSE)



