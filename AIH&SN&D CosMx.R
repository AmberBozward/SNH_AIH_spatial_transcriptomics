# loading the required libraries

#install.packages("scCustomize")
#install.packages("ggpubr")
#install.packages("ggplot2")
#install.packages("clustree")
#install.packages('devtools')
#devtools::install_github('immunogenomics/presto')
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("SingleR")
#devtools::install_github('satijalab/seurat-data')
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("clusterProfiler")
#install.packages("msigdbr")
#if (!requireNamespace('remotes', quietly = TRUE) {
#install.packages('remotes')
#remotes::install_github('satijalab/azimuth', ref = 'master')
#BiocManager::install("dittoSeq")
#install.packages("harmony")
#BiocManager::install("scDblFinder")
#BiocManager::install("GOSemSim")
#BiocManager::install("enrichplot")
#install.packages("ggupset")
#BiocManager::install("DESeq2")
  
library(dittoSeq)
library(harmony)
library(scDblFinder)
library(GOSemSim)
library(enrichplot)
library(ggupset)
library(DESeq2)
library(msigdbr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(SeuratData)
library(Seurat)
library(dplyr)
library(tidyr)
library(patchwork)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(scCustomize)
library(clustree)
library(presto)
library(ggrepel)
library(SingleR)

# working directory
setwd("/Users/bozwarag/Documents/Amber Bozward - Mac/CosMx data download/AIH & Seroneg matched/CosMX analysis/10.12.24/Plots made by me")

# dot plot for specific markers 
markers=c("CD68","CD4", "CD19")

subset_processing=function(seurat.object)
{
  # normalization
  seurat.object=NormalizeData(seurat.object)
  
  # Finding the most variable features (genes)
  seurat.object=FindVariableFeatures(seurat.object)
  
  # Scaling 
  seurat.object <- ScaleData(seurat.object)
  
  # Performing PCA 
  seurat.object <- RunPCA(seurat.object)
  
  return(seurat.object)
}

# function for running UMAP and clustering based on input PC dimensions
UMAP_clustering=function(seurat.object,cluster_name)
{
  seurat.object <- FindNeighbors(seurat.object, dims = 1:15, reduction="pca_unint")
  seurat.object <- FindClusters(seurat.object, cluster.name=cluster_name)
  seurat.object <- RunUMAP(seurat.object, dims = 1:15, reduction="pca_unint",reduction.name="umap_unint")
  
  return(seurat.object)
}

# harmony integration v5
integrate_harmony_v5=function(seurat.object,theta1,parameter)
{
  seurat.object <- IntegrateLayers(object = seurat.object, method = HarmonyIntegration,
                                   orig.reduction = "pca_unint", new.reduction = "integrated_pca", group.by.vars=c(parameter),
                                   theta=theta1)
  seurat.object <- FindNeighbors(seurat.object, reduction = "integrated_pca")
  seurat.object <- FindClusters(seurat.object, cluster.name = "harmony_clusters")
  seurat.object <- RunUMAP(seurat.object, reduction = "integrated_pca", dims = 1:15, reduction.name = "integrated_umap")
  
  return(seurat.object)
}


#### Plotting the UMAP for all cells ####
# loading the seurat object
object=readRDS("/Users/bozwarag/Documents/Amber Bozward - Mac/CosMx data download/AIH & Seroneg matched/CosMX analysis/10.12.24/final.RDS")
View(object@meta.data)

Idents(object)=object@meta.data$final_cell_types
DimPlot(object, reduction="integrated_umap",label=TRUE)

#### plotting genes of interest on the UMAP ####
ggp=FeaturePlot_scCustom(object, features = c("TIGIT"), reduction = "integrated_umap",label=TRUE)
ggp

ggsave("name.png",ggp,height = 2500,width = 2500,units = "px",dpi = 200) # for saving the plot

ggp + scale_color_gradient(low = "lightgrey", high = "red3") 

ggp = FeaturePlot_scCustom(
  object, 
  features = c("MX1"), 
  reduction = "integrated_umap", 
  label = TRUE,
  repel=TRUE,
  min.cutoff = "q30",  # Adjust minimum cutoff (e.g., 10th percentile),
  max.cutoff = "q95"
)

# dot plot for markers in "markers"
markers=c("CLEC2D")
plot=DotPlot_scCustom(object, features= markers) + NoLegend() +
  theme(axis.text.x = element_text(angle = 90), axis.text.y = element_text(size=9)) + xlab('Genes') + ylab('Clusters') +
  coord_flip() + ggtitle("") + theme(plot.title = element_text(hjust = 0.5, size = 14))

plot

# Plotting a marker between patient cohorts #

# CD63
colnames(object@meta.data)

cd63_data <- FetchData(object, vars = c("CD63", "sample_name")) 
cd63_data$sample_name <- as.factor(cd63_data$sample_name)

head(cd63_data)

ggplot(cd63_data, aes(x = sample_name, y = CD63, fill = sample_name)) +
  geom_boxplot() +
  labs(title = "CD63 Expression Across Conditions", 
       x = "sample_name", 
       y = "CD63 Expression Level") +
  theme_minimal() +
  theme(legend.position = "none")

# TNFSF14
TNFSF14_data <- FetchData(object, vars = c("TNFSF14", "sample_name")) 
TNFSF14_data$sample_name <- as.factor(TNFSF14_data$sample_name)

ggplot(TNFSF14_data, aes(x = sample_name, y = TNFSF14, fill = sample_name)) +
  geom_boxplot() +
  labs(title = "TNFSF14 Expression Across Conditions", 
       x = "sample_name", 
       y = "TNFSF14 Expression Level") +
  theme_minimal() +
  theme(legend.position = "none")

### Number of cells in each subset per samples ####

# number of cells per cell type per sample group
cts.sample=as.data.frame.matrix(table(Idents(object), object$sample_group), margin = 2)
cts.sample=cts.sample %>% 
  rownames_to_column("cluster")
write.table(cts.sample,file="Allcells_cellcounts_persamplegroup.tsv",sep="\t",row.names = FALSE)

cts.m=reshape2::melt(cts.sample)
cts.m$cluster=as.factor(cts.m$cluster)

# number of cells per cell type per sample
cts.sample=as.data.frame.matrix(table(Idents(object), object$sample_name), margin = 2)
cts.sample=cts.sample %>% 
  rownames_to_column("cluster")
write.table(cts.sample,file="Allcells_cellcounts_persample.tsv",sep="\t",row.names = FALSE)

# Create a stacked bar plot
ggp=ggplot(cts.m, aes(x = cluster, y = value, fill = variable)) +
  geom_bar(stat = "identity") + 
  labs(x = "Cell cluster", y = "Count", title = "Stacked Bar Plot of Cell Types") +
  my_theme

# Trying to put the sample name on the x axis ## This is the one used for the manuscript

cts.sample=as.data.frame.matrix(table(Idents(object), object$sample_name), margin = 2)
cts.sample=cts.sample %>% 
  rownames_to_column("cluster")
cts.m$sample <- as.factor(cts.m$sample)

cts.m <- cts.sample %>% 
  pivot_longer(cols = -cluster, names_to = "sample", values_to = "value")

ggp = ggplot(cts.m, aes(x = sample, y = value, fill = cluster)) +
  geom_bar(stat = "identity", position = "stack") + 
  labs(x = "Sample", y = "Cell Count", title = "Stacked Bar Plot of Cell Types Per Sample", fill = "Cell Type") +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels if needed

# Plotting % total cells - stacked plot 

cts.m <- cts.sample %>%
  pivot_longer(cols = -cluster, names_to = "sample", values_to = "count") %>%
  group_by(sample) %>%
  mutate(percentage = (count / sum(count)) * 100) %>%
  ungroup()

ggp = ggplot(cts.m, aes(x = sample, y = percentage, fill = cluster)) +
  geom_bar(stat = "identity", position = "stack") + 
  labs(x = "Sample", y = "Percentage of Total Cells", 
       title = "Proportion of Cell Types Per Sample", 
       fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels

#### Hepatocytes ####
# creating this subset from the final dataset

# Subset

#Hepatocytes = subset(object, subset = final_cell_types == "Hepatocyte")

#Hepatocytes <- NormalizeData(Hepatocytes)
#Hepatocytes <- FindVariableFeatures(Hepatocytes, selection.method = "vst", nfeatures = 2000)
#Hepatocytes <- ScaleData(Hepatocytes)
#Hepatocytes <- RunPCA(Hepatocytes, npcs = 30)
#ElbowPlot(Hepatocytes)
#Hepatocytes <- FindNeighbors(Hepatocytes, dims = 1:15)  # Adjust PCs based on ElbowPlot
#Hepatocytes <- FindClusters(Hepatocytes, resolution = 0.5)  # Tune resolution for cluster granularity
#Hepatocytes <- RunUMAP(Hepatocytes, dims = 1:15)

Hepatocytes=readRDS("/Users/bozwarag/Documents/Amber Bozward - Mac/CosMx data download/AIH & Seroneg matched/CosMX analysis/10.12.24/Plots made by me/Hepatocytes.RDS")

DimPlot(Hepatocytes, reduction = "umap", label = TRUE)

View(Hepatocytes@meta.data)

Idents(Hepatocytes) <- "seurat_clusters"


#plotting markers to define hepatocyte sub-clusters

pericentral=c("LGR5")

periportal=c("ARG1", "SCD", "HMGCS1", "ACSS2", "TM7SF2", "TMEM97", "CP", "CRP", "SLPI", "C2ORF82", "ACAT2", "TM4SF5", "MSMO1", "LEPR")


plot=DotPlot_scCustom(Hepatocytes, features= centralvenous) + NoLegend() +
  theme(axis.text.x = element_text(angle = 90), axis.text.y = element_text(size=9)) + xlab('Genes') + ylab('Clusters') +
  coord_flip() + ggtitle("title of the plot") + theme(plot.title = element_text(hjust = 0.5, size = 14))

plot

# find top 20 genes in each cluster

markers <- FindAllMarkers(Hepatocytes, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top20 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

top20_list <- split(top20$gene, top20$cluster)

write.csv(top20, "Top20_Genes_Per_Cluster.csv", row.names = FALSE)

# Annotations

new.cluster.ids=c("IS hepatocytes", "HMGCS1+ periportal", "IS hepatocytes", "HMGCS1+ periportal", "IS hepatocytes", "FN1+ hepatocytes", "MX1+ hepatocytes", "MX1+ hepatocytes", "IS hepatocytes", "FN1+ hepatocytes", "FN1+ hepatocytes", "LGR5+ pericentral", "LGR5+ pericentral", "IS hepatocytes")

names(new.cluster.ids)=levels(Hepatocytes)
Hepatocytes= RenameIdents(Hepatocytes,new.cluster.ids)

ggp=DimPlot(Hepatocytes, reduction="umap",label=TRUE,label.size = 5, pt.size = 1.5, repel=TRUE) + ggtitle("Hepatocytes") +
  my_theme + theme(plot.title = element_text(hjust = 0.5, size = 14)) + NoLegend()


# number of cells per cell type
cts.sample=as.data.frame.matrix(table(Idents(Hepatocytes), Hepatocytes$sample_group), margin = 2)
cts.sample=cts.sample %>% 
  rownames_to_column("cluster")
write.table(cts.sample,file="Hepatocytes_cellcounts_persample.tsv",sep="\t",row.names = FALSE)

cts.m=reshape2::melt(cts.sample)
cts.m$cluster=as.factor(cts.m$cluster)

# Create a stacked bar plot
ggp=ggplot(cts.m, aes(x = cluster, y = value, fill = variable)) +
  geom_bar(stat = "identity") + 
  labs(x = "Cell cluster", y = "Count", title = "Stacked Bar Plot of Cell Types") +
  my_theme

# Swapping so samples are on the x axis
cts.sample=as.data.frame.matrix(table(Idents(Hepatocytes), Hepatocytes$sample_group), margin = 2)
cts.sample=cts.sample %>% 
  rownames_to_column("cluster")

cts.m <- cts.sample %>% 
  pivot_longer(cols = -cluster, names_to = "sample", values_to = "value")

ggp = ggplot(cts.m, aes(x = sample, y = value, fill = cluster)) +
  geom_bar(stat = "identity", position = "stack") + 
  labs(x = "Sample", y = "Cell Count", title = "Stacked Bar Plot of Cell Types Per Sample", fill = "Cell Type") +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels if needed


# Plotting % total cells - stacked plot 

cts.m <- cts.sample %>%
  pivot_longer(cols = -cluster, names_to = "sample", values_to = "count") %>%
  group_by(sample) %>%
  mutate(percentage = (count / sum(count)) * 100) %>%
  ungroup()

ggp = ggplot(cts.m, aes(x = sample, y = percentage, fill = cluster)) +
  geom_bar(stat = "identity", position = "stack") + 
  labs(x = "Sample", y = "Percentage of Total Cells", 
       title = "Proportion of Cell Types Per Sample", 
       fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels

# Plotting a dot plot of the defining markers

genes.of.interest=c("LGR5","MX1", "FN1", "HMGCS1", "CD163", "TLR2")

plot=DotPlot_scCustom(Hepatocytes, features= genes.of.interest) +
  theme(axis.text.x = element_text(angle = 90), axis.text.y = element_text(size=9)) + xlab('Genes') + ylab('Clusters') +
  ggtitle("") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) + 
  my_theme

plot


saveRDS(Hepatocytes, file="/Users/bozwarag/Documents/Amber Bozward - Mac/CosMx data download/AIH & Seroneg matched/CosMX analysis/10.12.24/Plots made by me/Hepatocytes.RDS")


#### Macrophages ####

Macrophages=readRDS("/Users/bozwarag/Documents/Amber Bozward - Mac/CosMx data download/AIH & Seroneg matched/CosMX analysis/10.12.24/Plots made by me/Macrophages_edited.RDS")
View(Macrophages@meta.data)

DimPlot(Macrophages, reduction="integrated_umap",label=TRUE)

ggp=FeaturePlot_scCustom(Macrophages, features = c("CD63","C1QA"), reduction = "integrated_umap",label=TRUE, label.size = 3, pt.size = 0.5)
ggp


markers=c("CD63")
plot=DotPlot_scCustom(Macrophages, features= markers) + NoLegend() +
  theme(axis.text.x = element_text(angle = 90), axis.text.y = element_text(size=9)) + xlab('Genes') + ylab('Clusters') +
  coord_flip() + ggtitle("CD63 expression on macrophage subsets") + theme(plot.title = element_text(hjust = 0.5, size = 14))

plot


#getting the clustree
# resolution 0.2-1.6
for(i in c(1.6)){
  Macrophages <- Seurat::FindClusters(Macrophages, resolution = i)
}

View(Macrophages@meta.data)

# clustree
ggp=clustree(Macrophages, prefix = "R NA_snn_res.", 
             node_colour = "sc3_stability", use_core_edges = T, 
             node_label = NULL, node_label_aggr = NULL)   +
  scale_edge_color_continuous(low = "blue", high = "red") +
  guides(edge_colour = FALSE, edge_alpha = FALSE) +
  theme(legend.position = "bottom") 
ggsave("clustree.png",ggp,height = 2500,width = 2500,units = "px",dpi = 200)

Macrophages=SetIdent(Macrophages,value=Macrophages@meta.data$RNA_snn_res.1.2)

#plotting markers from Calli to define macrophages 

markers.myeloid=c("S100A8","S100A9","CD14","C1QA","C1QC","MARCO","CD68","MRC1","GPNMB","APOE","MERTK",
                  "SPP1","CD1C","CLEC10A","LAMP3","GZMB","IL1B","IL18",
                  "TNF","IL6","MX1","IFIT1","CXCL8","THBS1","NLRP3","NRG1", "LGALS3", "FABP5", "ANXA2", "CD63", "LYZ", "CD163", "TOX")

Mregs.markers=c("CD80", "PDCDL1", "MRC1", "CD86", "HLA-DRB1", "DHRS9")

granzyme.markers=c("GZMK", "GZMB", "GZMA")

inflammatory.markers=c("IGNG", "TLR4", "NOS2")

anti_inflammatory.markers=c("MRC1", "CD9", "IL4", "IL10", "IL13")

macrophage.markers=c("APOC1", "THBS1", "LYZ")

M1.markers=c("IL1A", "IL1B", "IL6", "TLR2", "TLR4", "CD80", "CD86")

M2.markers=c("PPARG", "ARG1", "CD163")

plot=DotPlot_scCustom(Macrophages, features= M2.markers) + NoLegend() +
  theme(axis.text.x = element_text(angle = 90), axis.text.y = element_text(size=9)) + xlab('Genes') + ylab('Clusters') +
  coord_flip() + ggtitle("title of the plot") + theme(plot.title = element_text(hjust = 0.5, size = 14))


plot=DotPlot_scCustom(Macrophages, features= inflammatory.markers) + NoLegend() +
  theme(axis.text.x = element_text(angle = 90), axis.text.y = element_text(size=9)) + xlab('Genes') + ylab('Clusters') +
  coord_flip() + ggtitle("pro-inflammatory") + theme(plot.title = element_text(hjust = 0.5, size = 14))


plot=DotPlot_scCustom(Macrophages, features= anti_inflammatory.markers) + NoLegend() +
  theme(axis.text.x = element_text(angle = 90), axis.text.y = element_text(size=9)) + xlab('Genes') + ylab('Clusters') +
  coord_flip() + ggtitle("anti-inflammatory") + theme(plot.title = element_text(hjust = 0.5, size = 14))

plot

# Annotations
new.cluster.ids=c("LYZhi macrophages", "Pro-inflammatory macrophages", "Mregs", "APOC1hi macrophages", "THSB1hi macrophages", "Mregs", "Anti-inflammatory macrophages", "Anti-inflammatory macrophages", "LYZhi macrophages", "soup", "THSB1hi macrophages")

names(new.cluster.ids)=levels(Macrophages)
Macrophages= RenameIdents(Macrophages,new.cluster.ids)

ggp=DimPlot(Macrophages, reduction="integrated_umap",label=TRUE,label.size = 5, pt.size = 1.5, repel=TRUE) + ggtitle("Macrophages") +
  my_theme + theme(plot.title = element_text(hjust = 0.5, size = 14)) + NoLegend()


# number of cells per cell type
cts.sample=as.data.frame.matrix(table(Idents(Macrophages), Macrophages$sample_group), margin = 2)
cts.sample=cts.sample %>% 
  rownames_to_column("cluster")
write.table(cts.sample,file="macrophage_cellcounts_persample.tsv",sep="\t",row.names = FALSE)

cts.m=reshape2::melt(cts.sample)
cts.m$cluster=as.factor(cts.m$cluster)

# Create a stacked bar plot
ggp=ggplot(cts.m, aes(x = cluster, y = value, fill = variable)) +
  geom_bar(stat = "identity") + 
  labs(x = "Cell cluster", y = "Count", title = "Stacked Bar Plot of Cell Types") +
  my_theme

# Swapping so samples are on the x axis
cts.sample=as.data.frame.matrix(table(Idents(Macrophages), Macrophages$sample_group), margin = 2)
cts.sample=cts.sample %>% 
  rownames_to_column("cluster")

cts.m <- cts.sample %>% 
  pivot_longer(cols = -cluster, names_to = "sample", values_to = "value")

ggp = ggplot(cts.m, aes(x = sample, y = value, fill = cluster)) +
  geom_bar(stat = "identity", position = "stack") + 
  labs(x = "Sample", y = "Cell Count", title = "Stacked Bar Plot of Cell Types Per Sample", fill = "Cell Type") +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels if needed


# Plotting % total cells - stacked plot 

cts.m <- cts.sample %>%
  pivot_longer(cols = -cluster, names_to = "sample", values_to = "count") %>%
  group_by(sample) %>%
  mutate(percentage = (count / sum(count)) * 100) %>%
  ungroup()

ggp = ggplot(cts.m, aes(x = sample, y = percentage, fill = cluster)) +
  geom_bar(stat = "identity", position = "stack") + 
  labs(x = "Sample", y = "Percentage of Total Cells", 
       title = "Proportion of Cell Types Per Sample", 
       fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels

# Plotting a dot plot of the defining markers

Macrophages@meta.data$new.cluster.ids=factor(Macrophages@meta.data$new.cluster.ids, levels=c("Hepatocyte", "IS Hepatocyte", "Hepatic stellate", "Epithelial", "Type 2 LSEC", "Type 1 LSEC", "Neutrophil", "Macrophage","Monocyte", "T cell", "NK cell", "Plasma cell", "B cell"))

Idents(Macrophages)=Macrophages@meta.data$new.cluster.ids

genes.of.interest=c("THBS1", "S100A9", "S100A8", "APOC1", "MARCO", "ANXA2", "FABP5", "LGALS3", "SPP1", "NLRP3", "CXCL8", "IFIT1", "MX1", "IL6", "TNF", "IL18", "IL1B", "CLEC10A", "MERTK", "LYZ", "APOE", "GPNMB", "MRC1", "C1QA", "CD14")

plot=DotPlot_scCustom(Macrophages, features= genes.of.interest) +
  theme(axis.text.x = element_text(angle = 90), axis.text.y = element_text(size=9)) + xlab('Genes') + ylab('Clusters') +
  ggtitle("") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) + 
  my_theme

plot

saveRDS(Macrophages, file="/Users/bozwarag/Documents/Amber Bozward - Mac/CosMx data download/AIH & Seroneg matched/CosMX analysis/10.12.24/Plots made by me/Macrophages_edited.RDS")


#### Renaming Macrophage subsets to M1 and M2 ####

Macrophages=readRDS("/Users/bozwarag/Documents/Amber Bozward - Mac/CosMx data download/AIH & Seroneg matched/CosMX analysis/10.12.24/Amber_CosMx_Nov_2024/subsets/Macrophage.RDS")

#getting the clustree
# resolution 0.2-1.6
#for(i in c(0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6)){
  object <- Seurat::FindClusters(object, resolution = i)
}


for(i in c(1.0)){
  Macrophages <- Seurat::FindClusters(Macrophages, resolution = i)
}

Macrophages=SetIdent(Macrophages,value=Macrophages@meta.data$RNA_snn_res.1)

M1.markers=c("IL1A", "IL1B", "IL6", "TLR2", "TLR4", "CD80", "CD86")
M2.markers=c("PPARG", "ARG1", "CD163")


plot=DotPlot_scCustom(Macrophages, features= M2.markers) + NoLegend() +
  theme(axis.text.x = element_text(angle = 90), axis.text.y = element_text(size=9)) + xlab('Genes') + ylab('Clusters') +
  coord_flip() + ggtitle("title of the plot") + theme(plot.title = element_text(hjust = 0.5, size = 14))


# Annotations
new.cluster.ids=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13")

new.cluster.ids=c("M1", "M1", "M1", "M1", "M2", "M1", "M1", "soup", "M1", "M2", "M2", "M2", "M2", "M2")


names(new.cluster.ids)=levels(Macrophages)
Macrophages= RenameIdents(Macrophages,new.cluster.ids)

ggp=DimPlot(Macrophages, reduction="integrated_umap",label=TRUE,label.size = 5, pt.size = 1.5, repel=TRUE) + ggtitle("Macrophages") +
  my_theme + theme(plot.title = element_text(hjust = 0.5, size = 14)) + NoLegend()


# number of cells per cell type
cts.sample=as.data.frame.matrix(table(Idents(Macrophages), Macrophages$sample_group), margin = 2)
cts.sample=cts.sample %>% 
  rownames_to_column("cluster")
write.table(cts.sample,file="macrophage_cellcounts_persample.tsv",sep="\t",row.names = FALSE)

cts.m=reshape2::melt(cts.sample)
cts.m$cluster=as.factor(cts.m$cluster)

# Create a stacked bar plot
ggp=ggplot(cts.m, aes(x = cluster, y = value, fill = variable)) +
  geom_bar(stat = "identity") + 
  labs(x = "Cell cluster", y = "Count", title = "Stacked Bar Plot of Cell Types") +
  my_theme

# Swapping so samples are on the x axis
cts.sample=as.data.frame.matrix(table(Idents(Macrophages), Macrophages$sample_group), margin = 2)
cts.sample=cts.sample %>% 
  rownames_to_column("cluster")

cts.m <- cts.sample %>% 
  pivot_longer(cols = -cluster, names_to = "sample", values_to = "value")

ggp = ggplot(cts.m, aes(x = sample, y = value, fill = cluster)) +
  geom_bar(stat = "identity", position = "stack") + 
  labs(x = "Sample", y = "Cell Count", title = "Stacked Bar Plot of Cell Types Per Sample", fill = "Cell Type") +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels if needed


# Plotting % total cells - stacked plot 

cts.m <- cts.sample %>%
  pivot_longer(cols = -cluster, names_to = "sample", values_to = "count") %>%
  group_by(sample) %>%
  mutate(percentage = (count / sum(count)) * 100) %>%
  ungroup()

ggp = ggplot(cts.m, aes(x = sample, y = percentage, fill = cluster)) +
  geom_bar(stat = "identity", position = "stack") + 
  labs(x = "Sample", y = "Percentage of Total Cells", 
       title = "Proportion of Cell Types Per Sample", 
       fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels


# To make a plot with each individual patient 

cts.sample=as.data.frame.matrix(table(Idents(Macrophages), Macrophages$sample_name), margin = 2)
cts.sample=cts.sample %>% 
  rownames_to_column("cluster")
write.table(cts.sample,file="macrophage_cellcounts_persample.tsv",sep="\t",row.names = FALSE)

cts.m=reshape2::melt(cts.sample)
cts.m$cluster=as.factor(cts.m$cluster)

# Create a stacked bar plot
ggp=ggplot(cts.m, aes(x = cluster, y = value, fill = variable)) +
  geom_bar(stat = "identity") + 
  labs(x = "Cell cluster", y = "Count", title = "Stacked Bar Plot of Cell Types") +
  my_theme

# Swapping so samples are on the x axis
cts.sample=as.data.frame.matrix(table(Idents(Macrophages), Macrophages$sample_name), margin = 2)
cts.sample=cts.sample %>% 
  rownames_to_column("cluster")

cts.m <- cts.sample %>% 
  pivot_longer(cols = -cluster, names_to = "sample", values_to = "value")

ggp = ggplot(cts.m, aes(x = sample, y = value, fill = cluster)) +
  geom_bar(stat = "identity", position = "stack") + 
  labs(x = "Sample", y = "Cell Count", title = "Stacked Bar Plot of Cell Types Per Sample", fill = "Cell Type") +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels if needed


# Plotting % total cells - stacked plot 

cts.m <- cts.sample %>%
  pivot_longer(cols = -cluster, names_to = "sample", values_to = "count") %>%
  group_by(sample) %>%
  mutate(percentage = (count / sum(count)) * 100) %>%
  ungroup()

ggp = ggplot(cts.m, aes(x = sample, y = percentage, fill = cluster)) +
  geom_bar(stat = "identity", position = "stack") + 
  labs(x = "Sample", y = "Percentage of Total Cells", 
       title = "Proportion of Cell Types Per Sample", 
       fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels

### To change the resolution - Macrophages ###

#get the top genes for the clusters 
## Marker genes

object_clean=JoinLayers(object_clean)
macrophage.markers=FindAllMarkers(object_clean, only.pos=TRUE, logfc.threshold = 0.25, min.pct=0.25)
top50.macrophage <- macrophage.markers %>%
  group_by(cluster) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 50)

# saving the marker genes table
write.table(top50.macrophage,file="top50_markers_cleaned_macrophagesres1.6.tsv",sep="\t",row.names = FALSE)


#highlight just one cluster on the UMAP 

ggp=Cluster_Highlight_Plot(object, cluster_name = c("14","1"), highlight_color = c("navy","green"), reduction="integrated_umap",
                           background_color = "lightgray") + my_theme

#removing soupy clusters 
object_clean=subset(object, idents="6", invert=TRUE)
View(object_clean@meta.data)
table(object_clean@meta.data$RNA_snn_res.1.6)

# cleaning up the seurat object
object_clean@meta.data=object_clean@meta.data[,-c(104:111)]

# Splitting the data into layers
object_clean[['RNA']]= split(object_clean[["RNA"]], f = object_clean$sample_name)

# re-scaling and re-clustering
object_clean=subset_processing(object_clean)
ElbowPlot(object_clean,ndims=40)

# clustering
object_clean=UMAP_clustering(object_clean,"fine_clusters")
DimPlot(object_clean, reduction="umap_unint",group.by="sample_name") 
View(object_clean@meta.data)

# integration by harmony
object_clean=integrate_harmony_v5(object_clean,2,'sample_name')
ggp=DimPlot(object_clean, reduction="integrated_umap", group.by="sample_name") + my_theme 
ggsave("integrated.png",ggp,height = 2500,width = 2500,units = "px",dpi = 200)

ggp=DimPlot(object_clean, reduction="integrated_umap", label=TRUE, label.size=6) + ggtitle("macrophages") + 
  my_theme + theme(plot.title = element_text(hjust = 0.5, size = 14)) 
ggsave("unannotated.png",ggp,height = 2500,width = 2500,units = "px",dpi = 200)

# resolution 0.2-1.6
for(i in c(0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4)){
  object_clean <- Seurat::FindClusters(object_clean, resolution = i)
}

ggp=clustree(object_clean, prefix = "RNA_snn_res.", 
             node_colour = "sc3_stability", use_core_edges = T, 
             node_label = NULL, node_label_aggr = NULL)   +
  scale_edge_color_continuous(low = "blue", high = "red") +
  guides(edge_colour = FALSE, edge_alpha = FALSE) +
  theme(legend.position = "bottom") 
ggsave("clustree.png",ggp,height = 2500,width = 2500,units = "px",dpi = 200)

object_clean=SetIdent(object_clean,value=object_clean@meta.data$RNA_snn_res.1.6)
Idents(object_clean)

ggp=DimPlot(object_clean, reduction="integrated_umap", label=TRUE, label.size=6) + ggtitle("macrophages") + 
  my_theme + theme(plot.title = element_text(hjust = 0.5, size = 14)) 
ggsave("unannotated.png",ggp,height = 2500,width = 2500,units = "px",dpi = 200)

# cell number by clusters
cts=as.data.frame.matrix(table(Idents(sc.data), sc.data$Condition), margin = 2)
cts=cts %>% 
  rownames_to_column("cluster")
cts.m=reshape2::melt(cts)

cts.m$cluster=as.factor(cts.m$cluster)
write.table(cts,file="cellcounts.tsv",sep="\t",row.names = FALSE)

# Create a stacked bar plot
ggp=ggplot(cts.m, aes(x = cluster, y = value, fill = variable)) +
  geom_bar(stat = "identity") + 
  labs(x = "Cell cluster", y = "Count", title = "Stacked Bar Plot of Cell Types") +
  my_theme

cts.sample=as.data.frame.matrix(table(Idents(sc.data), sc.data$sample_group), margin = 2)
cts.sample=cts.sample %>% 
  rownames_to_column("cluster")
write.table(cts.sample,file="cellcounts_persample.tsv",sep="\t",row.names = FALSE)

# marker gene plots
marker_genes(sc.data)

sc.data=JoinLayers(sc.data)
sc.data.markers=FindAllMarkers(sc.data, only.pos=TRUE, logfc.threshold = 0.25, min.pct=0.25)
top20.sc.data <- sc.data.markers %>%
  group_by(cluster) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 20)

top50.sc.data=sc.data.markers %>% group_by(cluster) %>% slice_max(n = 50, order_by = avg_log2FC)

# saving the marker genes table
write.table(top20.sc.data,file="top20_markers.tsv",sep="\t",row.names = FALSE)
write.table(top50.sc.data,file="top50_markers.tsv",sep="\t",row.names = FALSE)


#### B cells ####

Bcells=readRDS("/Users/bozwarag/Documents/Amber Bozward - Mac/CosMx data download/AIH & Seroneg matched/CosMX analysis/10.12.24/Plots made by me/Bcell_edited.RDS")
View(Bcells@meta.data)

DimPlot(Bcells, reduction="integrated_umap",label=TRUE)

# Split by samples 

DimPlot(subset(Bcells, subset = sample_group == "D"), reduction = "integrated_umap", label = TRUE)

# Displaying markers of interst

ggp=FeaturePlot_scCustom(Bcells, features = c("CD63","C1QA"), reduction = "integrated_umap",label=TRUE)
ggp

markers=c("CD63","C1QA", "MX1")
plot=DotPlot_scCustom(Bcells, features= markers) + NoLegend() +
  theme(axis.text.x = element_text(angle = 90), axis.text.y = element_text(size=9)) + xlab('Genes') + ylab('Clusters') +
  coord_flip() + ggtitle("title of the plot") + theme(plot.title = element_text(hjust = 0.5, size = 14))

plot


#getting the clustree
# resolution 0.2-1.6
for(i in c(0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2)){
  Bcells <- Seurat::FindClusters(Bcells, resolution = i)
}

View(Bcells@meta.data)

# clustree
ggp=clustree(Bcells, prefix = "RNA_snn_res.", 
             node_colour = "sc3_stability", use_core_edges = T, 
             node_label = NULL, node_label_aggr = NULL)   +
  scale_edge_color_continuous(low = "blue", high = "red") +
  guides(edge_colour = FALSE, edge_alpha = FALSE) +
  theme(legend.position = "bottom") 
ggsave("clustree.png",ggp,height = 2500,width = 2500,units = "px",dpi = 200)

Bcells=SetIdent(Bcells,value=Bcells@meta.data$RNA_snn_res.0.8)

Bcell.markers=c("CD27", "IGHD", "CD38", "CD9", "MZB1", "XBP1")

plot=DotPlot_scCustom(Bcells, features= Bcell.markers) + NoLegend() +
  theme(axis.text.x = element_text(angle = 90), axis.text.y = element_text(size=9)) + xlab('Genes') + ylab('Clusters') +
  coord_flip() + ggtitle("title of the plot") + theme(plot.title = element_text(hjust = 0.5, size = 14))

# Annotations
new.cluster.ids=c("Naive B cell", "Plasma cell", "Naive B cell", "Memory B cell", "Memory B cell", "Class switching B cell", "Naive B cell", "Naive B cell")

names(new.cluster.ids)=levels(Bcells)
Bcells= RenameIdents(Bcells,new.cluster.ids)

ggp=DimPlot(Bcells, reduction="integrated_umap",label=TRUE,label.size = 5, pt.size = 1.5, repel=TRUE) + ggtitle("B cells") +
  my_theme + theme(plot.title = element_text(hjust = 0.5, size = 14)) + NoLegend()

# number of cells per cell type
cts.sample=as.data.frame.matrix(table(Idents(Bcells), Bcells$sample_group), margin = 2)
cts.sample=cts.sample %>% 
  rownames_to_column("cluster")
write.table(cts.sample,file="Bcell_cellcounts_persample.tsv",sep="\t",row.names = FALSE)

cts.m=reshape2::melt(cts.sample)
cts.m$cluster=as.factor(cts.m$cluster)

# Create a stacked bar plot
ggp=ggplot(cts.m, aes(x = cluster, y = value, fill = variable)) +
  geom_bar(stat = "identity") + 
  labs(x = "Cell cluster", y = "Count", title = "Stacked Bar Plot of Cell Types") +
  my_theme

# Swapping so samples are on the x axis
cts.sample=as.data.frame.matrix(table(Idents(Bcells), Bcells$sample_group), margin = 2)
cts.sample=cts.sample %>% 
  rownames_to_column("cluster")

cts.m <- cts.sample %>% 
  pivot_longer(cols = -cluster, names_to = "sample", values_to = "value")

ggp = ggplot(cts.m, aes(x = sample, y = value, fill = cluster)) +
  geom_bar(stat = "identity", position = "stack") + 
  labs(x = "Sample", y = "Cell Count", title = "Stacked Bar Plot of Cell Types Per Sample", fill = "Cell Type") +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels if needed

# Plotting % total cells - stacked plot 

cts.m <- cts.sample %>%
  pivot_longer(cols = -cluster, names_to = "sample", values_to = "count") %>%
  group_by(sample) %>%
  mutate(percentage = (count / sum(count)) * 100) %>%
  ungroup()

ggp = ggplot(cts.m, aes(x = sample, y = percentage, fill = cluster)) +
  geom_bar(stat = "identity", position = "stack") + 
  labs(x = "Sample", y = "Percentage of Total Cells", 
       title = "Proportion of Cell Types Per Sample", 
       fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels


# Create a stacked bar plot for each sample individually


# number of cells per cell type
cts.sample=as.data.frame.matrix(table(Idents(Bcells), Bcells$sample_name), margin = 2)
cts.sample=cts.sample %>% 
  rownames_to_column("cluster")

cts.m=reshape2::melt(cts.sample)
cts.m$cluster=as.factor(cts.m$cluster)

# Swapping so samples are on the x axis
cts.sample=as.data.frame.matrix(table(Idents(Bcells), Bcells$sample_name), margin = 2)
cts.sample=cts.sample %>% 
  rownames_to_column("cluster")

cts.m <- cts.sample %>% 
  pivot_longer(cols = -cluster, names_to = "sample", values_to = "value")

ggp = ggplot(cts.m, aes(x = sample, y = value, fill = cluster)) +
  geom_bar(stat = "identity", position = "stack") + 
  labs(x = "Sample", y = "Cell Count", title = "Stacked Bar Plot of Cell Types Per Sample", fill = "Cell Type") +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels if needed

# Plotting % total cells - stacked plot 

cts.m <- cts.sample %>%
  pivot_longer(cols = -cluster, names_to = "sample", values_to = "count") %>%
  group_by(sample) %>%
  mutate(percentage = (count / sum(count)) * 100) %>%
  ungroup()

ggp = ggplot(cts.m, aes(x = sample, y = percentage, fill = cluster)) +
  geom_bar(stat = "identity", position = "stack") + 
  labs(x = "Sample", y = "Percentage of Total Cells", 
       title = "Proportion of Cell Types Per Sample", 
       fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels


#Dotplot for markers used for cell subsets

genes.of.interest=c("CD9","CD27", "XBP1", "CD38", "MZB1", "IGHD")

plot=DotPlot_scCustom(Bcells, features= genes.of.interest) +
  theme(axis.text.x = element_text(angle = 90), axis.text.y = element_text(size=9)) + xlab('Genes') + ylab('Clusters') +
  ggtitle("") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) + 
  my_theme

plot

saveRDS(Bcells, file="/Users/bozwarag/Documents/Amber Bozward - Mac/CosMx data download/AIH & Seroneg matched/CosMX analysis/10.12.24/Amber_CosMx_Nov_2024/subsets/Bcell_edited.RDS")


#### T cells ####

Tcells=readRDS("/Users/bozwarag/Documents/Amber Bozward - Mac/CosMx data download/AIH & Seroneg matched/CosMX analysis/10.12.24/Plots made by me/Tcell_edited.RDS")
View(Tcells@meta.data)

DimPlot(Tcells, reduction="integrated_umap",label=TRUE)

ggp=FeaturePlot_scCustom(Tcells, features = c("CD63","C1QA"), reduction = "integrated_umap",label=TRUE)
ggp

markers=c("CD63","C1QA", "MX1")
plot=DotPlot_scCustom(Bcells, features= markers) + NoLegend() +
  theme(axis.text.x = element_text(angle = 90), axis.text.y = element_text(size=9)) + xlab('Genes') + ylab('Clusters') +
  coord_flip() + ggtitle("title of the plot") + theme(plot.title = element_text(hjust = 0.5, size = 14))

plot


for(i in c(0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6)){
  Tcells <- Seurat::FindClusters(Tcells, resolution = i)
}

# clustree
ggp=clustree(Tcells, prefix = "RNA_snn_res.", 
             node_colour = "sc3_stability", use_core_edges = T, 
             node_label = NULL, node_label_aggr = NULL)   +
  scale_edge_color_continuous(low = "blue", high = "red") +
  guides(edge_colour = FALSE, edge_alpha = FALSE) +
  theme(legend.position = "bottom") 
ggsave("clustree.png",ggp,height = 2500,width = 2500,units = "px",dpi = 200)

Tcells=SetIdent(Tcells,value=Tcells@meta.data$RNA_snn_res.0.8)


Tcell.markers=c("CD3E", "CD4", "CD8A", "IL17A", "KLRB1", "FOXP3", "GZMB", "GZMA", "GZMK", "KLF2", "LAG3")


markers.T1=c("CD3D","CD3E","CD4","CD8A","CD8B","SELL","CCR7","LEF1","KLF2","IL7R","IFNG","IL2","TNF","CXCR5","CXCL13",
             "FOXP3","IL2RA","CTLA4","PDCD1","IKZF2","LAG3","IL17A","IL22","IL26","TWIST1","CD70","CCR9","ITGAE","CXCR3","CXCR6","KLRG1","CD38")

markers.T2=c("KLRB1","TRAV1-2","IL23R","CCR6","NKG7","GNLY","GZMK","GZMA","GZMH","GZMB","EGR1","KLRD1","AREG","FGFBP2","NCAM1",
             "FCGR3A","ZNF683","TRDC","TRGC1","TRDV1","TRDV2","TRDV3","TRGV9","MX1","IFITM3","FOS","CD69","MKI67","XCL1","XCL2","OASL","ITGA1",
             "KLRK1")

markers.T3=c("CISH","CD4","EGR1","IL6ST","IL18R1","SLAMF1","CTSB","CD40LG","TNF","IFNG","CXCR3","CXCR6","CCR6","KLRB1","IL4I1",
             "IL23R","SLC4A10","RORA","RORC","CCL20","ZBTB16","TRAV1-2")

top.markers=c("CCL5", "KLRK1", "NKG7")

plzf=c("ZBTB16", "CD3E", "CD4", "CD8A", "KLRB1")

plot=DotPlot_scCustom(Tcells, features= plzf) + NoLegend() +
  theme(axis.text.x = element_text(angle = 90), axis.text.y = element_text(size=9)) + xlab('Genes') + ylab('Clusters') +
  coord_flip() + ggtitle("T cells") + theme(plot.title = element_text(hjust = 0.5, size = 14))

plot

# Annotations
new.cluster.ids=c("CD4 T cell", "MAIT", "Cytotoxic CD8 T cell", "TMIC", "CD4 T cell", "NKG7+ CD8 T cell", 
                  "Treg", "Cytotoxic CD8 T cell", "CD8 T cell", "CD4 T cell", "CD4 T cell")
names(new.cluster.ids)=levels(Tcells)
Tcells= RenameIdents(Tcells,new.cluster.ids)

ggp=DimPlot(Tcells, reduction="integrated_umap",label=TRUE,label.size = 5, pt.size = 1.5, repel=TRUE) + ggtitle("T cells") +
  my_theme + theme(plot.title = element_text(hjust = 0.5, size = 14)) + NoLegend()

# number of cells per cell type
cts.sample=as.data.frame.matrix(table(Idents(Tcells), Tcells$sample_group), margin = 2)
cts.sample=cts.sample %>% 
  rownames_to_column("cluster")
write.table(cts.sample,file="Tcell_cellcounts_persample.tsv",sep="\t",row.names = FALSE)

cts.m=reshape2::melt(cts.sample)
cts.m$cluster=as.factor(cts.m$cluster)

# Create a stacked bar plot
ggp=ggplot(cts.m, aes(x = cluster, y = value, fill = variable)) +
  geom_bar(stat = "identity") + 
  labs(x = "Cell cluster", y = "Count", title = "Stacked Bar Plot of Cell Types") +
  my_theme

# Swapping so samples are on the x axis stacked plot
cts.sample=as.data.frame.matrix(table(Idents(Tcells), Tcells$sample_group), margin = 2)
cts.sample=cts.sample %>% 
  rownames_to_column("cluster")

cts.m <- cts.sample %>% 
  pivot_longer(cols = -cluster, names_to = "sample", values_to = "value")

ggp = ggplot(cts.m, aes(x = sample, y = value, fill = cluster)) +
  geom_bar(stat = "identity", position = "stack") + 
  labs(x = "Sample", y = "Cell Count", title = "Stacked Bar Plot of Cell Types Per Sample", fill = "Cell Type") +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels if needed

# Plotting % total cells - stacked plot 

cts.m <- cts.sample %>%
  pivot_longer(cols = -cluster, names_to = "sample", values_to = "count") %>%
  group_by(sample) %>%
  mutate(percentage = (count / sum(count)) * 100) %>%
  ungroup()

ggp = ggplot(cts.m, aes(x = sample, y = percentage, fill = cluster)) +
  geom_bar(stat = "identity", position = "stack") + 
  labs(x = "Sample", y = "Percentage of Total Cells", 
       title = "Proportion of Cell Types Per Sample", 
       fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels


#Stacked plot for each individual sample

# number of cells per cell type
cts.sample=as.data.frame.matrix(table(Idents(Tcells), Tcells$sample_name), margin = 2)
cts.sample=cts.sample %>% 
  rownames_to_column("cluster")
write.table(cts.sample,file="Tcell_cellcounts_persample.tsv",sep="\t",row.names = FALSE)

cts.m=reshape2::melt(cts.sample)
cts.m$cluster=as.factor(cts.m$cluster)

# Create a stacked bar plot
ggp=ggplot(cts.m, aes(x = cluster, y = value, fill = variable)) +
  geom_bar(stat = "identity") + 
  labs(x = "Cell cluster", y = "Count", title = "Stacked Bar Plot of Cell Types") +
  my_theme

# Swapping so samples are on the x axis stacked plot
cts.sample=as.data.frame.matrix(table(Idents(Tcells), Tcells$sample_name), margin = 2)
cts.sample=cts.sample %>% 
  rownames_to_column("cluster")

cts.m <- cts.sample %>% 
  pivot_longer(cols = -cluster, names_to = "sample", values_to = "value")

ggp = ggplot(cts.m, aes(x = sample, y = value, fill = cluster)) +
  geom_bar(stat = "identity", position = "stack") + 
  labs(x = "Sample", y = "Cell Count", title = "Stacked Bar Plot of Cell Types Per Sample", fill = "Cell Type") +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels if needed

# Plotting % total cells - stacked plot 

cts.m <- cts.sample %>%
  pivot_longer(cols = -cluster, names_to = "sample", values_to = "count") %>%
  group_by(sample) %>%
  mutate(percentage = (count / sum(count)) * 100) %>%
  ungroup()

ggp = ggplot(cts.m, aes(x = sample, y = percentage, fill = cluster)) +
  geom_bar(stat = "identity", position = "stack") + 
  labs(x = "Sample", y = "Percentage of Total Cells", 
       title = "Proportion of Cell Types Per Sample", 
       fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels


#Dotplot for markers used for cell subsets

genes.of.interest=c("CD8A","CD4", "FOXP3", "TIGIT", "ZBTB16", "NKG7", "GZMB", "GZMK", "GZMA", "CCR7")

plot=DotPlot_scCustom(Tcells, features= genes.of.interest) +
  theme(axis.text.x = element_text(angle = 90), axis.text.y = element_text(size=9)) + xlab('Genes') + ylab('Clusters') +
  ggtitle("") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) + 
  my_theme

plot

saveRDS(Tcells, file="/Users/bozwarag/Documents/Amber Bozward - Mac/CosMx data download/AIH & Seroneg matched/CosMX analysis/10.12.24/Amber_CosMx_Nov_2024/subsets/Tcell_edited.RDS")


#### Looking at SN, AIH and Donor UMAPS separately ####

cosmx_obj=readRDS("/Users/bozwarag/Documents/Amber Bozward - Mac/CosMx data download/AIH & Seroneg matched/CosMX analysis/10.12.24/cosmx_obj.RDS")


#cosmx_obj <- NormalizeData(object)   # Normalize counts
#cosmx_obj <- FindVariableFeatures(cosmx_obj)  # Identify variable genes
#cosmx_obj <- ScaleData(cosmx_obj)  # Scale the data

#cosmx_obj <- RunPCA(cosmx_obj, npcs = 20)  # Adjust `npcs` as needed
#print(cosmx_obj[["pca"]])

#head(cosmx_obj@meta.data)

saveRDS(cosmx_obj, file="/Users/bozwarag/Documents/Amber Bozward - Mac/CosMx data download/AIH & Seroneg matched/CosMX analysis/10.12.24/cosmx_obj.RDS")


## For SN ##

cosmx_subset_SN=readRDS("/Users/bozwarag/Documents/Amber Bozward - Mac/CosMx data download/AIH & Seroneg matched/CosMX analysis/10.12.24/cosmx_subset_SN.RDS")

#condition_of_interest <- "SN"
#cosmx_subset_SN <- subset(cosmx_obj, subset = sample_group == condition_of_interest)
#cosmx_subset_SN <- RunUMAP(cosmx_subset_SN, dims = 1:20)  

Idents(cosmx_subset_SN)=cosmx_subset_SN@meta.data$final_cell_types
DimPlot(cosmx_subset_SN, reduction = "umap") + ggtitle(paste("UMAP for", condition_of_interest)) + geom_text_repel()

#saveRDS(cosmx_subset_SN, file="/Users/bozwarag/Documents/Amber Bozward - Mac/CosMx data download/AIH & Seroneg matched/CosMX analysis/10.12.24/cosmx_subset_SN.RDS")


# plotting on the UMAP
ggp=FeaturePlot_scCustom(cosmx_subset_SN, features = c("CD63"), reduction = "umap",label=TRUE)
ggp + scale_color_gradient(low = "lightgrey", high = "red3") 

ggp = FeaturePlot_scCustom(
  cosmx_subset_SN, 
  features = c("CD63"), 
  reduction = "umap", 
  label = TRUE,
  repel = TRUE,
  min.cutoff = "q30",  # Adjust minimum cutoff (e.g., 10th percentile),
  max.cutoff = "q95"
)




## For AIH ##

cosmx_subset_AIH=readRDS("/Users/bozwarag/Documents/Amber Bozward - Mac/CosMx data download/AIH & Seroneg matched/CosMX analysis/10.12.24/cosmx_subset_AIH.RDS")

condition_of_interest <- "AIH"
cosmx_subset_AIH <- subset(cosmx_obj, subset = sample_group == condition_of_interest)
cosmx_subset_AIH <- RunUMAP(cosmx_subset_AIH, dims = 1:20)

Idents(cosmx_subset_AIH)=cosmx_subset_AIH@meta.data$final_cell_types
DimPlot(cosmx_subset_AIH, reduction = "umap") + ggtitle(paste("UMAP for", condition_of_interest))

#saveRDS(cosmx_subset_AIH, file="/Users/bozwarag/Documents/Amber Bozward - Mac/CosMx data download/AIH & Seroneg matched/CosMX analysis/10.12.24/cosmx_subset_AIH.RDS")

# plotting on the UMAP
ggp=FeaturePlot_scCustom(cosmx_subset_AIH, features = c("CD63"), reduction = "umap",label=TRUE)
ggp

ggp + scale_color_gradient(low = "lightgrey", high = "red3") 

ggp = FeaturePlot_scCustom(
  cosmx_subset_AIH, 
  features = c("CD63"), 
  reduction = "umap", 
  label = TRUE,
  repel = TRUE,
  min.cutoff = "q20",  # Adjust minimum cutoff (e.g., 10th percentile),
  max.cutoff = "q95"
)


## For Donor ##

cosmx_subset_D=readRDS("/Users/bozwarag/Documents/Amber Bozward - Mac/CosMx data download/AIH & Seroneg matched/CosMX analysis/10.12.24/cosmx_subset_D.RDS")

#condition_of_interest <- "D"
#cosmx_subset_D <- subset(cosmx_obj, subset = sample_group == condition_of_interest)
#cosmx_subset_D <- RunUMAP(cosmx_subset_D, dims = 1:20)  # Adjust dims if needed

Idents(cosmx_subset_D)=cosmx_subset_D@meta.data$final_cell_types
DimPlot(cosmx_subset_D, reduction = "umap") + ggtitle(paste("UMAP for", condition_of_interest))

saveRDS(cosmx_subset_D, file="/Users/bozwarag/Documents/Amber Bozward - Mac/CosMx data download/AIH & Seroneg matched/CosMX analysis/10.12.24/cosmx_subset_D.RDS")

ggp + scale_color_gradient(low = "lightgrey", high = "red") 

ggp = FeaturePlot_scCustom(
  cosmx_subset_D, 
  features = c("CD63"), 
  reduction = "umap", 
  label = TRUE,
  repel = TRUE,
  min.cutoff = "q20",  # Adjust minimum cutoff (e.g., 10th percentile),
  max.cutoff = "q99"
)

# plotting on the UMAP
ggp=FeaturePlot_scCustom(cosmx_subset_D, features = c("CD63"), reduction = "umap",label=TRUE)
ggp


#### Plotting markers on the separate cohorts ####

FeaturePlot(cosmx_subset_D, features = c("TNFSF14", "TNFRSF14", "CD63", "SPP1", "CD44", "LTBR", "MERTK", "PDGFRA"),label=TRUE, ncol=8)

# dot plot for markers in "markers"

markers=c("TNFSF14", "TNFRSF14", "CD63", "SPP1", "CD44", "LTBR", "MERTK", "PDGFRA")
DotPlot(cosmx_subset_SN, features= markers) + NoLegend() +
  theme(axis.text.x = element_text(angle = 90), axis.text.y = element_text(size=9)) + xlab('Genes') + ylab('Clusters') +
  coord_flip() + ggtitle("title of the plot") + theme(plot.title = element_text(hjust = 0.5, size = 14))

# Plotting for viral markers which are in the top 100 upregulated in SN compred to both AIH and Donor

FeaturePlot(object, features = c("MX1", "OAS1", "OAS2", "OAS3"),label=TRUE, ncol=4)
FeaturePlot(object, features = c("ISG15", "STAT1", "IFIT1", "IFIT3", "IFI6", "IFI44L"),label=TRUE, ncol=3)

#### Cell chat ####

#install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")
#remotes::install_github("sqjin/CellChat")

#install.packages("devtools")
#devtools::install_github("sqjin/CellChat", force = TRUE)

library(CellChat)

CellChat_NP_SN=readRDS("/Users/bozwarag/Documents/Amber Bozward - Mac/CosMx data download/AIH & Seroneg matched/CosMX analysis/10.12.24/Amber_CosMx_Nov_2024/cellchat/CellChat_NP_SN.RDS")
CellChat_NP_AIH=readRDS("/Users/bozwarag/Documents/Amber Bozward - Mac/CosMx data download/AIH & Seroneg matched/CosMX analysis/10.12.24/Amber_CosMx_Nov_2024/cellchat/CellChat_NP_AIH.RDS")
CellChat_NP_D=readRDS("/Users/bozwarag/Documents/Amber Bozward - Mac/CosMx data download/AIH & Seroneg matched/CosMX analysis/10.12.24/Amber_CosMx_Nov_2024/cellchat/CellChat_NP_D.RDS")


gene_of_interest <- "LIGHT"

# Check available slots in your CellChat object
slotNames(CellChat_NP_SN)

lr_database <- CellChat_NP_SN@DB$interaction
head(lr_database)  # Show the first few rows

unique(lr_database$ligand)   # List all ligands
unique(lr_database$receptor) # List all receptors


# Extract the ligand-receptor interaction database
lr_pairs <- CellChat_NP_SN@DB$interaction

# Find interactions where GeneX is a ligand
ligand_pairs <- lr_pairs[grep(gene_of_interest, lr_pairs$ligand), ]

# Find interactions where GeneX is a receptor
receptor_pairs <- lr_pairs[grep(gene_of_interest, lr_pairs$receptor), ]

# Print results
print(ligand_pairs)
print(receptor_pairs)

netVisual_bubble(CellChat_NP_SN, signaling = gene_of_interest)

## Cant find TNFSF14 in the list so troubleshooting ##
#This is because the pathway is called LIGHT

# Check all signaling pathways present in your CellChat object
unique(CellChat_NP_SN@netP$pathways)

## Changing the cell types in the plots - currently not working ##

unique(CellChat_NP_SN@idents)  # List all cell types

netVisual_bubble(CellChat_NP_SN, 
                 signaling = gene_of_interest, 
                 sources.use = "Macrophage")  # Replace with your cell type name

### Create a circular plot for the specific signaling pathway

# Visualize the contribution of cell types

netVisual_chord_cell(CellChat_NP_SN, signaling = "CD45")
netVisual_chord_cell(CellChat_NP_AIH, signaling = "GAS")
netVisual_chord_cell(CellChat_NP_D, signaling = "SPP1")

head(CellChat_NP_SN@net$LRs)

netVisual_chord_gene(CellChat_NP_SN, signaling = "SPP1")

# Visualize pathway activity using heatmaps
netVisual_heatmap(CellChat_NP_SN)

# Visualize pathway activity with specific pathways
netVisual_activity(CellChat_NP_SN, pathway = "LIGHT")

ls("package:CellChat")



#### Plotting a dotplot of markers for each cell type ####

## run this first to set up the plotting theme. We can customise later if needed:

# my_theme
my_theme = theme(
  panel.grid = element_blank(),
  panel.background = element_blank(),
  panel.border = element_rect(colour = "black", fill=NA, size=1),
  plot.background = element_blank(),
  legend.background = element_rect(fill="transparent", colour=NA),
  legend.key = element_rect(fill="transparent", colour=NA),
  plot.title = element_text(size=12, margin = margin(b = 5),hjust=0,vjust=0.5, family="Arial", face="bold"),
  title = element_text(size = 12, margin = margin(b = 5),hjust=0,vjust=0.5, family="Arial", face="bold"),
  axis.text.y = element_text(size = 11, margin = margin(r = 5),hjust=1,vjust=0.5, family="Arial", face="bold",colour="black"),
  axis.text.x = element_text(size = 11, margin = margin(t = 5),hjust=1,vjust=0.5, family="Arial", face="bold",colour="black"),
  axis.title.y = element_text(size = 12, margin = margin(r = 10),angle = 90,hjust=0.5,vjust=0.5, family="Arial", face="bold"),
  axis.title.x = element_text(size = 12, margin = margin(t = 10),hjust=0.5,vjust=1, family="Arial", face="bold"),
  legend.text=element_text(size=12, family="Arial", face="bold"),
  legend.key.size=unit(1,"line"),
  plot.margin=unit(c(0.4,0.4,0.4,0.4), "cm"),
  strip.text.x = element_text(size = 12, family="Arial", face="bold", vjust=1),
  panel.spacing = unit(1, "lines")
)

# set the identity of seurat object to the column that has the annotations if it has not been already set
Idents(object)= object$final_cell_types

## For default scale (-2 to 2)

#install.packages("scCustomize")

library(scCustomize)

# ordering the cell subsets on the Y axis 

object@meta.data$final_cell_types=factor(object@meta.data$final_cell_types, levels=c("Hepatocyte", "IS Hepatocyte", "Hepatic stellate", "Epithelial", "Type 2 LSEC", "Type 1 LSEC", "Neutrophil", "Macrophage","Monocyte", "T cell", "NK cell", "Plasma cell", "B cell"))

Idents(object)=object@meta.data$final_cell_types

genes.of.interest=c("MZB1", "JCHAIN", "IGHA1", "IGHG1", "IGHG2", "IGKC", "CD27", "NCAM1", "NKG7", "CD2", "CD3E", "CD8A", "CD4", "CCL5", "CD69", "CD14", "CD68", "CD163", "LYZ", "MARCO", "ITGAM", "CD33", "FCGR3A/B", "LYVE1", "CD36", "ICAM1", "HSPG2", "FLT1", "LDB2", "CDH1", "EPCAM", "KRT19", "SOX9", "CD24", "ACTA2", "THBS1", "COL1A1", "PDGFRA", "TCF7", "SERPINA1", "APOA1", "ARG1", "APOC1")
genes=c("MX1", "OAS1", "OAS2", "ISG15", "STAT1", "IFIT3", "IFI6", "IFI44L")

plot=DotPlot_scCustom(object, features= genes) +
  theme(axis.text.x = element_text(angle = 90), axis.text.y = element_text(size=9)) + xlab('Genes') + ylab('Clusters') +
  ggtitle("") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.5) + 
  my_theme

plot

## For scale between 0 and 1 ## I did not use this

# Min-max normalization function (scales between 0 and 1)

min_max_scale <- function(x)
{
  (x - min(x)) / (max(x) - min(x))
}

# plot
pct_cells <- DotPlot_scCustom(object, features = genes, group.by = "final_cell_types")$data #feature = the list of genes in use

colnames(pct_cells)=c("avg.exp","pct.cells","Genes","Cluster","Mean.exp")

pct_cells$Mean.exp=min_max_scale(pct_cells$Mean.exp)

ggplot(pct_cells, aes(x = Genes, y = Cluster, size = pct.cells, color = Mean.exp)) +
  geom_point() +
  scale_color_viridis(option = "D", limits = c(0, 1)) +  # Ensure 0 to 1 range
  my_theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_size_continuous(range = c(0, 8)) # range for the size of dots


#### Plotting molecules ####


# get markers
DefaultAssay(object) = 'RNA'
object = SetIdent(object, value = "final_cell_types")
markers_RNA = FindAllMarkers(object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers_RNA %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC) -> top25_RNA

# get pseudo by cell type
pseudo = AggregateExpression(object, group.by = c("final_cell_types"), return.seurat = TRUE)

# Perform log-normalization
pseudo = NormalizeData(pseudo)

# get expression matrix
pseudo_em = data.frame(pseudo@assays$RNA$data)

# filter for markers
pseudo_em = pseudo_em[top25_RNA$gene,]

# get annotations
pseudo_meta = pseudo$final_cell_types

# tidy up
rm(object, pseudo)
gc()


#### Slide 4 ####
## Donor sample
# load
sc.data = readRDS("/Users/bozwarag/Documents/Amber Bozward - Mac/CosMx data download/AIH & Seroneg matched/RDS/Segmentation and molecules/seuratObject_Amber.Slide.5.RDS")

# Extract expression information
exp.ma = sc.data@assays$RNA$data

# Get predictions
predictions_internal = SingleR(test = exp.ma, ref = pseudo_em , labels = pseudo_meta)

# add predictions to metadata
sc.data = AddMetaData(sc.data, predictions_internal$labels, col.name = "predictions_internal")

# set ident
sc.data = SetIdent(sc.data, value = "predictions_internal")

# switch from centroids to segmentation (its in the object, under Images)
DefaultBoundary(sc.data[["Amber.Slide.4"]]) = "segmentation"

# plot molecules onto segmentation
png("Slide_4_segmentation_molecules.png", height = 5000, width = 5000)
print(ImageDimPlot(sc.data, axes = FALSE, size = 2, alpha = 1, mols.size = 1, mols.alpha = 1, coord.fixed = FALSE, molecules = c("CD63"), nmols = 10000000) +
  xlim(9.5,9.8) + ylim(4,4.3))


#### Slide 1 #### - not working as file cannot export with cell segmentation and molecules

# load
sc.data = readRDS("/Users/bozwarag/Desktop/Naomi_slide_1_02_04_2025_15_18_43_657/seuratObject_Naomi.slide.1.RDS")

# Extract expression information
exp.ma = sc.data@assays$RNA$data

# Get predictions
predictions_internal = SingleR(test = exp.ma, ref = pseudo_em , labels = pseudo_meta)

# add predictions to metadata
sc.data = AddMetaData(sc.data, predictions_internal$labels, col.name = "predictions_internal")

# set ident
sc.data = SetIdent(sc.data, value = "predictions_internal")

# switch from centroids to segmentation (its in the object, under Images)
DefaultBoundary(sc.data[["Amber.slide.1"]]) = "segmentation"

# plot molecules onto segmentation
png("Slide_4_segmentation_molecules.png", height = 5000, width = 5000)
print(ImageDimPlot(sc.data, axes = TRUE, size = 2, alpha = 1, mols.size = 1, mols.alpha = 1, coord.fixed = FALSE, molecules = c("CD63"), nmols = 10000000) +
        xlim(5.5,6) + ylim(5.5,6))
dev.off()

view(sc.data)


#### AIH vs SN - molecules images ####
## SN_606
# load
sc.data = readRDS("/Users/bozwarag/Documents/Amber Bozward - Mac/CosMx data download/AIH & Seroneg matched/AIH_SN_Donor_RDS_molecules/SN_606.RDS")

#Normalise
sc.data <- NormalizeData(sc.data)

# Extract expression information
exp.ma = sc.data@assays$RNA$data

# Get predictions
predictions_internal = SingleR(test = exp.ma, ref = pseudo_em , labels = pseudo_meta)

# add predictions to metadata
sc.data = AddMetaData(sc.data, predictions_internal$labels, col.name = "predictions_internal")

# set ident
sc.data = SetIdent(sc.data, value = "predictions_internal")

# switch from centroids to segmentation (its in the object, under Images)
DefaultBoundary(sc.data[["SN_606"]]) = "segmentation"

# plot molecules onto segmentation
png("SN_606_segmentation_molecules.png", height = 5000, width = 5000)
print(ImageDimPlot(sc.data, axes = FALSE, size = 2, alpha = 1, mols.size = 1, mols.alpha = 1, coord.fixed = FALSE, molecules = c("CD63"), nmols = 10000000) +
        xlim(9.5,9.8) + ylim(4,4.3))




#### PBC data for grant application ####

# load
sc.data = readRDS("/Users/bozwarag/Documents/Amber Bozward - Mac/SJT/2025/CosMx/RDS objects/seuratObject_Naomi.Slide.2.RDS")

# Extract expression information
exp.ma = sc.data@assays$RNA$data

# Get predictions
predictions_internal = SingleR(test = exp.ma, ref = pseudo_em , labels = pseudo_meta)

# add predictions to metadata
sc.data = AddMetaData(sc.data, predictions_internal$labels, col.name = "predictions_internal")

# set ident
sc.data = SetIdent(sc.data, value = "predictions_internal")

# switch from centroids to segmentation (its in the object, under Images)
DefaultBoundary(sc.data[["Naomi.Slide.2"]]) = "segmentation"

# plot molecules onto segmentation
png("Slide_4_segmentation_molecules.png", height = 5000, width = 5000)
print(ImageDimPlot(sc.data, axes = TRUE, size = 2, alpha = 1, mols.size = 0.5, mols.alpha = 1, coord.fixed = FALSE, molecules = c("TCF7"), nmols = 10000000) +
        xlim(4,5) + ylim(4,5))


# test FOV's 

test = subset(sc.data, subset=fov %in% c(1,2))

ImageDimPlot(test, fov="Naomi.Slide.2", axes = TRUE, size = 2, alpha = 1, mols.size = 0.5, mols.alpha = 1, coord.fixed = FALSE, molecules = c("TCF7"), nmols = 10000000 )



#### Retrieving the top 100 genes ####
## SN vs AIH&donor

object=readRDS("/Users/bozwarag/Documents/Amber Bozward - Mac/CosMx data download/AIH & Seroneg matched/CosMX analysis/10.12.24/final.RDS")

colnames(object@meta.data)

table(object$sample_group)

# Create a new grouping variable
object$group_compare <- ifelse(object$sample_group == "SN", "SN", "Other")
table(object$group_compare)

# Run differential expression
markers <- FindMarkers(
  object,
  ident.1 = "SN",
  ident.2 = "Other",
  group.by = "group_compare",
  logfc.threshold = 0.25,     # You can adjust this
  min.pct = 0.1,
  test.use = "wilcox"         # Default; or "MAST", "DESeq2", etc.
)

# Filter for upregulated genes only (avg_log2FC > 0), sort by FC
upregulated_genes <- markers[markers$avg_log2FC > 0, ]
top100_up <- head(upregulated_genes[order(-upregulated_genes$avg_log2FC), ], 100)

# Write to CSV
write.csv(top100_up, "top100_upregulated_SN_vs_Other.csv")


## SN vs AIH

object=readRDS("/Users/bozwarag/Documents/Amber Bozward - Mac/CosMx data download/AIH & Seroneg matched/CosMX analysis/10.12.24/final.RDS")


# Subset the Seurat object to only SN and AIH samples
sub_obj <- subset(object, subset = group %in% c("SN", "AIH"))

Idents(object) <- object$sample_group


# Run differential expression
markers_sn_vs_aih <- FindMarkers(
  object,
  ident.1 = "SN",
  ident.2 = "AIH",
  logfc.threshold = 0.25,  # Adjust as needed
  min.pct = 0.1,
  test.use = "wilcox"      # or "MAST" for more robust stats
)

# Filter for upregulated genes (positive log2FC)
up_sn <- markers_sn_vs_aih[markers_sn_vs_aih$avg_log2FC > 0, ]

# Sort and take top 100
top100_sn_vs_aih <- head(up_sn[order(-up_sn$avg_log2FC), ], 100)

# View or export
View(top100_sn_vs_aih)
write.csv(top100_sn_vs_aih, "top100_SN_vs_AIH.csv")


#### Boxplot for markers of interest ####

head(seurat_obj@meta.data)

# If condition is not yet set, you may need to assign it based on sample or group ID:
# seurat_obj$condition <- seurat_obj$sample_group

# Violin plot with box overlay
VlnPlot(
  object = object,
  features = "MX1",
  group.by = "sample_group",
  pt.size = 0.1
) + geom_boxplot(width = 0.1, fill = NA, outlier.shape = NA)

# Multiple genes at once 

# List of genes you're interested in
genes <- c("MX1", "OAS1", "OAS2", "ISG15", "STAT1", "IFIT3", "IFI6", "IFI44L")  # Replace with your genes
genes1 <- c("TNFSF14", "TNFRSF14", "LTBR")

gene_order <- c("MX1", "OAS1", "OAS2", "ISG15", "STAT1", "IFIT3", "IFI6", "IFI44L")

# Fetch expression and condition info
expr_data <- FetchData(object, vars = c(genes1, "sample_group"))

# Reshape data from wide to long format
expr_long <- expr_data %>%
  pivot_longer(cols = all_of(genes1), names_to = "gene", values_to = "expression")

# Make sure 'gene' is a factor in that order
expr_long$gene <- factor(expr_long$gene, levels = gene_order)

# Plot with facets
plot_obj <- ggplot(expr_long, aes(x = sample_group, y = expression, fill = sample_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  facet_wrap(~ gene, scales = "free_y", ncol = 4) +
  theme_minimal(base_size = 12) +
  labs(title = "Gene Expression Across Conditions",
       x = "sample_group", y = "Expression") +
  scale_fill_manual(values = c(
    "SN" = "orange",
    "D" = "blue",
    "AIH" = "red"
  ))

print(plot_obj)

# Step 1: Define your gene of interest
gene <- "TIGIT"

# Step 2: Extract expression values + condition metadata
df <- FetchData(object, vars = c(gene, "sample_group"))

# Step 3: Rename columns for clarity
colnames(df) <- c("expression", "sample_group")

# Step 4: Plot using ggplot2
ggplot(df, aes(x = sample_group, y = expression, fill = sample_group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  theme_minimal(base_size = 14) +
  labs(title = paste("Expression of", gene, "Across Conditions"),
       x = "Condition", y = "Expression Level") +
  scale_fill_brewer(palette = "Set2")



#### Plotting a boxplot of a gene in cell subsets ####

# Use the already made seurat objects for each cell type
# Fetch expression and condition info
expr_data <- FetchData(Macrophages, vars = c(genes, "sample_group"))

# Reshape data from wide to long format
expr_long <- expr_data %>%
  pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "expression")

# Make sure 'gene' is a factor in that order
expr_long$gene <- factor(expr_long$gene, levels = gene_order)

# Plot with facets
plot_obj <- ggplot(expr_long, aes(x = sample_group, y = expression, fill = sample_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  facet_wrap(~ gene, scales = "free_y", ncol = 4) +
  theme_minimal(base_size = 12) +
  labs(title = "Viral marker expression in Macrophages",
       x = "sample_group", y = "Expression") +
  scale_fill_manual(values = c(
    "SN" = "orange",
    "D" = "blue",
    "AIH" = "red"
  ))

print(plot_obj)


#### Plotting genes on the CosMx tissue image ####

# Visualize TNFSF14 expression across the tissue
SpatialFeaturePlot(object, features = "TNFSF14", pt.size.factor = 1.6)


plot_obj <- ggplot(object@meta.data, aes(x = CenterX_global_px, y = CenterY_global_px, color = TNFSF14_expr)) +
  geom_point(size = 0.5) +
  scale_color_viridis_c(na.value = "grey90") +
  scale_y_reverse() +  # Optional: flips Y axis to match tissue orientation
  coord_fixed() +
  theme_void() +
  ggtitle("TNFSF14 Expression Across Tissue")

print(plot_obj)


#### Creating raw count files from SN and AIH samples ####

pseudo =AggregateExpression(object, group.by= 'sample_name', return.seurat=TRUE)
pseudo =data.frame(pseudo@assays$RNA$counts)
write.table(pseudo,file="count_matrix_per_sample.tsv",sep="\t",row.names = TRUE, col.names=NA)


# Final rds object doesnt have AIH_2709 so need to do this separately
AIH_2709 = readRDS("/Users/bozwarag/Documents/Amber Bozward - Mac/CosMx data download/AIH & Seroneg matched/AIH_SN_Donor_RDS_molecules/AIH_2709.RDS")
pseudo =AggregateExpression(AIH_2709, group.by= 'orig.ident', return.seurat=TRUE)
pseudo =data.frame(pseudo@assays$RNA$counts)
write.table(pseudo,file="count_matrix_per_sample_AIH_2709.tsv",sep="\t",row.names = TRUE, col.names=NA)


