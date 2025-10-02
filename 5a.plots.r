#### Load Libraries ####

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(devtools)
library(hdf5r)
library(Biobase)
library(RColorBrewer)
library(viridis)
library(amap)
library(reshape2)
library(stringr)

#### Guides and Manuals ####

## https://bioconductor.org/packages/release/bioc/vignettes/granulator/inst/doc/granulator.html#reference-profiles




#### Helper Functions ####

save_plot = function(ggp,plot_height,plot_width,plot_path)
{
  png(plot_path, height=plot_height, width=plot_width, pointsize=5)
  print(ggp)
  dev.off()
  
  #svg(gsub(".png", ".svg", plot_path), height=plot_height/90, width=plot_width/90)
  #print(ggp)
  #dev.off()
  
  # clears all devices - as a safety measure
  while (dev.cur()>1) dev.off()
}


make_heatmap = function(hm.matrix)
{
  
  # does the y clustering
  y.dist = Dist(hm.matrix, method="spearman")
  y.cluster = hclust(y.dist, method="average")
  y.dd = as.dendrogram(y.cluster)
  y.dd.reorder = reorder(y.dd,0,FUN="average")
  y.order = order.dendrogram(y.dd.reorder)
  
  # does the y clustering
  x.dist = Dist(t(hm.matrix), method="spearman")
  x.cluster = hclust(x.dist, method="average")
  x.dd = as.dendrogram(x.cluster)
  x.dd.reorder = reorder(x.dd,0,FUN="average")
  x.order = order.dendrogram(x.dd.reorder)
  
  hm.matrix_clustered = hm.matrix[y.order,x.order]
  
  # melt and plot
  hm.matrix_clustered = melt(hm.matrix_clustered)
  
  # colour palette
  colours = c("purple","black","yellow")
  palette = colorRampPalette(colours)(100)
  
  # plot
  ggp = ggplot(hm.matrix_clustered, aes(x=Var2, y=Var1, fill=value)) + 
    geom_tile() + 
    scale_fill_gradientn(colours = palette) + 
    labs(x="", y="") + 
    theme(legend.position="right", legend.title = element_blank(), legend.spacing.x = unit(0.25, 'cm'),axis.text.x = element_text(angle = 90, hjust = 1), axis.ticks=element_blank())
    
  return(ggp)
}




#### Load data ####

## wd
setwd("C:/Users/JohnJ/Desktop/analysis/Amber_Bozward_Visium_2023/analysis")

## load
liver_st = readRDS(file = "liver_integrated_light_psudo_decon_linneages.rds")
liver_sc = readRDS(file = "public/GSE158723/single_cell_integrated.rds")

liver_st_cluster_markers = readRDS(file = "2a_integration_light/st_clusters_all_markers.drs")
liver_sc_cluster_markers = read.table(file = "public/GSE158723/all_markers.csv")

decon = readRDS(file = "3d.pseudobulk_deconvolution/pseudo_decon.rds")
liver_st_cluster_pseduo = readRDS(file = "3d.pseudobulk_deconvolution/pseudo_st.rds")
liver_sc_cluster_pseduo = readRDS(file = "3d.pseudobulk_deconvolution/pseudo_sc.rds")

decon_linneages = readRDS(file = "3e.pseudobulk_linneages_deconvolution/pseudo_decon.rds")
liver_st_cluster_pseduo_linneages = readRDS(file = "3e.pseudobulk_linneages_deconvolution/pseudo_st.rds")
liver_sc_cluster_pseduo_linneages = readRDS(file = "3e.pseudobulk_linneages_deconvolution/pseudo_sc.rds")



#### To Do ####

# Differential by Ambers regions
# Differential more generally
# Deconvolute without sub-cell types, just overall linneages


################## 1. MARKERS ################## 



#### Albumin ####

## Markers give nice granularity

# Albumin spatial
plot = SpatialFeaturePlot(liver_st, features = "ALB", stroke = NA, pt.size.factor = 7, ncol=5)
png("5a_plots/Alb_Spatial.png", height = 2000, width = 2000)
print(plot)
dev.off()

# Albumin UMAP
plot = FeaturePlot(liver_st, reduction = "umap", features = "ALB", pt.size = 3) & scale_color_gradientn(colors = plasma(n = 10, direction = -1), limits = c(0, 5))
png("5a_plots/Alb_UMAP.png", height = 1000, width = 1000)
print(plot)
dev.off()



################## 2. CLUSTERS ################## 

#### Cluster UMAP ####

# Plot UMAP
ggp = DimPlot(liver_st, reduction = "umap", group.by = c("ident"), pt.size = 2)
save_plot(ggp, 1000,1000,"5a_plots/clusters_overview/clusters_UMAP.png")


#### Cluster Spatial ####

# Plot spatial clusters
ggp = SpatialDimPlot(liver_st, pt.size.factor = 6, stroke = NA, ncol=6)
save_plot(ggp, 1000, 7500, "5a_plots/clusters_overview/clusters_spatial.png")


#### Cluster Markers Heatmap ####

# get the top 10 for each cluster
liver_st_cluster_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> liver_st_cluster_markers_top10

# plot heamtap
ggp = DoHeatmap(liver_st, features = liver_st_cluster_markers_top10$gene) + NoLegend() + theme(text = element_text(size = 20))
save_plot(ggp, 2000, 2000, "5a_plots/cluster_markers/clusters_heatmap.png")


#### Cluster Markers UMAP ####

liver_st_cluster_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> liver_st_cluster_markers_top5


for (cluster_number in 0:(length(unique(liver_st_cluster_markers_top5$cluster))-1))
{
  markers = subset(liver_st_cluster_markers_top5, cluster == cluster_number)$gene
  ggp = FeaturePlot(liver_st, reduction = "umap", features = markers, ncol=5,pt.size = 2) & scale_color_gradientn(colors = plasma(n = 10, direction = -1), limits = c(0, 8))
  save_plot(ggp, 1000,2500,paste0("5a_plots/cluster_markers/cluster_markers_UMAP_",cluster_number,".png"))
}


#### Cluster Markers Spatial ####

for (cluster_number in 0:(length(unique(liver_st_cluster_markers_top5$cluster))-1))
{
  markers = subset(liver_st_cluster_markers_top5, cluster == cluster_number)$gene
  ggp = SpatialFeaturePlot(liver_st, features = markers, stroke = NA, pt.size.factor = 8, ncol=11)
  save_plot(ggp, 2500,2500,paste0("5a_plots/cluster_markers/cluster_markers_Spatial_",cluster_number,".png"))
}




#### Cluster Frequency by Disease ####

# setup the table
cf_by_disease = data.frame(matrix(0,nrow=length(unique(pt$Var2)),ncol=length(unique(pt$Var1))))
names(cf_by_disease) = unique(pt$Var1)
row.names(cf_by_disease) = unique(pt$Var2)

# fill the table
for (row_index in 1:nrow(pt))
{
  cluster = as.character(pt[row_index,"Var1"])
  sample = as.character(pt[row_index,"Var2"])
  frequency = pt[row_index,"Freq"]
  cf_by_disease[sample,cluster] = frequency
}

# Normalise by sample
cf_by_disease = cf_by_disease / rowSums(cf_by_disease)

# add groups
cf_by_disease$groups = str_split_fixed(row.names(cf_by_disease), "_", 2)[,1]

# p-vals
for (col_index in 1:ncol(cf_by_disease))
{
  AH = cf_by_disease[1:5,col_index]
  SN = cf_by_disease[6:11,col_index]
  p = t.test(AH, SN)$p.value
  print(p)
}

# melt
cf_by_disease = melt(cf_by_disease,id.vars="groups")

# plot
ggp = ggplot(cf_by_disease, aes(x = variable, y = value, fill = groups)) + 
  geom_boxplot(alpha = 0.25,  width = 0.5) + 
  labs(x = "Spatial cluster", y = "Fraction of Spots Per Sample (scaled)")
ggp

save_plot(ggp, 500,750,"5a_plots/clusters_overview/cluster_frequency_by_disease.png")




#### Cluster Frequency by Disease - Scaled ####

# setup the table
cf_by_disease = data.frame(matrix(0,nrow=length(unique(pt$Var2)),ncol=length(unique(pt$Var1))))
names(cf_by_disease) = unique(pt$Var1)
row.names(cf_by_disease) = unique(pt$Var2)

# fill the table
for (row_index in 1:nrow(pt))
{
  cluster = as.character(pt[row_index,"Var1"])
  sample = as.character(pt[row_index,"Var2"])
  frequency = pt[row_index,"Freq"]
  cf_by_disease[sample,cluster] = frequency
}

# Normalise by sample
cf_by_disease = cf_by_disease / rowSums(cf_by_disease)

# Scale
cf_by_disease = data.frame(scale(cf_by_disease))

# add groups
cf_by_disease$groups = str_split_fixed(row.names(cf_by_disease), "_", 2)[,1]

# melt
cf_by_disease = melt(cf_by_disease,id.vars="groups")

# plot
ggp = ggplot(cf_by_disease, aes(x = variable, y = value, fill = groups)) + 
  geom_boxplot(alpha = 0.25,  width = 0.5) + 
  labs(x = "Spatial cluster", y = "Fraction of Spots Per Sample (scaled)")
ggp

save_plot(ggp, 500,750,"5a_plots/clusters_overview/cluster_frequency_by_disease_scaled.png")



#### Cluster Frequency by Sample ####

# makes a nice table of the number of spots in each cluster, by sample and condition
pt = table(Idents(liver_st), liver_st$orig.ident)
pt = as.data.frame(pt)
pt$Var1 = as.character(pt$Var1)
pt$Group = str_split_fixed(pt$Var2, "_", 2)[,1]
pt$Var1 = factor(pt$Var1, levels = c(0:14))
levels(pt$Var1)

# plots per sample frequencies as a stacked bar chart
ggp = ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  theme(legend.title = element_blank(), axis.text.x = element_text(angle=45, hjust=1))

save_plot(ggp, 500,750,"5a_plots/clusters_overview/cluster_frequency_by_sample.png")




################## 3. SINGLE CELL ################## 


#### Single Cell Cluster UMAP ####

# TODO


################## 4. CLUSTER DECONVOLUTION ################## 

#### Cluster Deconvolution Spatial ####

# get the fractions - we use SVR.
fractions = decon$proportions$svr_sig1

ggp = SpatialFeaturePlot(liver_st, features = names(fractions), stroke = NA, pt.size.factor = 6)
save_plot(ggp, 12000,5000,"5a_plots/cell_types_overview/cluster_deconvolution_spatial.png")


#### Cluster Deconvolution UMAP ####

ggp = FeaturePlot(liver_st, reduction = "umap", features = names(fractions), ncol=5)
save_plot(ggp, 1500,1500,"5a_plots/cell_types_overview/cluster_deconvolution_UMAP.png")


#### Cluster Deconvolution Heatmap ####

hm.data = decon$proportions$svr_sig1
hm.data = scale(hm.data)
row.names(hm.data) = paste0("cluster_", row.names(hm.data))
ggp = make_heatmap(hm.data)
save_plot(ggp, 500,500,"5a_plots/cell_types_overview/cluster_deconvolution_heatmap.png")


#### Cluster Deconvolution Barchart - by cluster ####

bc.data = decon$proportions$svr_sig1
bc.data[bc.data<0] <- 0
bc.data$cluster = paste0("cluster_",row.names(bc.data))
bc.data = melt(bc.data, id.vars="cluster")
ggp = ggplot(bc.data, aes(x = value, y = variable)) + geom_bar(stat = "identity") + facet_wrap(~cluster)
save_plot(ggp, 1500,1500,"5a_plots/cell_types_overview/cluster_deconvolution_barchart_by_cluster.png")

bc.data = decon$proportions$svr_sig1
bc.data = data.frame(scale(bc.data))
bc.data$cluster = paste0("cluster_",row.names(bc.data))
bc.data = melt(bc.data, id.vars="cluster")
ggp = ggplot(bc.data, aes(x = value, y = variable)) + geom_bar(stat = "identity") + facet_wrap(~cluster) + labs(x = "scaled abundance", y = "cell cluster")
save_plot(ggp, 1500,1500,"5a_plots/cell_types_overview/cluster_deconvolution_scaled_barchart_by_cluster.png")


#### Cluster Deconvolution Barchart - by cell ####

bc.data = decon$proportions$svr_sig1
bc.data[bc.data<0] <- 0
bc.data$cluster = paste0("cluster_",row.names(bc.data))
bc.data = melt(bc.data, id.vars="cluster")
ggp = ggplot(bc.data, aes(x = value, y = cluster )) + geom_bar(stat = "identity") + facet_wrap(~variable)
save_plot(ggp, 1500,1500,"5a_plots/cell_types_overview/cluster_deconvolution_barchart_by_cell.png")

bc.data = decon$proportions$svr_sig1
bc.data = data.frame(scale(bc.data))
bc.data$cluster = paste0("cluster_",row.names(bc.data))
bc.data = melt(bc.data, id.vars="cluster")
ggp = ggplot(bc.data, aes(x = value, y = cluster )) + geom_bar(stat = "identity") + facet_wrap(~variable)
save_plot(ggp, 1500,1500,"5a_plots/cell_types_overview/cluster_deconvolution_scaled_barchart_by_cell.png")





#### Cell Type Marker Genes UMAP ####


liver_sc_cluster_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> liver_sc_cluster_markers_top5

for (cluster_number in 1:length(unique(liver_sc_cluster_markers_top5$cluster)))
{
  cluster_name = unique(liver_sc_cluster_markers_top5$cluster)[cluster_number]
  markers = subset(liver_sc_cluster_markers_top5, cluster == cluster_name)$gene
  ggp = FeaturePlot(liver_st, reduction = "umap", features = markers, ncol=5,pt.size = 2) & scale_color_gradientn(colors = plasma(n = 10, direction = -1), limits = c(0, 8))
  save_plot(ggp, 1000,3000,paste0("5a_plots/cell_type_markers/cluster_markers_UMAP_",cluster_name,".png"))
}


#### Cell Type Marker Genes Spatial ####


for (cluster_number in 1:length(unique(liver_sc_cluster_markers_top5$cluster)))
{
  cluster_name = unique(liver_sc_cluster_markers_top5$cluster)[cluster_number]
  markers = subset(liver_sc_cluster_markers_top5, cluster == cluster_name)$gene
  ggp = SpatialFeaturePlot(liver_st, features = markers, stroke = NA, pt.size.factor = 8, ncol=11)
  save_plot(ggp, 5000,2500,paste0("5a_plots/cell_type_markers/cluster_markers_Spatial_",cluster_name,".png"))
}





################## 5. CLUSTER DECONVOLUTION - LINNEAGES ################## 


#### Linneages - sc UMAP ####

## manually match
new.cluster.ids = c("Hep", "Hep", "Hep", "Hep", "Hep", "Hep", "Hep", "MP", "MP", "Hep", "Hep", "CC", "EC", "Hep", "LP", "LP", "EC", "MP", "MP", "HSC", "HSC", "LP", "LP", "LP", "MP")

# rename
names(new.cluster.ids) = levels(liver_sc)
liver_sc = RenameIdents(liver_sc, new.cluster.ids)

# check
ggp = DimPlot(liver_sc, reduction = "umap", label = TRUE, pt.size = 2)
save_plot(ggp, 1000,1000,"5a_plots/linneages_overview/sc_linneages_UMAP.png")



#### Linneages Cluster Deconvolution Spatial ####

# get the fractions - we use SVR.
fractions = decon_linneages$proportions$svr_sig1

ggp = SpatialFeaturePlot(liver_st, features = names(fractions), stroke = NA, pt.size.factor = 6)
save_plot(ggp, 4000,5000,"5a_plots/linneages_overview/cluster_deconvolution_spatial.png")


#### Linneages Cluster Deconvolution UMAP ####

ggp = FeaturePlot(liver_st, reduction = "umap", features = names(fractions), ncol=3)
save_plot(ggp, 1500,2000,"5a_plots/linneages_overview/cluster_deconvolution_UMAP.png")


#### Linneages Cluster Deconvolution Heatmap ####

hm.data = decon_linneages$proportions$svr_sig1
hm.data = scale(hm.data)
row.names(hm.data) = paste0("cluster_", row.names(hm.data))
ggp = make_heatmap(hm.data)
save_plot(ggp, 500,500,"5a_plots/linneages_overview/cluster_deconvolution_heatmap.png")



#### Linneages Cluster Deconvolution Barchart - by cluster ####

bc.data = decon_linneages$proportions$svr_sig1
bc.data$cluster = paste0("cluster_",row.names(bc.data))
bc.data = melt(bc.data, id.vars="cluster")
ggp = ggplot(bc.data, aes(x = value, y = variable)) + geom_bar(stat = "identity") + facet_wrap(~cluster)
save_plot(ggp, 1500,1500,"5a_plots/linneages_overview/cluster_deconvolution_barchart_by_cluster.png")

bc.data = decon_linneages$proportions$svr_sig1
bc.data = data.frame(scale(bc.data))
bc.data$cluster = paste0("cluster_",row.names(bc.data))
bc.data = melt(bc.data, id.vars="cluster")
ggp = ggplot(bc.data, aes(x = value, y = variable)) + geom_bar(stat = "identity") + facet_wrap(~cluster) + labs(x = "scaled abundance", y = "cell cluster")
save_plot(ggp, 1500,1500,"5a_plots/linneages_overview/cluster_deconvolution_scaled_barchart_by_cluster.png")


#### Linneages Cluster Deconvolution Barchart - by cell ####

bc.data = decon_linneages$proportions$svr_sig1
bc.data$cluster = paste0("cluster_",row.names(bc.data))
bc.data = melt(bc.data, id.vars="cluster")
ggp = ggplot(bc.data, aes(x = value, y = cluster )) + geom_bar(stat = "identity") + facet_wrap(~variable)
save_plot(ggp, 1500,1500,"5a_plots/linneages_overview/cluster_deconvolution_barchart_by_cell.png")

bc.data = decon_linneages$proportions$svr_sig1
bc.data = data.frame(scale(bc.data))
bc.data$cluster = paste0("cluster_",row.names(bc.data))
bc.data = melt(bc.data, id.vars="cluster")
ggp = ggplot(bc.data, aes(x = value, y = cluster )) + geom_bar(stat = "identity") + facet_wrap(~variable)
save_plot(ggp, 1500,1500,"5a_plots/linneages_overview/cluster_deconvolution_scaled_barchart_by_cell.png")



#### Linneages Cell Type Marker Genes UMAP ####

liver_sc = PrepSCTFindMarkers(liver_sc)
liver_sc_cluster_markers_linneages = FindAllMarkers(liver_sc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
liver_sc_cluster_markers_linneages %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> liver_sc_cluster_markers_linneages_top5

for (cluster_number in 1:length(unique(liver_sc_cluster_markers_linneages_top5$cluster)))
{
  cluster_name = unique(liver_sc_cluster_markers_linneages_top5$cluster)[cluster_number]
  markers = subset(liver_sc_cluster_markers_linneages_top5, cluster == cluster_name)$gene
  ggp = FeaturePlot(liver_st, reduction = "umap", features = markers, ncol=5,pt.size = 2) & scale_color_gradientn(colors = plasma(n = 10, direction = -1), limits = c(0, 10))
  save_plot(ggp, 1000,3000,paste0("5a_plots/linneages_markers/cluster_markers_UMAP_",cluster_name,".png"))
}


#### Linneages Cell Type Marker Genes Spatial ####


for (cluster_number in 1:length(unique(liver_sc_cluster_markers_linneages_top5$cluster)))
{
  cluster_name = unique(liver_sc_cluster_markers_linneages_top5$cluster)[cluster_number]
  markers = subset(liver_sc_cluster_markers_linneages_top5, cluster == cluster_name)$gene
  ggp = SpatialFeaturePlot(liver_st, features = markers, stroke = NA, pt.size.factor = 8, ncol=11)
  save_plot(ggp, 5000,2500,paste0("5a_plots/linneages_markers/cluster_markers_Spatial_",cluster_name,".png"))
}


