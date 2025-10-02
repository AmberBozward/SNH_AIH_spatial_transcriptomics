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
de_folder = "analysis/combined_all/immune/"

# load sc object
sc.data = readRDS("analysis/combined_all/all_integrated_named.RDS")


#### Subset for immune ####

# get immune cells
sc.data.sub = subset(sc.data, subset = named_clusters_cell_type == "Macrophage" | named_clusters_cell_type == "T cell" | named_clusters_cell_type == "B cell")

# Perform log-normalization
sc.data.sub = NormalizeData(sc.data.sub)
sc.data.sub = ScaleData(sc.data.sub)

# Run the PCA 
sc.data.sub = RunPCA(object = sc.data.sub, assay ="RNA", reduction.name='pca_unint')
max_dims = 15

# Integrate - using harmony
sc.data.sub = RunHarmony(sc.data.sub,group.by.vars = c("sample_name"),theta = c(2),reduction = "pca_unint", assay.use = "RNA", reduction.save = "integrated_pca")

# re-cluster after integration
sc.data.sub = FindNeighbors(sc.data.sub, reduction = "integrated_pca", dims = 1:max_dims)
sc.data.sub = FindClusters(sc.data.sub, resolution = 0.75, verbose =TRUE, cluster.name = "integrated_clusters" )
sc.data.sub = RunUMAP(sc.data.sub, spread = 6, min.dist = 0.25, dims = 1:max_dims, reduction = "integrated_pca", reduction.name ="integrated_umap")

# Integrated UMAP
p1 = DimPlot(sc.data.sub, reduction = "integrated_umap", group.by = "sample_name", pt.size = 1)
p2 = DimPlot(sc.data.sub, reduction = "integrated_umap", group.by = "sample_group", pt.size = 1)
p3 = DimPlot(sc.data.sub, reduction = "integrated_umap", group.by = "integrated_clusters", pt.size = 1, label = TRUE)
p4 = DimPlot(sc.data.sub, reduction = "integrated_umap", group.by = "named_clusters_cell_type", pt.size = 1, label = TRUE)

ggp = p1 + p2 + p3 + p4
png(paste0(de_folder, "/Immune_UMAP.png"), height = 4000, width = 4000)
print(ggp)
dev.off()


#### Immune markers ####

# set to RNA
DefaultAssay(sc.data.sub) = 'RNA'

# get markers
sc.data.sub = JoinLayers(sc.data.sub)
markers_RNA = FindAllMarkers(sc.data.sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers_RNA %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC) -> top15_RNA

# save markers
write.table(markers_RNA, paste0(de_folder, "/Non_immune_cluster_markers.csv"), quote=FALSE, sep = "\t")
write.table(top15_RNA, paste0(de_folder, "/Non_immune_cluster_markers_top15.csv"), quote=FALSE, sep = "\t")

# Plot clusters heatmap
png(paste0(de_folder, "/Non_immune_cluster_markers_heatmap.png"), height = 4000, width = 2000)
DoHeatmap(sc.data.sub, features = top15_RNA$gene) + NoLegend() + theme(text = element_text(size = 20))
dev.off()



#### Immune SingleR ####

# Download avaialable human primary cell atlas markers
hpca.se = celldex::HumanPrimaryCellAtlasData()

# Extract expression information
exp.ma = sc.data.sub@assays$RNA$data

# Get predictions
predictions_broad = SingleR(test = exp.ma, ref = hpca.se,labels = hpca.se$label.main)

# add predictions to metadata
sc.data.sub = AddMetaData(sc.data.sub, predictions_broad$labels, col.name = "predictions_broad")

# UMAP
p1 = DimPlot(sc.data.sub, reduction = "integrated_umap", group.by = "sample_name", pt.size = 1)
p2 = DimPlot(sc.data.sub, reduction = "integrated_umap", group.by = "sample_group", pt.size = 1)
p3 = DimPlot(sc.data.sub, reduction = "integrated_umap", group.by = "predictions_broad", pt.size = 1, label = TRUE)
p4 = DimPlot(sc.data.sub, reduction = "integrated_umap", group.by = "named_clusters_cell_type", pt.size = 1, label = TRUE)

ggp = p1 + p2 + p3 + p4
png(paste0(de_folder, "/Immune_UMAP_singleR.png"), height = 4000, width = 4000)
print(ggp)
dev.off()

# get markers
sc.data.sub = JoinLayers(sc.data.sub)
markers_RNA = FindAllMarkers(sc.data.sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers_RNA %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC) -> top15_RNA

# Plot clusters heatmap
sc.data.sub = SetIdent(sc.data.sub, value = "predictions_broad")
png(paste0(de_folder, "/Immune_cluster_markers_heatmap_singleR.png"), height = 4000, width = 2000)
DoHeatmap(sc.data.sub, features = top15_RNA$gene) + NoLegend() + theme(text = element_text(size = 20))
dev.off()


# cell types vector
unique(sc.data.sub$predictions_broad)

# plot cell
type = "NK_cell"
cells = Cells(subset(sc.data.sub, subset = predictions_broad == type))
png(paste0(de_folder, "/UMAP_",type,".png"), height = 4000, width = 4000)
DimPlot(sc.data.sub, reduction = "integrated_umap", cells.highlight = cells, cols.highlight = "darkblue", cols= "grey")
dev.off()

# plot cell
type = "Monocyte"
cells = Cells(subset(sc.data.sub, subset = predictions_broad == type))
png(paste0(de_folder, "/UMAP_",type,".png"), height = 4000, width = 4000)
DimPlot(sc.data.sub, reduction = "integrated_umap", cells.highlight = cells, cols.highlight = "darkblue", cols= "grey")
dev.off()

# plot cell
type = "Macrophage"
cells = Cells(subset(sc.data.sub, subset = predictions_broad == type))
png(paste0(de_folder, "/UMAP_",type,".png"), height = 4000, width = 4000)
DimPlot(sc.data.sub, reduction = "integrated_umap", cells.highlight = cells, cols.highlight = "darkblue", cols= "grey")
dev.off()

# plot cell
type = "DC"
cells = Cells(subset(sc.data.sub, subset = predictions_broad == type))
png(paste0(de_folder, "/UMAP_",type,".png"), height = 4000, width = 4000)
DimPlot(sc.data.sub, reduction = "integrated_umap", cells.highlight = cells, cols.highlight = "darkblue", cols= "grey")
dev.off()

# plot cell
type = "Neutrophils"
cells = Cells(subset(sc.data.sub, subset = predictions_broad == type))
png(paste0(de_folder, "/UMAP_",type,".png"), height = 4000, width = 4000)
DimPlot(sc.data.sub, reduction = "integrated_umap", cells.highlight = cells, cols.highlight = "darkblue", cols= "grey")
dev.off()

# plot cell
type = "T_cells"
cells = Cells(subset(sc.data.sub, subset = predictions_broad == type))
png(paste0(de_folder, "/UMAP_",type,".png"), height = 4000, width = 4000)
DimPlot(sc.data.sub, reduction = "integrated_umap", cells.highlight = cells, cols.highlight = "darkblue", cols= "grey")
dev.off()

# plot cell
type = "B_cell"
cells = Cells(subset(sc.data.sub, subset = predictions_broad == type))
png(paste0(de_folder, "/UMAP_",type,".png"), height = 4000, width = 4000)
DimPlot(sc.data.sub, reduction = "integrated_umap", cells.highlight = cells, cols.highlight = "darkblue", cols= "grey")
dev.off()


#### Save immune RDS ####

#saveRDS(sc.data.sub, file="analysis/combined_all/immune_integrated.RDS")
sc.data.sub = readRDS(file="analysis/combined_all/immune_integrated.RDS")


#### Assign immune cluster names ####


# active ident
sc.data.sub = SetIdent(sc.data.sub, value = "seurat_clusters")

# Rename idents
new.cluster.ids = c("Macrophage",
                    "T cell",
                    "Hepatocyte",
                    "Monocyte",
                    "Macrophage",
                    "T cell",
                    "NK cell",
                    "Plasma cell",
                    "Macrophage",
                    "Hepatocyte",
                    "Macrophage",
                    "Hepatocyte",
                    "B cell",
                    "Neutrophil")

names(new.cluster.ids) = levels(sc.data.sub)
sc.data.sub = RenameIdents(sc.data.sub, new.cluster.ids)

# get table of immune cell names
names_immune = data.frame(sc.data.sub@active.ident)
names_immune$name = row.names(names_immune)
names(names_immune) = c("cell_type","cell_id")

# save table
write.table(names_immune, paste0(de_folder, "/Immune_cell_names.csv"), quote=FALSE, sep = "\t")



#### Subset for non-immune ####

# get immune cells
sc.data.sub = subset(sc.data, subset = named_clusters_cell_type != "Macrophage" & named_clusters_cell_type != "T cell" & named_clusters_cell_type != "B cell")

# Perform log-normalization
sc.data.sub = NormalizeData(sc.data.sub)
sc.data.sub = ScaleData(sc.data.sub)

# Run the PCA 
sc.data.sub = RunPCA(object = sc.data.sub, assay ="RNA", reduction.name='pca_unint')
max_dims = 15

# Integrate - using harmony
sc.data.sub = RunHarmony(sc.data.sub,group.by.vars = c("sample_name"),theta = c(2),reduction = "pca_unint", assay.use = "RNA", reduction.save = "integrated_pca")

# re-cluster after integration
sc.data.sub = FindNeighbors(sc.data.sub, reduction = "integrated_pca", dims = 1:max_dims)
sc.data.sub = FindClusters(sc.data.sub, resolution = 0.75, verbose =TRUE, cluster.name = "integrated_clusters" )
sc.data.sub = RunUMAP(sc.data.sub, spread = 6, min.dist = 0.25, dims = 1:max_dims, reduction = "integrated_pca", reduction.name ="integrated_umap")

# Integrated UMAP
p1 = DimPlot(sc.data.sub, reduction = "integrated_umap", group.by = "sample_name", pt.size = 1)
p2 = DimPlot(sc.data.sub, reduction = "integrated_umap", group.by = "sample_group", pt.size = 1)
p3 = DimPlot(sc.data.sub, reduction = "integrated_umap", group.by = "integrated_clusters", pt.size = 1, label = TRUE)
p4 = DimPlot(sc.data.sub, reduction = "integrated_umap", group.by = "named_clusters_cell_type", pt.size = 1, label = TRUE)

ggp = p1 + p2 + p3 + p4
png(paste0(de_folder, "/Non_immune_UMAP.png"), height = 4000, width = 4000)
print(ggp)
dev.off()


#### Non immune markers ####

# set to RNA
DefaultAssay(sc.data.sub) = 'RNA'

# get markers
sc.data.sub = JoinLayers(sc.data.sub)
markers_RNA = FindAllMarkers(sc.data.sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers_RNA %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC) -> top15_RNA

# save markers
write.table(markers_RNA, paste0(de_folder, "/Non_immune_cluster_markers.csv"), quote=FALSE, sep = "\t")
write.table(top15_RNA, paste0(de_folder, "/Non_immune_cluster_markers_top15.csv"), quote=FALSE, sep = "\t")

# Plot clusters heatmap
png(paste0(de_folder, "/Non_immune_cluster_markers_heatmap.png"), height = 4000, width = 2000)
DoHeatmap(sc.data.sub, features = top15_RNA$gene) + NoLegend() + theme(text = element_text(size = 20))
dev.off()



#### Non immune SingleR ####

# Download avaialable human primary cell atlas markers
hpca.se = celldex::HumanPrimaryCellAtlasData()

# Extract expression information
exp.ma = sc.data.sub@assays$RNA$data

# Get predictions
predictions_broad = SingleR(test = exp.ma, ref = hpca.se,labels = hpca.se$label.main)

# add predictions to metadata
sc.data.sub = AddMetaData(sc.data.sub, predictions_broad$labels, col.name = "predictions_broad")

# UMAP
p1 = DimPlot(sc.data.sub, reduction = "integrated_umap", group.by = "sample_name", pt.size = 1)
p2 = DimPlot(sc.data.sub, reduction = "integrated_umap", group.by = "sample_group", pt.size = 1)
p3 = DimPlot(sc.data.sub, reduction = "integrated_umap", group.by = "predictions_broad", pt.size = 1, label = TRUE)
p4 = DimPlot(sc.data.sub, reduction = "integrated_umap", group.by = "named_clusters_cell_type", pt.size = 1, label = TRUE)

ggp = p1 + p2 + p3 + p4
png(paste0(de_folder, "/Non_immune_UMAP_singleR.png"), height = 4000, width = 4000)
print(ggp)
dev.off()

# get markers
sc.data.sub = JoinLayers(sc.data.sub)
markers_RNA = FindAllMarkers(sc.data.sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers_RNA %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC) -> top15_RNA

# Plot clusters heatmap
sc.data.sub = SetIdent(sc.data.sub, value = "predictions_broad")
png(paste0(de_folder, "/Non_immune_cluster_markers_heatmap_singleR.png"), height = 4000, width = 2000)
DoHeatmap(sc.data.sub, features = top15_RNA$gene) + NoLegend() + theme(text = element_text(size = 20))
dev.off()


#### Save non immune RDS ####

#saveRDS(sc.data.sub, file="analysis/combined_all/non_immune_integrated.RDS")
sc.data.sub = readRDS(file="analysis/combined_all/non_immune_integrated.RDS")

#### Assign non immune cluster names ####


# active ident
sc.data.sub = SetIdent(sc.data.sub, value = "seurat_clusters")

# Rename idents

new.cluster.ids = c("Hepatocyte",
                    "Hepatic stellate",
                    "Hepatocyte",
                    "IS Hepatocyte",
                    "Hepatocyte",
                    "Type 2 LSEC",
                    "Hepatocyte",
                    "Hepatocyte",
                    "Epithelial",
                    "Hepatic stellate",
                    "Type 1 LSEC",
                    "Hepatocyte",
                    "Hepatic stellate",
                    "Epithelial",
                    "Hepatocyte",
                    "IS Hepatocyte")

names(new.cluster.ids) = levels(sc.data.sub)
sc.data.sub = RenameIdents(sc.data.sub, new.cluster.ids)

# get table of non immune cell names
names_non_immune = data.frame(sc.data.sub@active.ident)
names_non_immune$name = row.names(names_non_immune)
names(names_non_immune) = c("cell_type","cell_id")

# save table
write.table(names_non_immune, paste0(de_folder, "/Non_immune_cell_names.csv"), quote=FALSE, sep = "\t")




#### Update final names ####

# bind immune and non immune
names_all = rbind(names_non_immune, names_immune)
final_names = names_all$cell_type
names(final_names) = names_all$cell_id

# add to meta data
sc.data = AddMetaData(object = sc.data, metadata = final_names, col.name = 'final_cell_types')

# UMAP
png(paste0(de_folder, "/UMAP_FINAL_NAMES.png"), height = 4000, width = 4000)
DimPlot(sc.data, reduction = "integrated_umap", group.by = "final_cell_types", pt.size = 1, label = TRUE)
dev.off()



#### Colours ####

gg_color_hue = function(n)
{
  # default colour blind palettes
  cblind_palette = c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
  
  # append ggplot colours if groups goes above colour blind palette length
  if (n <= length(cblind_palette))
  {
    return(cblind_palette[1:n])
  }
  else
  {
    nn = n - length(cblind_palette)
    hues = seq(15, 375, length = nn + 1)
    return(append(cblind_palette,hcl(h = hues, l = 65, c = 100)[1:nn]))
  }
}

cell_colours = gg_color_hue(13)




#### Spatial ####

# active ident
sc.data = SetIdent(sc.data, value = "final_cell_types")


# plots
sc.data.sub = subset(sc.data, orig.ident == "AIH_1607")
png(paste0(de_folder, "/AIH_1607_segmentation.png"), width = 4000, height = 4000)
ImageDimPlot(sc.data.sub) + scale_fill_manual(values = cell_colours)
dev.off()

sc.data.sub = subset(sc.data, orig.ident == "AIH_6079")
png(paste0(de_folder, "/AIH_6079_segmentation.png"), width = 4000, height = 4000)
ImageDimPlot(sc.data.sub) + scale_fill_manual(values = cell_colours)
dev.off()

sc.data.sub = subset(sc.data, orig.ident == "SN_606")
png(paste0(de_folder, "/SN_606_segmentation.png"), width = 4000, height = 4000)
ImageDimPlot(sc.data.sub) + scale_fill_manual(values = cell_colours)
dev.off()

sc.data.sub = subset(sc.data, orig.ident == "SN_1042")
png(paste0(de_folder, "/SN_1042_segmentation.png"), width = 4000, height = 4000)
ImageDimPlot(sc.data.sub) + scale_fill_manual(values = cell_colours)
dev.off()

sc.data.sub = subset(sc.data, orig.ident == "SN_2739")
png(paste0(de_folder, "/SN_2739_segmentation.png"), width = 4000, height = 4000)
ImageDimPlot(sc.data.sub) + scale_fill_manual(values = cell_colours)
dev.off()

sc.data.sub = subset(sc.data, orig.ident == "D_7678")
png(paste0(de_folder, "/D_7678_segmentation.png"), width = 4000, height = 4000)
ImageDimPlot(sc.data.sub) + scale_fill_manual(values = cell_colours)
dev.off()

sc.data.sub = subset(sc.data, orig.ident == "D_6414")
png(paste0(de_folder, "/D_6414_segmentation.png"), width = 4000, height = 4000)
ImageDimPlot(sc.data.sub) + scale_fill_manual(values = cell_colours)
dev.off()

sc.data.sub = subset(sc.data, orig.ident == "D_6446")
png(paste0(de_folder, "/D_6446_segmentation.png"), width = 4000, height = 4000)
ImageDimPlot(sc.data.sub) + scale_fill_manual(values = cell_colours)
dev.off()




#### Save final RDS ####

saveRDS(sc.data, file="analysis/combined_all/final.RDS")



