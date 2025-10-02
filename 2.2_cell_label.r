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
library(dittoSeq)


#### Functions ####

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

cell_colours = gg_color_hue(8)





#### Parameters ####


# set memory available
options(future.globals.maxSize = 100000 * 1024^2)


#### Load data ####


# folders
setwd("C:/Users/JohnJ/Desktop/analysis/Ye_Oo/Amber_B_CosMx_2024/")
markers_folder = "analysis/combined_all/markers/"

# load sc object
sc.data = readRDS("analysis/combined_all/all_integrated.RDS")


#### Assign cluster names ####


# active ident
sc.data = SetIdent(sc.data, value = "seurat_clusters")

# Rename idents
new.cluster.ids = c("Hepatocyte",
                    "T cell",
                    "Hepatocyte",
                    "Hepatic stellate 2",
                    "Hepatocyte",
                    "Hepatocyte",
                    "Macrophage",
                    "Endothelial",
                    "Macrophage",
                    "Epithelial",
                    "Hepatic stellate 1",
                    "Hepatic stellate 2",
                    "Macrophage",
                    "Endothelial",
                    "B cell",
                    "Hepatic stellate 2",
                    "Hepatocyte")

names(new.cluster.ids) = levels(sc.data)
sc.data = RenameIdents(sc.data, new.cluster.ids)

# add to meta data
sc.data = AddMetaData(object = sc.data, metadata = sc.data@active.ident, col.name = 'named_clusters_cell_type')



# UMAP
p1 = DimPlot(sc.data, reduction = "integrated_umap", group.by = "sample_name", pt.size = 1)
p2 = DimPlot(sc.data, reduction = "integrated_umap", group.by = "sample_group", pt.size = 1)
p3 = DimPlot(sc.data, reduction = "integrated_umap", group.by = "integrated_clusters", pt.size = 1, label = TRUE)
p4 = DimPlot(sc.data, reduction = "integrated_umap", group.by = "named_clusters_cell_type", pt.size = 1, label = TRUE)

ggp = p1 + p2 + p3 + p4
png(paste0(markers_folder, "/RNA_named_UMAP.png"), height = 4000, width = 4000)
print(ggp)
dev.off()


#### Spatial ####

# plots
sc.data.sub = subset(sc.data, orig.ident == "AIH_1607")
png(paste0(markers_folder, "/AIH_1607_segmentation.png"), width = 4000, height = 4000)
ImageDimPlot(sc.data.sub) + scale_fill_manual(values = cell_colours)
dev.off()

sc.data.sub = subset(sc.data, orig.ident == "AIH_6079")
png(paste0(markers_folder, "/AIH_6079_segmentation.png"), width = 4000, height = 4000)
ImageDimPlot(sc.data.sub) + scale_fill_manual(values = cell_colours)
dev.off()

sc.data.sub = subset(sc.data, orig.ident == "SN_606")
png(paste0(markers_folder, "/SN_606_segmentation.png"), width = 4000, height = 4000)
ImageDimPlot(sc.data.sub) + scale_fill_manual(values = cell_colours)
dev.off()

sc.data.sub = subset(sc.data, orig.ident == "SN_1042")
png(paste0(markers_folder, "/SN_1042_segmentation.png"), width = 4000, height = 4000)
ImageDimPlot(sc.data.sub) + scale_fill_manual(values = cell_colours)
dev.off()

sc.data.sub = subset(sc.data, orig.ident == "SN_2739")
png(paste0(markers_folder, "/SN_2739_segmentation.png"), width = 4000, height = 4000)
ImageDimPlot(sc.data.sub) + scale_fill_manual(values = cell_colours)
dev.off()

sc.data.sub = subset(sc.data, orig.ident == "D_7678")
png(paste0(markers_folder, "/D_7678_segmentation.png"), width = 4000, height = 4000)
ImageDimPlot(sc.data.sub) + scale_fill_manual(values = cell_colours)
dev.off()

sc.data.sub = subset(sc.data, orig.ident == "D_6414")
png(paste0(markers_folder, "/D_6414_segmentation.png"), width = 4000, height = 4000)
ImageDimPlot(sc.data.sub) + scale_fill_manual(values = cell_colours)
dev.off()

sc.data.sub = subset(sc.data, orig.ident == "D_6446")
png(paste0(markers_folder, "/D_6446_segmentation.png"), width = 4000, height = 4000)
ImageDimPlot(sc.data.sub) + scale_fill_manual(values = cell_colours)
dev.off()


#### Frequencies ####

# frequencies 
freqs = table(sc.data@active.ident, sc.data@meta.data$orig.ident)
write.table(freqs, paste0(markers_folder, "/cell_frequencies.csv"), quote=FALSE, sep = "\t")

# plot
png(paste0(markers_folder, "/cell_frequency.png"), height = 500, width = 500)
dittoBarPlot(object = sc.data,var = "named_clusters_cell_type",group.by = "orig.ident")
dev.off()




#### Save RDS ####

saveRDS(sc.data, file="analysis/combined_all/all_integrated_named.RDS")







