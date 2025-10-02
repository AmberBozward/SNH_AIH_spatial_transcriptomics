library(Seurat)
library(spatstat)
library(grid)
library(ggplotify)
library(sf)
library(ggplot2)
library(magick)


#### Parameters ####


# set memory available
options(future.globals.maxSize = 100000 * 1024^2)


#### Load data ####


# folders
setwd("C:/Users/JohnJ/Desktop/analysis/Ye_Oo/Amber_B_CosMx_2024/")

# polygons folder
poly_folder = "analysis/combined_all/polygons/"

# load sc object
sc.data = readRDS("analysis/combined_all/final.RDS")


#### Get data ####

# get points
xy_cells = data.frame(cbind(sc.data$CenterX_global_px,
                            sc.data$CenterY_global_px,
                            sc.data$final_cell_types, 
                            sc.data$fov))

xy_cells$sample = sc.data$orig.ident
names(xy_cells) = c("X","Y","CELL","FOV","SAMPLE")

# save 
write.table(xy_cells, "analysis/combined_all/polygons/centroids.csv", quote=FALSE, sep = "\t", col.names = TRUE)


#### Load cells ####


# load file
selection = read.table("analysis/combined_all/polygons/polygons_out.csv", sep = "\t", header = TRUE, row.names = 1)

# update meta-data
sc.data = AddMetaData(sc.data, selection, col.name = "P_VS_NP")


# active ident
sc.data = SetIdent(sc.data, value = "P_VS_NP")


# save final RDS
saveRDS(sc.data, file="analysis/combined_all/final.RDS")




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





#### plots ####



# plots
sc.data.sub = subset(sc.data, orig.ident == "AIH_1607")
png(paste0(poly_folder, "/AIH_1607_segmentation.png"), width = 4000, height = 4000)
ImageDimPlot(sc.data.sub) + scale_fill_manual(values = cell_colours)
dev.off()

sc.data.sub = subset(sc.data, orig.ident == "AIH_6079")
png(paste0(poly_folder, "/AIH_6079_segmentation.png"), width = 4000, height = 4000)
ImageDimPlot(sc.data.sub) + scale_fill_manual(values = cell_colours)
dev.off()

sc.data.sub = subset(sc.data, orig.ident == "SN_606")
png(paste0(poly_folder, "/SN_606_segmentation.png"), width = 4000, height = 4000)
ImageDimPlot(sc.data.sub) + scale_fill_manual(values = cell_colours)
dev.off()

sc.data.sub = subset(sc.data, orig.ident == "SN_1042")
png(paste0(poly_folder, "/SN_1042_segmentation.png"), width = 4000, height = 4000)
ImageDimPlot(sc.data.sub) + scale_fill_manual(values = cell_colours)
dev.off()

sc.data.sub = subset(sc.data, orig.ident == "SN_2739")
png(paste0(poly_folder, "/SN_2739_segmentation.png"), width = 4000, height = 4000)
ImageDimPlot(sc.data.sub) + scale_fill_manual(values = cell_colours)
dev.off()

sc.data.sub = subset(sc.data, orig.ident == "D_7678")
png(paste0(poly_folder, "/D_7678_segmentation.png"), width = 4000, height = 4000)
ImageDimPlot(sc.data.sub) + scale_fill_manual(values = cell_colours)
dev.off()

sc.data.sub = subset(sc.data, orig.ident == "D_6414")
png(paste0(poly_folder, "/D_6414_segmentation.png"), width = 4000, height = 4000)
ImageDimPlot(sc.data.sub) + scale_fill_manual(values = cell_colours)
dev.off()

sc.data.sub = subset(sc.data, orig.ident == "D_6446")
png(paste0(poly_folder, "/D_6446_segmentation.png"), width = 4000, height = 4000)
ImageDimPlot(sc.data.sub) + scale_fill_manual(values = cell_colours)
dev.off()



