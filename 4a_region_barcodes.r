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
library(MAST)
library(ggrepel)

#### Guides and Manuals ####

## https://bioconductor.org/packages/release/bioc/vignettes/granulator/inst/doc/granulator.html#reference-profiles





#### Helper Functions ####

theme_SL2 <- function () { 
  theme_bw() %+replace% 
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
      plot.background = element_blank(), 
      legend.background = element_rect(fill="transparent", colour=NA),
      legend.key = element_rect(fill="transparent", colour=NA),
      plot.title = element_text(size=12, margin = margin(b = 5),hjust=0,vjust=0.5, face="bold"),
      title = element_text(size = 12, margin = margin(b = 5),hjust=0,vjust=0.5, face="bold"),
      axis.text.y = element_text(size = 10, margin = margin(r = 5),hjust=1,vjust=0.5, face="bold",colour="black"),
      axis.text.x = element_text(size = 10, margin = margin(t = 5),hjust=0.5,vjust=1, face="bold",colour="black"),
      axis.title.y = element_text(size = 11, margin = margin(r = 10),angle = 90,hjust=0.5,vjust=0.5, face="bold"),
      axis.title.x = element_text(size = 11, margin = margin(t = 10),hjust=0.5,vjust=1, face="bold"),
      legend.text=element_text(size=11, face="bold"),
      legend.title=element_blank(), 
      legend.key.size=unit(2.5,"line"),
      plot.margin=unit(c(0.4,0.4,0.4,0.4), "cm")
    )
}


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


make_volcano_plot <- function(de_annotated, sig_up_colour,sig_down_colour, non_sig_colour, dot_size, dot_transparency, sig_up_name, sig_down_name, non_sig_name, x_axis_label, y_axis_label,legend_position, topn, label_size) 
{
  
  # create a direction column
  de_volcano = de_annotated
  de_volcano$direction = "c"
  
  de_sig_up = subset(de_volcano, avg_log2FC > 0 & p_val_adj < 0.05)
  de_sig_down = subset(de_volcano, avg_log2FC < 0 & p_val_adj < 0.05)
  de_non_sig = subset(de_volcano, p_val_adj >= 0.05)
  
  if (nrow(de_sig_up) > 0)
  {
    de_sig_up$direction = "a"
  }
  if (nrow(de_sig_down) > 0)
  {
    de_sig_down$direction = "b"
  }
  
  de_volcano = rbind(de_non_sig, de_sig_down, de_sig_up)
  
  
  # get genes to label
  de_topn = de_volcano[order(de_volcano$p_val),]
  de_topn = de_topn[1:topn,]

  
  # get the limits
  limits = c(-max(abs(de_volcano$avg_log2FC)), max(abs(de_volcano$avg_log2FC)))
  
  # make the plot
  ggp = ggplot(data=de_volcano, aes(x=avg_log2FC, y=-log10(p_val), colour=direction)) +
    geom_point(size=dot_size,alpha=dot_transparency) +
    geom_label_repel(data=de_topn, aes(label=rownames(de_topn)), size=label_size, force = 1,box.padding = 1, show.legend = FALSE, colour = "black") +
    scale_color_manual(breaks=c("a","b","c"), values=c(sig_up_colour ,sig_down_colour, non_sig_colour), labels=c(sig_up_name, sig_down_name, non_sig_name)) +
    labs(x=x_axis_label, y=y_axis_label) +
    xlim(limits) +
    theme_bw() + 
    theme(legend.position=legend_position, legend.title = element_blank())
  
  return(ggp)
}





#### Load data ####

## wd
setwd("C:/Users/JohnJ/Desktop/analysis/Amber_Bozward_Visium_2023/analysis")

## load
liver_st = readRDS(file = "liver_integrated_light_psudo_decon_linneages.rds")
barcodes_deep = read.table("4a_region_barcodes/Deep_Parenchyma_vs_Deep_NonParenchyma.csv", header = TRUE, sep=",", row.names=1)
barcodes_all = read.table("4a_region_barcodes/Parenchyma_vs_NonParenchyma.csv", header = TRUE, sep=",", row.names=1)



#### Convert ABs Slice IDs to JCs ####


#AB1 = JC6
#AB2 = JC7
#AB3 = JC8
#AB4 = JC1
#AB5 = JC2
#AB6 = JC9
#AB7 = JC10
#AB8 = JC11
#AB9 = JC4
#AB10 = JC5
#AB11 = JC3


# barcodes_deep
row.names(barcodes_deep) = gsub("_10","_JC5",row.names(barcodes_deep))
row.names(barcodes_deep) = gsub("_11","_JC3",row.names(barcodes_deep))
row.names(barcodes_deep) = gsub("_1","_JC6",row.names(barcodes_deep))
row.names(barcodes_deep) = gsub("_2","_JC7",row.names(barcodes_deep))
row.names(barcodes_deep) = gsub("_3","_JC8",row.names(barcodes_deep))
row.names(barcodes_deep) = gsub("_4","_JC1",row.names(barcodes_deep))
row.names(barcodes_deep) = gsub("_5","_JC2",row.names(barcodes_deep))
row.names(barcodes_deep) = gsub("_6","_JC9",row.names(barcodes_deep))
row.names(barcodes_deep) = gsub("_7","_JC10",row.names(barcodes_deep))
row.names(barcodes_deep) = gsub("_8","_JC11",row.names(barcodes_deep))
row.names(barcodes_deep) = gsub("_9","_JC4",row.names(barcodes_deep))

row.names(barcodes_deep) = gsub("_JC6","_6",row.names(barcodes_deep))
row.names(barcodes_deep) = gsub("_JC7","_7",row.names(barcodes_deep))
row.names(barcodes_deep) = gsub("_JC8","_8",row.names(barcodes_deep))
row.names(barcodes_deep) = gsub("_JC1","_1",row.names(barcodes_deep))
row.names(barcodes_deep) = gsub("_JC2","_2",row.names(barcodes_deep))
row.names(barcodes_deep) = gsub("_JC9","_9",row.names(barcodes_deep))
row.names(barcodes_deep) = gsub("_JC10","_10",row.names(barcodes_deep))
row.names(barcodes_deep) = gsub("_JC11","_11",row.names(barcodes_deep))
row.names(barcodes_deep) = gsub("_JC4","_4",row.names(barcodes_deep))
row.names(barcodes_deep) = gsub("_JC5","_5",row.names(barcodes_deep))
row.names(barcodes_deep) = gsub("_JC3","_3",row.names(barcodes_deep))


# barcodes_all
row.names(barcodes_all) = gsub("_10","_JC5",row.names(barcodes_all))
row.names(barcodes_all) = gsub("_11","_JC3",row.names(barcodes_all))
row.names(barcodes_all) = gsub("_1","_JC6",row.names(barcodes_all))
row.names(barcodes_all) = gsub("_2","_JC7",row.names(barcodes_all))
row.names(barcodes_all) = gsub("_3","_JC8",row.names(barcodes_all))
row.names(barcodes_all) = gsub("_4","_JC1",row.names(barcodes_all))
row.names(barcodes_all) = gsub("_5","_JC2",row.names(barcodes_all))
row.names(barcodes_all) = gsub("_6","_JC9",row.names(barcodes_all))
row.names(barcodes_all) = gsub("_7","_JC10",row.names(barcodes_all))
row.names(barcodes_all) = gsub("_8","_JC11",row.names(barcodes_all))
row.names(barcodes_all) = gsub("_9","_JC4",row.names(barcodes_all))


row.names(barcodes_all) = gsub("_JC6","_6",row.names(barcodes_all))
row.names(barcodes_all) = gsub("_JC7","_7",row.names(barcodes_all))
row.names(barcodes_all) = gsub("_JC8","_8",row.names(barcodes_all))
row.names(barcodes_all) = gsub("_JC1","_1",row.names(barcodes_all))
row.names(barcodes_all) = gsub("_JC2","_2",row.names(barcodes_all))
row.names(barcodes_all) = gsub("_JC9","_9",row.names(barcodes_all))
row.names(barcodes_all) = gsub("_JC10","_10",row.names(barcodes_all))
row.names(barcodes_all) = gsub("_JC11","_11",row.names(barcodes_all))
row.names(barcodes_all) = gsub("_JC4","_4",row.names(barcodes_all))
row.names(barcodes_all) = gsub("_JC5","_5",row.names(barcodes_all))
row.names(barcodes_all) = gsub("_JC3","_3",row.names(barcodes_all))



#### Add to meta data ####

# tidy up the groups - deep
for (row_index in 1:nrow(barcodes_deep))
{
  group = barcodes_deep[row_index,1]
  if (group == "Deep parenchyma"){barcodes_deep[row_index,1] = "Deep parenchyma"}
  else if (group == "Deep non-parenchyma"){barcodes_deep[row_index,1] = "Deep non-parenchyma"}
  else {barcodes_deep[row_index,1] = "Other"}
  
}
names(barcodes_deep) = "barcodes_deep"


# tidy up the groups - all
for (row_index in 1:nrow(barcodes_all))
{
  group = barcodes_all[row_index,1]
  if (group == "Parenchyma"){barcodes_all[row_index,1] = "Parenchyma"}
  else if (group == "Non-parenchyma"){barcodes_all[row_index,1] = "Non-parenchyma"}
  else {barcodes_all[row_index,1] = "Other"}
}
names(barcodes_all) = "barcodes_all"

# add to the seurat object
liver_st = AddMetaData(liver_st, barcodes_deep)
liver_st = AddMetaData(liver_st, barcodes_all)


# cells in each group
cells_deep_parenchyma = row.names(subset(barcodes_deep,barcodes_deep=="deep parenchyma"))
cells_deep_non_parenchyma = row.names(subset(barcodes_deep,barcodes_deep=="deep non-parenchyma"))
cells_parenchyma = row.names(subset(barcodes_all,barcodes_all=="Parenchyma"))
cells_non_parenchyma = row.names(subset(barcodes_all,barcodes_all=="Non-parenchyma"))



#### Check regions ####

# Plot UMAP
ggp = DimPlot(liver_st, reduction = "umap", group.by = "barcodes_deep", cells.highlight = list(cells_deep_parenchyma, cells_deep_non_parenchyma), sizes.highlight = 2, pt.size = 1) + scale_colour_manual(labels=c("other", "deep non-parenchyma","deep parenchyma"), values=c("grey", "blue", "red"))
save_plot(ggp, 1000, 1000, "4a_region_barcodes/barcodes_UMAP_deep.png")

ggp = DimPlot(liver_st, reduction = "umap", group.by = "barcodes_all", cells.highlight = list(cells_parenchyma, cells_non_parenchyma), sizes.highlight = 2, pt.size = 1) + scale_colour_manual(labels=c("other", "non-parenchyma","parenchyma"), values=c("grey", "blue", "red"))
save_plot(ggp, 1000, 1000, "4a_region_barcodes/barcodes_UMAP_all.png")


# Plots spatial
ggp = SpatialDimPlot(liver_st, group.by = "barcodes_deep", pt.size.factor = 6, stroke = NA, ncol=11)
save_plot(ggp, 1500, 6000, "4a_region_barcodes/barcodes_Spatial_deep.png")

ggp = SpatialDimPlot(liver_st, group.by = "barcodes_all", pt.size.factor = 6, stroke = NA, ncol=11)
save_plot(ggp, 1500, 6000, "4a_region_barcodes/barcodes_Spatial_all.png")



#### Differential Expression - MAST ####

## cells by disease
cells_by_sample = liver_st$orig.ident
cells_by_disease = cells_by_sample

# get disease status
for (index in 1:length(cells_by_disease))
{
  sample = as.character(cells_by_disease[index])
  sample = unlist(strsplit(sample, "_", fixed=TRUE))[1]
  cells_by_disease[index] = sample
}

cells_by_disease = data.frame(cells_by_disease)
names(cells_by_disease) = "disease"

# add to Seurat
liver_st = AddMetaData(liver_st, cells_by_disease)


## Parenchyma / Non-Parenchyma by disease

# deep
barcodes_deep_by_disease = barcodes_deep
for (index in 1:nrow(barcodes_deep))
{
  cell = row.names(barcodes_deep)[index]
  disease = cells_by_disease[cell,1]
  region = barcodes_deep_by_disease[cell,1]
  barcodes_deep_by_disease[cell,1] = paste(disease,region,sep=" ")
}
names(barcodes_deep_by_disease) = "barcodes_deep_by_disease"
liver_st = AddMetaData(liver_st, barcodes_deep_by_disease)

# all
barcodes_all_by_disease = barcodes_all
for (index in 1:nrow(barcodes_deep))
{
  cell = row.names(barcodes_all)[index]
  disease = cells_by_disease[cell,1]
  region = barcodes_all_by_disease[cell,1]
  barcodes_all_by_disease[cell,1] = paste(disease,region,sep=" ")
}
names(barcodes_all_by_disease) = "barcodes_all_by_disease"
liver_st = AddMetaData(liver_st, barcodes_all_by_disease)

## Do DE
liver_st = SetIdent(liver_st, value = "barcodes_deep")
de_deep = FindMarkers(liver_st, ident.1="Deep parenchyma", ident.2="Deep non-parenchyma", test.use = "MAST")

liver_st = SetIdent(liver_st, value = "barcodes_all")
de_all = FindMarkers(liver_st, ident.1="Parenchyma", ident.2="Non-parenchyma", test.use = "MAST")

liver_st = SetIdent(liver_st, value = "barcodes_deep_by_disease")
de_deep_parenchyma = FindMarkers(liver_st, ident.1="SN Deep parenchyma", ident.2="AH Deep parenchyma", test.use = "MAST")
de_deep_nonparenchyma = FindMarkers(liver_st, ident.1="SN Deep non-parenchyma", ident.2="AH Deep non-parenchyma", test.use = "MAST")

liver_st = SetIdent(liver_st, value = "barcodes_all_by_disease")
de_all_parenchyma = FindMarkers(liver_st, ident.1="SN Parenchyma", ident.2="AH Parenchyma", test.use = "MAST")
de_all_nonparenchyma = FindMarkers(liver_st, ident.1="SN Non-parenchyma", ident.2="AH Non-parenchyma", test.use = "MAST")

# save the DE results
write.table(de_deep, "4a_region_barcodes/DE_deep_parenchyma_VS_deep_non_parenchyma.csv", sep="\t", quote=FALSE)
write.table(de_all, "4a_region_barcodes/DE_parenchyma_VS_non_parenchyma.csv", sep="\t", quote=FALSE)
write.table(de_deep_parenchyma, "4a_region_barcodes/DE_deep_parenchyma_SN_VS_AL.csv", sep="\t", quote=FALSE)
write.table(de_deep_nonparenchyma, "4a_region_barcodes/DE_deep_non_parenchyma_SN_VS_AL.csv", sep="\t", quote=FALSE)
write.table(de_all_parenchyma, "4a_region_barcodes/DE_parenchyma_SN_VS_AL.csv", sep="\t", quote=FALSE)
write.table(de_all_nonparenchyma, "4a_region_barcodes/DE_non_parenchyma_SN_VS_AL.csv", sep="\t", quote=FALSE)



#### Volcano - MAST ####

height = 750
width = 750
sig_up_colour = "coral3"
sig_down_colour = "cornflowerblue"
non_sig_colour = "grey75"
dot_size = 1.5
dot_transparency = 1
sig_up_name = "upregulated"
sig_down_name = "downregulated"
non_sig_name = "non-significant"
topn = 100 # option is an integer or "ALL_SIG"
label_size = 3.5
x_axis_label = expression("log"[2]* " fold change")
y_axis_label = expression("-log"[10]* " p-value")
legend_position = "bottom"

ggp = make_volcano_plot(de_deep, sig_up_colour,sig_down_colour, non_sig_colour, dot_size, dot_transparency, sig_up_name, sig_down_name, non_sig_name, x_axis_label, y_axis_label,legend_position, topn, label_size) 
save_plot(ggp, height, width, "4a_region_barcodes/Volcano_deep_parenchyma_VS_deep_non_parenchyma.png")

ggp = make_volcano_plot(de_all, sig_up_colour,sig_down_colour, non_sig_colour, dot_size, dot_transparency, sig_up_name, sig_down_name, non_sig_name, x_axis_label, y_axis_label,legend_position, topn, label_size) 
save_plot(ggp, height, width, "4a_region_barcodes/Volcano_parenchyma_VS_non_parenchyma.png")

ggp = make_volcano_plot(de_deep_parenchyma, sig_up_colour,sig_down_colour, non_sig_colour, dot_size, dot_transparency, sig_up_name, sig_down_name, non_sig_name, x_axis_label, y_axis_label,legend_position, topn, label_size) 
save_plot(ggp, height, width, "4a_region_barcodes/Volcano_deep_parenchyma_SN_VS_AL.png")

ggp = make_volcano_plot(de_deep_nonparenchyma, sig_up_colour,sig_down_colour, non_sig_colour, dot_size, dot_transparency, sig_up_name, sig_down_name, non_sig_name, x_axis_label, y_axis_label,legend_position, topn, label_size) 
save_plot(ggp, height, width, "4a_region_barcodes/Volcano_deep_non_parenchyma_SN_VS_AL.png")

ggp = make_volcano_plot(de_all_parenchyma, sig_up_colour,sig_down_colour, non_sig_colour, dot_size, dot_transparency, sig_up_name, sig_down_name, non_sig_name, x_axis_label, y_axis_label,legend_position, topn, label_size) 
save_plot(ggp, height, width, "4a_region_barcodes/Volcano_parenchyma_SN_VS_AL.png")

ggp = make_volcano_plot(de_all_nonparenchyma, sig_up_colour,sig_down_colour, non_sig_colour, dot_size, dot_transparency, sig_up_name, sig_down_name, non_sig_name, x_axis_label, y_axis_label,legend_position, topn, label_size) 
save_plot(ggp, height, width, "4a_region_barcodes/Volcano_non_parenchyma_SN_VS_AL.png")



#### Differential Expression - PSEUDO BULK ####

saveRDS(liver_st, file = "liver_integrated_light_psudo_decon_linneages_regions.rds")
