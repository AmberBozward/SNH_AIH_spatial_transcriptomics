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



#### Load data ####

## wd
setwd("C:/Users/JohnJ/Desktop/analysis/Amber_Bozward_Visium_2023/analysis")

## load
liver_st = readRDS(file = "liver_integrated_light_psudo_decon_linneages_regions.rds")


#### Pseudo ####

# gets the cell types
cell_groups_deep = data.frame(liver_st$barcodes_deep)
cell_groups_all = data.frame(liver_st$barcodes_all)
cell_samples = data.frame(liver_st$orig.ident)

names(cell_groups_deep) = "cell_type"
names(cell_groups_all) = "cell_type"
names(cell_samples) = "sample"

# get the groups
groups = c("Deep parenchyma","Deep non-parenchyma","Parenchyma","Non-parenchyma")
samples = c("AH_1607", "AH_2709", "AH_7538", "AH_6079", "AH_7014", "SN_606", "SN_827", "SN_1042", "SN_2739", "SN_3303", "SN_4393")
samples_decon = paste(expand.grid(a=samples, b=groups)$a, expand.grid(a=samples, b=groups)$b,sep="_")

# gets the sc EM
st_counts = liver_st@assays$Spatial@counts

# create a table to store the pseudo counts
st_pseudo_counts = data.frame(matrix(0,ncol=length(samples_decon), nrow=nrow(st_counts)))
names(st_pseudo_counts) = samples_decon
row.names(st_pseudo_counts) = row.names(st_counts)

## gets pseudobulk
counter = 0
for (col_index in 1:ncol(st_counts))
{
  # progress update
  if (counter %% 100 == 0) {print(counter)}
  counter = counter + 1
  
  # get all the right data
  cell_data = st_counts[,col_index]
  cell_id = gsub("\\.","-",dimnames(st_counts)[[2]][col_index])
  cell_group_deep = as.character(cell_groups_deep[cell_id,])
  cell_group_all = as.character(cell_groups_all[cell_id,])
  cell_sample = as.character(cell_samples[cell_id,])
  sample_decon_deep = paste(cell_sample, cell_group_deep, sep="_")
  sample_decon_all = paste(cell_sample, cell_group_all, sep="_")
  
  # update the counts
  if (cell_group_deep %in% c("Deep parenchyma","Deep non-parenchyma","Parenchyma", "Non-parenchyma"))
  {
    st_pseudo_counts[,sample_decon_deep] = st_pseudo_counts[,sample_decon_deep] + cell_data
  }
  
  if (cell_group_all %in% c("Deep parenchyma","Deep non-parenchyma","Parenchyma", "Non-parenchyma"))
  {
    st_pseudo_counts[,sample_decon_all] = st_pseudo_counts[,sample_decon_all] + cell_data
  }
}

write.table(st_pseudo_counts, "4b_region_barcodes_pseudo/counts_pseudo.csv", sep="\t", quote=FALSE)
