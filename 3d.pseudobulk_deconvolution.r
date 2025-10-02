#### Load Libraries ####

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(devtools)
library(hdf5r)
library(Biobase)
library(DESeq2)
library(granulator)
library(RColorBrewer)
library(stringr)

#### Guides and Manuals ####

## https://bioconductor.org/packages/release/bioc/vignettes/granulator/inst/doc/granulator.html#reference-profiles



#### Function to run DESEq2 ####

run_deseq2 = function(countdata, sample_sheet)
{
  

  # Convert to matrix
  countdata = as.matrix(countdata)
  
  # Filter for genes with a median >= 0 read
  countdata = subset(countdata,apply(countdata, 1, mean) >= 1)
  
  # Assign condition
  condition = factor(sample_sheet$GROUP)
  
  ## Analysis with DESeq2 ----------------------------------------------------
  
  # Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
  coldata = data.frame(row.names=colnames(countdata), condition)
  dds = DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
  dds
  
  # Run the DESeq pipeline
  dds = DESeq(dds)
  
  #Get the expression matrix:
  
  # Get the EM
  resdata = round(as.data.frame(counts(dds, normalized=TRUE)),2)

  return(resdata)
  
}






#### Load integrated tables for further analysis ####

## wd
setwd("C:/Users/JohnJ/Desktop/analysis/Amber_Bozward_Visium_2023/analysis")

## load
liver_st = readRDS(file = "liver_integrated_light.rds")
liver_sc = readRDS(file = "public/GSE158723/single_cell_integrated.rds")

# check
DimPlot(liver_sc, reduction = "umap", label = TRUE, pt.size = 0.5)



#### Pseudobulk of public single cell data ####

# gets the cell types
cell_types = data.frame(liver_sc@active.ident)
names(cell_types) = "cell_type"

# gets the sc EM
sc_counts = liver_sc@assays$RNA@counts

# makes the pseudobulk output
sc_pseudo_counts = data.frame(matrix(0,ncol=length(unique(cell_types$cell_type)), nrow=nrow(sc_counts)))
names(sc_pseudo_counts) = unique(cell_types$cell_type)
row.names(sc_pseudo_counts) = row.names(sc_counts)
  
  
# gets pseudobulk
for (col_index in 1:ncol(sc_counts))
{
  cell_data = sc_counts[,col_index]
  cell_id = gsub("\\.","-",dimnames(sc_counts)[[2]][col_index])
  cell_type = as.character(cell_types[cell_id,"cell_type"])
  sc_pseudo_counts[,cell_type] = sc_pseudo_counts[,cell_type] + cell_data
}

## normalises pseudobulk using DESEq2

# make the sample sheet
pseudo_ss = data.frame(names(sc_pseudo_counts))
names(pseudo_ss) = "SAMPLE"
row.names(pseudo_ss) = pseudo_ss$names.sc_pseudo_counts.
pseudo_ss$GROUP = sample.int(2,nrow(pseudo_ss), replace = TRUE)

# deseq
sc_pseudo_em = run_deseq2(sc_pseudo_counts, pseudo_ss)

# tidy up
rm(sc_pseudo_counts, liver_sc, sc_counts)
gc()




#### Pseudobulk of spatial data ####

# gets the cell types
cell_types = data.frame(liver_st@active.ident)
names(cell_types) = "cell_type"

# gets the sc EM
st_counts = liver_st@assays$Spatial@counts

# makes the pseudobulk output
st_pseudo_counts = data.frame(matrix(0,ncol=length(unique(cell_types$cell_type)), nrow=nrow(st_counts)))
names(st_pseudo_counts) = unique(cell_types$cell_type)
row.names(st_pseudo_counts) = row.names(st_counts)


# gets pseudobulk
for (col_index in 1:ncol(st_counts))
{
  cell_data = st_counts[,col_index]
  cell_id = gsub("\\.","-",dimnames(st_counts)[[2]][col_index])
  cell_type = as.character(cell_types[cell_id,"cell_type"])
  st_pseudo_counts[,cell_type] = st_pseudo_counts[,cell_type] + cell_data
}


## normalises pseudobulk using DESEq2

# make the sample sheet
pseudo_ss = data.frame(names(st_pseudo_counts))
names(pseudo_ss) = "SAMPLE"
row.names(pseudo_ss) = pseudo_ss$names.st_pseudo_counts.
pseudo_ss$GROUP = sample.int(2,nrow(pseudo_ss), replace = TRUE)

# deseq
st_pseudo_em = run_deseq2(st_pseudo_counts, pseudo_ss)

# tidy up
rm(st_pseudo_counts, liver_st, st_counts)
gc()



#### Deconvolute the Spatial Data ####

# plot the similarity matrix
plot = plot_similarity(sigMatrix=as.matrix(sc_pseudo_em))
png("3d.pseudobulk_deconvolution/similarity_matrix.png", height = 750, width = 750)
print(plot)
dev.off()

# do deconvolution
decon = deconvolute(m = as.matrix(st_pseudo_em), sigMatrix = as.matrix(sc_pseudo_em))

plot = plot_proportions(deconvoluted = decon, method = 'svr')
plot = plot_deconvolute(deconvoluted = decon, scale = TRUE, labels = TRUE)


# save RDS
#saveRDS(st_pseudo_em, file = "3d.pseudobulk_deconvolution/pseudo_st.rds")
#saveRDS(sc_pseudo_em, file = "3d.pseudobulk_deconvolution/pseudo_sc.rds")
#saveRDS(decon, file = "3d.pseudobulk_deconvolution/pseudo_decon.rds")


## add cell fractions to spatial data

# get the fractions - we use SVR.
fractions = decon$proportions$svr_sig1

# get the cell types
cell_types = data.frame(liver_st@active.ident)
names(cell_types) = "cell_type"

# get the spot table
spot_cell_fractions = data.frame(matrix(0,nrow=nrow(cell_types),ncol=ncol(fractions)))
names(spot_cell_fractions) = names(fractions)
row.names(spot_cell_fractions) = row.names(cell_types)

# fill the table
for (row_index in 1:nrow(spot_cell_fractions))
{
  spot_id = row.names(spot_cell_fractions)[row_index]
  spot_cluster = as.character(cell_types[row_index,"cell_type"])
  row_fraction = fractions[spot_cluster,]
  spot_cell_fractions[spot_id,] = row_fraction
}

# Add the fractions to the seurat metadata
liver_st = AddMetaData(liver_st, spot_cell_fractions)

# save
#saveRDS(liver_st,"liver_integrated_light_psudo_decon.rds")




#### Spatial by Cell Type ####


plot = SpatialFeaturePlot(liver_st, features = names(spot_cell_fractions), stroke = NA, pt.size.factor = 6)
png("3d.pseudobulk_deconvolution/Light_Pseudo_Decon.png", height = 12000, width = 5000)
print(plot)
dev.off()


#### UMAP by Cell Type ####

# Plot UMAP
plot = FeaturePlot(liver_st, reduction = "umap", features = names(fractions), ncol=5)
png("3d.pseudobulk_deconvolution/Light_Pseudo_Decon_UMAP.png", height = 1500, width = 1500)
print(plot)
dev.off()



