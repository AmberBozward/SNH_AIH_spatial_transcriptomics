#### Load Libraries ####

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(devtools)
library(hdf5r)
library(SCDC)
library(Biobase)

#### Guides and Manuals ####

# https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_07_spatial.html
# https://satijalab.org/seurat/articles/spatial_vignette.html
# https://www.livercellatlas.org/download.php
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04574-5



#### Load integrated tables for further analysis ####

## wd
setwd("C:/Users/JohnJ/Desktop/analysis/Amber_Bozward_Visium_2023/analysis")

## load
liver_st = readRDS(file = "liver_integrated_light.rds")
liver_sc = readRDS(file = "public/GSE158723/single_cell_integrated.rds")

# check
DimPlot(liver_sc, reduction = "umap", label = TRUE, pt.size = 0.5)



#### Convert Sub-Cell Types to Linneages ####

## manually match
new.cluster.ids = c("Hep", "Hep", "Hep", "Hep", "Hep", "Hep", "Hep", "MP", "MP", "Hep", "Hep", "CC", "EC", "Hep", "LP", "LP", "EC", "MP", "MP", "HSC", "HSC", "LP", "LP", "LP", "MP")

# rename
names(new.cluster.ids) = levels(liver_sc)
liver_sc = RenameIdents(liver_sc, new.cluster.ids)

# check
DimPlot(liver_sc, reduction = "umap", label = TRUE, pt.size = 0.5)






#### Spatial Deconvolution ####


## Get the markers
markers_sc = FindAllMarkers(liver_sc, only.pos = TRUE, logfc.threshold = 0.1,
                            test.use = "wilcox", min.pct = 0.05, min.diff.pct = 0.1, max.cells.per.ident = 200,
                            return.thresh = 0.05, assay = "RNA")

# Filter for genes that are also present in the ST data
markers_sc = markers_sc[markers_sc$gene %in% rownames(liver_st), ]


# Select top 20 genes per cluster, select top by first p-value, then absolute
# diff in pct, then quota of pct.
markers_sc$pct.diff = markers_sc$pct.1 - markers_sc$pct.2
markers_sc$log.pct.diff = log2((markers_sc$pct.1 * 99 + 1)/(markers_sc$pct.2 * 99 +
                                                               1))
markers_sc %>%
  group_by(cluster) %>%
  top_n(-2500, p_val) %>%
  top_n(1250, pct.diff) %>%
  top_n(500, log.pct.diff) -> top20
m_feats = unique(as.character(top20$gene))

# add cell types to SC metadata
liver_sc = AddMetaData(liver_sc, liver_sc@active.ident, col.name = "cell_types")


## Deconvolution

## Create expression sets
eset_SC = ExpressionSet(assayData = as.matrix(liver_sc@assays$RNA@counts[m_feats,]), phenoData = AnnotatedDataFrame(liver_sc@meta.data))
eset_ST = ExpressionSet(assayData = as.matrix(liver_st@assays$Spatial@counts[m_feats,]), phenoData = AnnotatedDataFrame(liver_st@meta.data))

## deconvolve
deconvolution_crc = SCDC::SCDC_prop(bulk.eset = eset_ST, sc.eset = eset_SC, ct.varname = "cell_types",ct.sub = as.character(unique(eset_SC$cell_types)))

# check
head(deconvolution_crc$prop.est.mvw)

# add deconvolution results to spatial data
liver_st@assays[["SCDC"]] = CreateAssayObject(data = t(deconvolution_crc$prop.est.mvw))




# Seems to be a bug in SeuratData package that the key is not set and any
# plotting function etc. will throw an error.
if (length(liver_st@assays$SCDC@key) == 0) {
  liver_st@assays$SCDC@key = "scdc_"
}
DefaultAssay(liver_st) <- "SCDC"



# plot
plot = SpatialFeaturePlot(liver_st, features = c("Hep.1", "Hep.2", "Hep.3", "Hep.4", "Hep.5", "Hep.6", "Hep.7", "Hep.8", "Hep.9","Hep.10"), pt.size.factor = 6, ncol = 9,crop = TRUE, image.alpha = 0, stroke = NA)
png("3f_spatial_deconvolution_linneages_light/light_integrated_deonvolution_HEP_spatial.png", height = 6000, width = 4000)
print(plot)
dev.off()

plot = SpatialFeaturePlot(liver_st, features = c("MP.1", "MP.2", "MP.3", "MP.4", "CC.1", "EC.1", "EC.2", "EC.3", "LP.1", "LP.2", "LP.3", "LP.4", "LP.5", "HSC.1", "HSC.2"), pt.size.factor = 6, ncol = 9,crop = TRUE, image.alpha = 0, stroke = NA)
png("3f_spatial_deconvolution_linneages_light/light_integrated_deonvolution_Kuppher_spatial.png", height = 6000, width = 4000)
print(plot)
dev.off()




#### Save table for analysis ####

saveRDS(liver_st, file = "liver_integrated_light_decon.rds")
