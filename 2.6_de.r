library(ggplot2)
library(Seurat)
library(SeuratData)
library(dplyr)
library(patchwork)
library(org.Hs.eg.db)
library(clusterProfiler)
library(msigdbr)
library(Azimuth)
library(dittoSeq)
library(harmony)
library(scDblFinder)
library(EnhancedVolcano)
library(GOSemSim)
library(enrichplot)
library(ggupset)
library(DESeq2)
library(CellChat)


#### Parameters ####


# set memory available
options(future.globals.maxSize = 100000 * 1024^2)


#### Load data ####


# folders
setwd("C:/Users/JohnJ/Desktop/analysis/Ye_Oo/Amber_B_CosMx_2024/analysis/combined_all/")


# load sc object
sc.data = readRDS("final.RDS")

# Remove underscore
sc.data$P.VS.NP = as.factor(sc.data$P_VS_NP)
sc.data$final.cell.types = as.factor(sc.data$final_cell_types)

# set idents
Idents(sc.data) = sc.data$sample_group


#### DE function ####

do_de = function(out.path, seurat.dat, g1, g2, group)
{
  
  # make out folder
  if (!file.exists(out.path))
  {
    dir.create(out.path)
  }
  

  # get de
  seurat.dat = NormalizeData(seurat.dat)
  seurat.dat = ScaleData(seurat.dat)
  de.sub = FindMarkers(seurat.dat, ident.1 = g1, ident.2 = g2)
  
  # save
  write.table(de.sub, paste0(out.path, "/DE_",g1,"_vs_",g2,".csv"), quote=FALSE, sep = "\t")
  
  # volcano
  ggp = EnhancedVolcano(de.sub , rownames(de.sub),x ="avg_log2FC",y ="p_val_adj", FCcutoff = 1, pCutoff = 0.05)
  png(paste0(out.path, "/Volcano_",g1,"_vs_",g2,".png"), height = 1000, width = 1000)
  print(ggp)
  dev.off()
  
  
  ## heatmap
  de.sub.sig = subset(de.sub, p_val_adj < 0.05 & abs(avg_log2FC) > 1)
  ggp = DoHeatmap(seurat.dat, features = row.names(de.sub.sig), group.by = group)
  png(paste0(out.path, "/Heatmap_",g1,"_vs_",g2,".png"), height = 1000, width = 1000)
  print(ggp)
  dev.off()
  
  
  ## violins
  de.sub.sig.up = subset(de.sub, p_val_adj < 0.05 & avg_log2FC > 1)
  de.sub.sig.dn = subset(de.sub, p_val_adj < 0.05 & avg_log2FC < 1)
  
  ggp = VlnPlot(seurat.dat, features = row.names(de.sub.sig.up)[1:10], split.by = group, pt.size = 0, ncol = 5)
  png(paste0(out.path, "/Violin_up_",g1,"_vs_",g2,".png"), height = 1000, width = 1000)
  print(ggp)
  dev.off()
  
  ggp = VlnPlot(seurat.dat, features = row.names(de.sub.sig.dn)[1:10], split.by = group, pt.size = 0, ncol = 5)
  png(paste0(out.path, "/Violin_dn_",g1,"_vs_",g2,".png"), height = 1000, width = 1000)
  print(ggp)
  dev.off()

}

#### DE handler function ####

de_handler_sample_group = function(cell.type, region)
{
  sc.data.sub = subset(sc.data, P.VS.NP == region & final.cell.types == cell.type & (sample_group == "AIH" | sample_group == "D"))
  do_de(paste0("DE/",region,"/",cell.type,"/"), sc.data.sub, "AIH", "D", "sample_group")
  sc.data.sub = subset(sc.data, P.VS.NP == region & final.cell.types == cell.type & (sample_group == "SN" | sample_group == "D"))
  do_de(paste0("DE/",region,"/",cell.type,"/"), sc.data.sub, "SN", "D", "sample_group")
  sc.data.sub = subset(sc.data, P.VS.NP == region & final.cell.types == cell.type & (sample_group == "SN" | sample_group == "AIH"))
  do_de(paste0("DE/",region,"/",cell.type,"/"), sc.data.sub, "SN", "AIH", "sample_group")
}

#### NP - sample group ####

levels(sc.data$final.cell.types)

# run de
de_handler_sample_group("Hepatocyte", "NP")
de_handler_sample_group("Hepatic stellate", "NP")
de_handler_sample_group("IS Hepatocyte", "NP")
de_handler_sample_group("Type 2 LSEC", "NP")
de_handler_sample_group("Epithelial", "NP")
de_handler_sample_group("Type 1 LSEC", "NP")
de_handler_sample_group("Macrophage", "NP")
de_handler_sample_group("T cell", "NP")
de_handler_sample_group("Monocyte", "NP")
de_handler_sample_group("NK cell", "NP")
de_handler_sample_group("Plasma cell", "NP")
de_handler_sample_group("B cell", "NP")
de_handler_sample_group("Neutrophil", "NP")


#### P - sample group ####

levels(sc.data$final.cell.types)

# run de
de_handler_sample_group("Hepatocyte", "P")
de_handler_sample_group("Hepatic stellate", "P")
de_handler_sample_group("IS Hepatocyte", "P")
de_handler_sample_group("Type 2 LSEC", "P")
de_handler_sample_group("Epithelial", "P")
de_handler_sample_group("Type 1 LSEC", "P")
de_handler_sample_group("Macrophage", "P")
de_handler_sample_group("T cell", "P")
de_handler_sample_group("Monocyte", "P")
de_handler_sample_group("NK cell", "P")
de_handler_sample_group("Plasma cell", "P")
de_handler_sample_group("B cell", "P")
de_handler_sample_group("Neutrophil", "P")




#### Cell frequencies by region ####


## Cell frequencies

# frequencies 
freqs = data.frame(table(sc.data@meta.data$orig.ident, sc.data$final.cell.types, sc.data$P.VS.NP))
names(freqs) = c("sample","predominant.cell","tissue","count")
write.table(freqs, "DE/Frequences.csv", quote=FALSE, sep = "\t", row.names = FALSE)





#### Pseudobulk ####


# Aggregate for all cell types at once
pseudo = AggregateExpression(sc.data, group.by = c("orig.ident", "final.cell.types", "P.VS.NP"), return.seurat = TRUE)
pseudo = data.frame(pseudo@assays$RNA$counts)

write.table(pseudo, "pseudo_bulk/counts.csv", sep = "\t", quote = FALSE, row.names = TRUE)






#### CellChat function ####


do_cellchat = function(seurat.dat, cell_types, assay)
{
  future::plan("multisession", workers = 4)
  cellchat_data = createCellChat(object = seurat.dat, group.by = cell_types, assay = assay)
  cellchat_DB = CellChatDB.human
  cellchat_DB = subsetDB(cellchat_DB)
  cellchat_data@DB = cellchat_DB
  cellchat_data = subsetData(cellchat_data)
  cellchat_data = identifyOverExpressedGenes(cellchat_data)
  cellchat_data = identifyOverExpressedInteractions(cellchat_data)
  cellchat_data = computeCommunProb(cellchat_data , type = "triMean")
  cellchat_data = filterCommunication(cellchat_data, min.cells = 10)
  cellchat_data = computeCommunProbPathway(cellchat_data)
  cellchat_data = aggregateNet(cellchat_data)
  cellchat_data = netAnalysis_computeCentrality(cellchat_data , slot.name = "netP")
  return(cellchat_data)
}

do_cellchat_plots = function(cellchat_data, out.path, principal.cell, total.types)
{
  # make out folder
  if (!file.exists(out.path))
  {
    dir.create(out.path)
  }
  
  
  # circle plots
  groupSize = as.numeric(table(cellchat_data@idents)) 
  
  ggp = netVisual_circle(cellchat_data@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions", vertex.label.cex = 2)
  png(paste0(out.path, "/global_count_circos.png"), height = 1000, width = 1000)
  print(ggp)
  dev.off()
  
  ggp = netVisual_circle(cellchat_data@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", vertex.label.cex = 2)
  png(paste0(out.path, "/global_strength_circos.png"), height = 1000, width = 1000)
  print(ggp)
  dev.off()
  
  # heatmaps
  ggp = netVisual_heatmap(cellchat_data, font.size = 15)
  png(paste0(out.path, "/global_scount_heatmap.png"), height = 1000, width = 1000)
  print(ggp)
  dev.off()
  
  ggp = netVisual_heatmap(cellchat_data, measure = "weight", font.size = 15)
  png(paste0(out.path, "/global_strength_heatmap.png"), height = 1000, width = 1000)
  print(ggp)
  dev.off()
  
  
  # principal cell
  ggp = netVisual_bubble(cellchat_data, sources.use = principal.cell, targets.use = c(1:total.types), font.size = 20)
  png(paste0(out.path, "/",principal.cell,"_outgoing_bubble.png"), height = 750, width = 750)
  print(ggp)
  dev.off()
  
  ggp = netVisual_bubble(cellchat_data, sources.use = c(1:total.types), targets.use = principal.cell, font.size = 20)
  png(paste0(out.path, "/",principal.cell,"_incoming_bubble.png"), height = 750, width = 750)
  print(ggp)
  dev.off()
  
  # Signal scatter
  ggp = netAnalysis_signalingRole_scatter(cellchat_data, label.size = 5)
  png(paste0(out.path, "/global_scatter.png"), height = 750, width = 750)
  print(ggp)
  dev.off()
  
}




compare_cellchat = function(g1, g2, cc1, cc2)
{
  cellchat.list = list("g1" = cc1, "g2" = cc2) 
  names(cellchat.list)[1] = g1
  names(cellchat.list)[2] = g2
  
  merged_cellchat = mergeCellChat(cellchat.list, add.names = names(cellchat.list))
  return(merged_cellchat)
}


do_compare_cellchat_plots = function(merged_cellchat, out.path, principal.cell, total.types)
{
  # make out folder
  if (!file.exists(out.path))
  {
    dir.create(out.path)
  }
  
  # number of differential interactions
  gg1 = compareInteractions(merged_cellchat, show.legend = F, group = c(1,2))
  gg2 = compareInteractions(merged_cellchat, show.legend = F, group = c(1,2), measure = "weight")
  png(paste0(out.path, "/barcharts.png"), height = 500, width = 500)
  print(gg1 + gg2)
  dev.off()
  
  
  # differential interaction circos
  ggp = netVisual_diffInteraction(merged_cellchat , weight.scale = T) 
  png(paste0(out.path, "/diff_count_circos.png"), height = 500, width = 500)
  print(ggp)
  dev.off()
  
  ggp = netVisual_diffInteraction(merged_cellchat , weight.scale = T, measure = "weight")
  png(paste0(out.path, "/diff_strength_circos.png"), height = 500, width = 500)
  print(ggp)
  dev.off()
  
  
  # differential interaction heatmap
  g1 = netVisual_heatmap(merged_cellchat, font.size = 15)
  g2 = netVisual_heatmap(merged_cellchat , measure = "weight", font.size = 15) 
  png(paste0(out.path, "/diff_heatmap.png"), height = 750, width = 1000)
  print(g1 + g2)
  dev.off()
  
  
  # principal cell
  ggp = netVisual_bubble(merged_cellchat, sources.use = principal.cell, targets.use = c(1:total.types),  comparison = c(1, 2), font.size = 20)
  png(paste0(out.path, "/",principal.cell,"_outgoing_bubble.png"), height = 750, width = 750)
  print(ggp)
  dev.off()
  
  ggp = netVisual_bubble(merged_cellchat, sources.use = c(1:total.types), targets.use = principal.cell,  comparison = c(1, 2), font.size = 20)
  png(paste0(out.path, "/",principal.cell,"_incoming_bubble.png"), height = 750, width = 750)
  print(ggp)
  dev.off()
  
}


#### CellChat  ####


## NP

# get CCI
sc.data.sub = subset(sc.data, P.VS.NP == "NP" & sample_group == "AIH")
cell_chat_1 = do_cellchat(sc.data.sub, "final.cell.types", "RNA")
gc()
sc.data.sub = subset(sc.data, P.VS.NP == "NP" & sample_group == "SN")
cell_chat_2 = do_cellchat(sc.data.sub, "final.cell.types", "RNA")
gc()
sc.data.sub = subset(sc.data, P.VS.NP == "NP" & sample_group == "D")
cell_chat_3 = do_cellchat(sc.data.sub, "final.cell.types", "RNA")
gc()

# plot
do_cellchat_plots(cell_chat_1, "cellchat/NP_AIH", "B cell", 13)
do_cellchat_plots(cell_chat_1, "cellchat/NP_AIH", "Macrophage", 13)
do_cellchat_plots(cell_chat_2, "cellchat/NP_SN", "B cell", 13)
do_cellchat_plots(cell_chat_2, "cellchat/NP_SN", "Macrophage", 13)
do_cellchat_plots(cell_chat_3, "cellchat/NP_D", "B cell", 13)
do_cellchat_plots(cell_chat_3, "cellchat/NP_D", "Macrophage", 13)

# save objects
saveRDS(cell_chat_1, file="cellchat/CellChat_NP_AIH.RDS")
saveRDS(cell_chat_2, file="cellchat/CellChat_NP_SN.RDS")
saveRDS(cell_chat_3, file="cellchat/CellChat_NP_D.RDS")


# differential
cellchat_c = compare_cellchat("AIH", "D", cell_chat_1, cell_chat_3)
gc()
do_compare_cellchat_plots(cellchat_c, "cellchat/NP_AIH_vs_D", "B cell", 13)
do_compare_cellchat_plots(cellchat_c, "cellchat/NP_AIH_vs_D", "Macrophage", 13)

cellchat_c = compare_cellchat("SN", "D", cell_chat_2, cell_chat_3)
gc()
do_compare_cellchat_plots(cellchat_c, "cellchat/NP_SN_vs_D", "B cell", 13)
do_compare_cellchat_plots(cellchat_c, "cellchat/NP_SN_vs_D", "Macrophage", 13)

cellchat_c = compare_cellchat("SN", "AIH", cell_chat_2, cell_chat_1)
gc()
do_compare_cellchat_plots(cellchat_c, "cellchat/NP_SN_vs_AIH", "B cell", 13)
do_compare_cellchat_plots(cellchat_c, "cellchat/NP_SN_vs_AIH", "Macrophage", 13)





## P

# get CCI
sc.data.sub = subset(sc.data, P.VS.NP == "P" & sample_group == "AIH")
cell_chat_1 = do_cellchat(sc.data.sub, "final.cell.types", "RNA")
gc()
sc.data.sub = subset(sc.data, P.VS.NP == "P" & sample_group == "SN")
cell_chat_2 = do_cellchat(sc.data.sub, "final.cell.types", "RNA")
gc()
sc.data.sub = subset(sc.data, P.VS.NP == "P" & sample_group == "D")
cell_chat_3 = do_cellchat(sc.data.sub, "final.cell.types", "RNA")
gc()

# plot
do_cellchat_plots(cell_chat_1, "cellchat/P_AIH", "B cell", 13)
do_cellchat_plots(cell_chat_1, "cellchat/P_AIH", "Macrophage", 13)
do_cellchat_plots(cell_chat_2, "cellchat/P_SN", "B cell", 13)
do_cellchat_plots(cell_chat_2, "cellchat/P_SN", "Macrophage", 13)
do_cellchat_plots(cell_chat_3, "cellchat/P_D", "B cell", 13)
do_cellchat_plots(cell_chat_3, "cellchat/P_D", "Macrophage", 13)

# save objects
saveRDS(cell_chat_1, file="cellchat/CellChat_NP_AIH.RDS")
saveRDS(cell_chat_2, file="cellchat/CellChat_NP_SN.RDS")
saveRDS(cell_chat_3, file="cellchat/CellChat_NP_D.RDS")


# differential
cellchat_c = compare_cellchat("AIH", "D", cell_chat_1, cell_chat_3)
gc()
do_compare_cellchat_plots(cellchat_c, "cellchat/P_AIH_vs_D", "B cell", 13)
do_compare_cellchat_plots(cellchat_c, "cellchat/P_AIH_vs_D", "Macrophage", 13)

cellchat_c = compare_cellchat("SN", "D", cell_chat_2, cell_chat_3)
gc()
do_compare_cellchat_plots(cellchat_c, "cellchat/P_SN_vs_D", "B cell", 13)
do_compare_cellchat_plots(cellchat_c, "cellchat/P_SN_vs_D", "Macrophage", 13)

cellchat_c = compare_cellchat("SN", "AIH", cell_chat_2, cell_chat_1)
gc()
do_compare_cellchat_plots(cellchat_c, "cellchat/P_SN_vs_AIH", "B cell", 13)
do_compare_cellchat_plots(cellchat_c, "cellchat/P_SN_vs_AIH", "Macrophage", 13)




#### Niches ####

poly_folder = "polygons/"
  
# plots
sc.data.sub = subset(sc.data, orig.ident == "AIH_1607")
sc.data.sub = BuildNicheAssay(sc.data.sub,"Nano","final.cell.types",assay = "niche",cluster.name = "niches",neighbors.k = 20,niches.k = 10)
png(paste0(poly_folder, "/AIH_1607_Niches.png"), width = 4000, height = 4000)
ImageDimPlot(sc.data.sub, group.by = "niches")
dev.off()

sc.data.sub = subset(sc.data, orig.ident == "AIH_6079")
sc.data.sub = BuildNicheAssay(sc.data.sub,"Nano.2","final.cell.types",assay = "niche",cluster.name = "niches",neighbors.k = 20,niches.k = 10)
png(paste0(poly_folder, "/AIH_6079_Niches.png"), width = 4000, height = 4000)
ImageDimPlot(sc.data.sub, group.by = "niches")
dev.off()

sc.data.sub = subset(sc.data, orig.ident == "SN_606")
sc.data.sub = BuildNicheAssay(sc.data.sub,"Nano.3","final.cell.types",assay = "niche",cluster.name = "niches",neighbors.k = 20,niches.k = 10)
png(paste0(poly_folder, "/SN_606_Niches.png"), width = 4000, height = 4000)
ImageDimPlot(sc.data.sub, group.by = "niches")
dev.off()

sc.data.sub = subset(sc.data, orig.ident == "SN_1042")
sc.data.sub = BuildNicheAssay(sc.data.sub,"Nano.4","final.cell.types",assay = "niche",cluster.name = "niches",neighbors.k = 20,niches.k = 10)
png(paste0(poly_folder, "/SN_1042_Niches.png"), width = 4000, height = 4000)
ImageDimPlot(sc.data.sub, group.by = "niches")
dev.off()

sc.data.sub = subset(sc.data, orig.ident == "SN_2739")
sc.data.sub = BuildNicheAssay(sc.data.sub,"Nano.5","final.cell.types",assay = "niche",cluster.name = "niches",neighbors.k = 20,niches.k = 10)
png(paste0(poly_folder, "/SN_2739_Niches.png"), width = 4000, height = 4000)
ImageDimPlot(sc.data.sub, group.by = "niches")
dev.off()

sc.data.sub = subset(sc.data, orig.ident == "D_7678")
sc.data.sub = BuildNicheAssay(sc.data.sub,"Nano.6","final.cell.types",assay = "niche",cluster.name = "niches",neighbors.k = 20,niches.k = 10)
png(paste0(poly_folder, "/D_7678_Niches.png"), width = 4000, height = 4000)
ImageDimPlot(sc.data.sub, group.by = "niches")
dev.off()

sc.data.sub = subset(sc.data, orig.ident == "D_6414")
sc.data.sub = BuildNicheAssay(sc.data.sub,"Nano.7","final.cell.types",assay = "niche",cluster.name = "niches",neighbors.k = 20,niches.k = 10)
png(paste0(poly_folder, "/D_6414_Niches.png"), width = 4000, height = 4000)
ImageDimPlot(sc.data.sub, group.by = "niches")
dev.off()

sc.data.sub = subset(sc.data, orig.ident == "D_6446")
sc.data.sub = BuildNicheAssay(sc.data.sub,"Nano.8","final.cell.types",assay = "niche",cluster.name = "niches",neighbors.k = 20,niches.k = 10)
png(paste0(poly_folder, "/D_6446_Niches.png"), width = 4000, height = 4000)
ImageDimPlot(sc.data.sub, group.by = "niches")
dev.off()







#### CDKN1A Spatial feature plot ####



poly_folder = "de/spatial_plots/"

# plots
sc.data.sub = subset(sc.data, orig.ident == "AIH_1607")
png(paste0(poly_folder, "/AIH_1607_CDKN1A.png"), width = 4000, height = 4000)
ImageFeaturePlot(sc.data.sub, features = "CDKN1A", axes = TRUE, cols = c("blue","yellow","red"))
dev.off()

sc.data.sub = subset(sc.data, orig.ident == "AIH_6079")
png(paste0(poly_folder, "/AIH_6079_CDKN1A.png"), width = 4000, height = 4000)
ImageFeaturePlot(sc.data.sub, features = "CDKN1A", axes = TRUE, cols = c("blue","yellow","red"))
dev.off()

sc.data.sub = subset(sc.data, orig.ident == "SN_606")
png(paste0(poly_folder, "/SN_606_CDKN1A.png"), width = 4000, height = 4000)
ImageFeaturePlot(sc.data.sub, features = "CDKN1A", axes = TRUE, cols = c("blue","yellow","red"))
dev.off()

sc.data.sub = subset(sc.data, orig.ident == "SN_1042")
png(paste0(poly_folder, "/SN_1042_CDKN1A.png"), width = 4000, height = 4000)
ImageFeaturePlot(sc.data.sub, features = "CDKN1A", axes = TRUE, cols = c("blue","yellow","red"))
dev.off()

sc.data.sub = subset(sc.data, orig.ident == "SN_2739")
png(paste0(poly_folder, "/SN_2739_CDKN1A.png"), width = 4000, height = 4000)
ImageFeaturePlot(sc.data.sub, features = "CDKN1A", axes = TRUE, cols = c("blue","yellow","red"))
dev.off()

sc.data.sub = subset(sc.data, orig.ident == "D_7678")
png(paste0(poly_folder, "/D_7678_CDKN1A.png"), width = 4000, height = 4000)
ImageFeaturePlot(sc.data.sub, features = "CDKN1A", axes = TRUE, cols = c("blue","yellow","red"))
dev.off()

sc.data.sub = subset(sc.data, orig.ident == "D_6414")
png(paste0(poly_folder, "/D_6414_CDKN1A.png"), width = 4000, height = 4000)
ImageFeaturePlot(sc.data.sub, features = "CDKN1A", axes = TRUE, cols = c("blue","yellow","red"))
dev.off()

sc.data.sub = subset(sc.data, orig.ident == "D_6446")
png(paste0(poly_folder, "/D_6446_CDKN1A.png"), width = 4000, height = 4000)
ImageFeaturePlot(sc.data.sub, features = "CDKN1A", axes = TRUE, cols = c("blue","yellow","red"))
dev.off()


