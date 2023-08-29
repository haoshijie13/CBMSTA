library(Seurat)
library(harmony)
library(patchwork)
library(data.table)
library(dplyr)

options(future.globals.maxSize = 1000 * 1024^2)

#region <- read.csv("macaque_ROI_label.csv",header = TRUE)
#scRNA_harmony$region <- as.character(region[1,])

scRNA_harmony <- readRDS("macaque_all_gene_eachlobule500_3D_noharmony_noL0.rds")

scRNA_harmony <- NormalizeData(scRNA_harmony, normalization.method = "LogNormalize", scale.factor = 10000)

scRNA_harmony <- FindVariableFeatures(scRNA_harmony, selection.method = "vst", nfeatures = 2000)
scRNA_harmony <- ScaleData(scRNA_harmony)
scRNA_harmony <- RunPCA(scRNA_harmony)

scRNA_harmony <- FindNeighbors(scRNA_harmony, dims = 1:20)
scRNA_harmony <- FindClusters(scRNA_harmony, resolution = 0.5)
scRNA_harmony <- RunUMAP(scRNA_harmony, dims = 1:20)
scRNA_harmony <- RunTSNE(scRNA_harmony, dims = 1:20)

saveRDS(scRNA_harmony,"macaque_all_gene_eachlobule500_3D_noharmony_noL0_nor.rds")

##plot
#group_by_cluster
#plot1 = DimPlot(scRNA_harmony, reduction = "umap", label=T)
#group_by_sample
#plot2 = DimPlot(scRNA_harmony, reduction = "umap", group.by='region_2')

#plot3 = DimPlot(scRNA_harmony, reduction = "tsne", label=T)
#group_by_sample
#plot4 = DimPlot(scRNA_harmony, reduction = "tsne", group.by='region_2')

#combinate
#pdf("macaque_all_gene_3D_harmony_eachlobule500.cluster.pdf",width=8,height=6)
#plot1
#plot2
#plot3
#plot4
#dev.off()

