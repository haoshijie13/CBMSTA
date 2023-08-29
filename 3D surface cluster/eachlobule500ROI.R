library(Seurat)
library(harmony)
library(patchwork)
library(data.table)
library(dplyr)

options(future.globals.maxSize = 1000 * 1024^2)

df <- fread("macaque_all_gene.csv",sep = ",",header = TRUE)
gene <- read.csv("Macaque_allslice_common_genes.csv",header = TRUE,row.names = 1)

dfmiss <- fread("macaque_miss_gene.csv",sep = ",",header = FALSE)
genemiss <- read.csv("Macaque_allslice_miss_genes.csv",header = TRUE,row.names = 1)

dim(df)
length(gene$gene)


df <- as.data.frame(df)
dfmiss <- as.data.frame(dfmiss)
dim(dfmiss)

rownames(df) <- gene$gene
colnames(dfmiss) <- genemiss$gene

label <- read.csv("macaque_gene_eachlobe_500label.csv",header=FALSE)
names(label) <- "label"
dim(label)

dfmisst <- as.data.frame(t(dfmiss))
colnames(dfmisst) <- colnames(df)

dim(dfmisst)
dfa <- rbind(df,dfmisst)

dft <- as.data.frame(t(dfa))
dft <- cbind(label,dft)


result <- aggregate(.~label,data=dft,mean)
rownames(result) <- result$label
resultt <- as.data.frame(t(result))
resultt <- resultt[-1,]

dim(resultt)

scRNA_harmony <- CreateSeuratObject(resultt,project ="macaque_5H")
saveRDS(scRNA_harmony,"macaque_all_gene_3D_eachlobule500_intail.rds")

#rownames(region) <- region$label
#scRNA_harmony@meta.data$region <- region[rownames(scRNA_harmony@meta.data),]$region
#scRNA_harmony$region_2 <- paste0("region_",scRNA_harmony$region)
#region <- read.csv("macaque_ROI_label.csv",header = TRUE)
#scRNA_harmony$region <- as.character(region[1,])

scRNA_harmony <- FindVariableFeatures(scRNA_harmony, selection.method = "vst", nfeatures = 2000)
scRNA_harmony <- ScaleData(scRNA_harmony)
scRNA_harmony <- RunPCA(scRNA_harmony)

scRNA_harmony <- FindNeighbors(scRNA_harmony, dims = 1:20)
scRNA_harmony <- FindClusters(scRNA_harmony, resolution = 0.5)
scRNA_harmony <- RunUMAP(scRNA_harmony, dims = 1:20)
scRNA_harmony <- RunTSNE(scRNA_harmony, dims = 1:20)

saveRDS(scRNA_harmony,"macaque_all_gene_eachlobule500_3D_noharmony.rds")

##作图
#group_by_cluster
plot1 = DimPlot(scRNA_harmony, reduction = "umap", label=T)
#group_by_sample
#plot2 = DimPlot(scRNA_harmony, reduction = "umap", group.by='region_2')

plot3 = DimPlot(scRNA_harmony, reduction = "tsne", label=T)
#group_by_sample
#plot4 = DimPlot(scRNA_harmony, reduction = "tsne", group.by='region_2')

#combinate
pdf("macaque_all_gene_3D_harmony_eachlobule500.cluster.pdf",width=8,height=6)
plot1
#plot2
plot3
#plot4
dev.off()

