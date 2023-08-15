library(Seurat)  
library(ggsci)
library(ggplot2)
library(patchwork)
library(cowplot)
library(pheatmap)
library(psych)

args<-commandArgs(T)

macaqueobj <- readRDS(args[1])
marmosetobj <- readRDS(args[2])
miceobj <- readRDS(args[3])

result <- read.csv(args[4],header = TRUE,row.names = 1)
homo <- read.table(args[5],header = TRUE)
markerdf <- merge(homo,result,by.x = "mouseGene",by.y = "gene")

macaqueobj_sub <- subset(macaqueobj,cluster%in%c("Astrocyte","Bergmann","Choroid","Endothelial_stalk","Fibroblast","Golgi","Granule","Microglia","MLI1","MLI2","ODC","OPC","PLI","Purkinje","UBC"))
marmosetobj_sub <- subset(marmosetobj,cluster%in%c("Astrocyte","Bergmann","Choroid","Endo_stalk","Fibroblast","Golgi","Granule","Mic","MLI1","MLI2","OLG","OPC","PLI","Purkinje","UBC"))
miceobj_sub <- subset(miceobj,celltype%in%c("Astrocyte","Bergmann","Choroid","Endothelial","Fibroblast","Golgi","Granule","Microglia","MLI1","MLI2","ODC","OPC","PLI","Purkinje","UBC"))

ct1 <- c("Astrocyte","Bergmann","Choroid","Endo","Fibroblast","Golgi","Granule","Microglia","MLI1","MLI2","ODC","OPC","PLI","Purkinje","UBC")

genelist1 <- c()
genelist2 <- c()
genelist3 <- c()

for (cell in ct1) {
  gene11 <- markerdf$mouseGene[which(markerdf$type=="common"&markerdf$celltype==cell)]
  gene12 <- markerdf$macaqueGene[which(markerdf$type=="common"&markerdf$celltype==cell)]
  gene13 <- markerdf$marmosetGene[which(markerdf$type=="common"&markerdf$celltype==cell)]
  genelist1 <- c(genelist1,gene11)
  genelist2 <- c(genelist2,gene12)
  genelist3 <- c(genelist3,gene13)
  
}

for (cell in ct1) {
  gene21 <- markerdf$mouseGene[which(markerdf$type=="macaque_uniq"&markerdf$celltype==cell)]
  gene22 <- markerdf$macaqueGene[which(markerdf$type=="macaque_uniq"&markerdf$celltype==cell)]
  gene23 <- markerdf$marmosetGene[which(markerdf$type=="macaque_uniq"&markerdf$celltype==cell)]
  genelist1 <- c(genelist1,gene21)
  genelist2 <- c(genelist2,gene22)
  genelist3 <- c(genelist3,gene23)
}

for (cell in ct1) {
  gene31 <- markerdf$mouseGene[which(markerdf$type=="marmoset_uniq"&markerdf$celltype==cell)]
  gene32 <- markerdf$macaqueGene[which(markerdf$type=="marmoset_uniq"&markerdf$celltype==cell)]
  gene33 <- markerdf$marmosetGene[which(markerdf$type=="marmoset_uniq"&markerdf$celltype==cell)]
  genelist1 <- c(genelist1,gene31)
  genelist2 <- c(genelist2,gene32)
  genelist3 <- c(genelist3,gene33)
}

for (cell in ct1) {
  gene41 <- markerdf$mouseGene[which(markerdf$type=="mice_uniq"&markerdf$celltype==cell)]
  gene42 <- markerdf$macaqueGene[which(markerdf$type=="mice_uniq"&markerdf$celltype==cell)]
  gene43 <- markerdf$marmosetGene[which(markerdf$type=="mice_uniq"&markerdf$celltype==cell)]
  genelist1 <- c(genelist1,gene41)
  genelist2 <- c(genelist2,gene42)
  genelist3 <- c(genelist3,gene43)
}


macaqueobj_sub <- ScaleData(macaqueobj_sub, features = genelist2)
macaqueobj_sub$cluster <- factor(macaqueobj_sub$cluster,levels = c("Astrocyte","Bergmann","Choroid","Endothelial_stalk","Fibroblast","Golgi","Granule","Microglia","MLI1","MLI2","ODC","OPC","PLI","Purkinje","UBC"))
Idents(macaqueobj_sub) <- "cluster"
p1 <- DoHeatmap(macaqueobj_sub, features = genelist2,angle = 0,hjust = 0.5,label = FALSE,
                group.colors = c("#DC050C","#FB8072","#1965B0","#7BAFDE","#FF7F00","#FDB462","#E7298A","#33A02C","#B2DF8A","#55A1B1","#8DD3C7","#A6761D","#7570B3","#BEAED4","#6a3d9a"))&
  scale_fill_gradientn(colors = c("white", "#fde0dd", "#c51b8a"))&theme(axis.text.y = element_blank())

marmosetobj_sub <- ScaleData(marmosetobj_sub, features = genelist3)
marmosetobj_sub$cluster <- factor(marmosetobj_sub$cluster,levels = c("Astrocyte","Bergmann","Choroid","Endo_stalk","Fibroblast","Golgi","Granule","Mic","MLI1","MLI2","OLG","OPC","PLI","Purkinje","UBC"))
Idents(marmosetobj_sub) <- "cluster"
p2 <- DoHeatmap(marmosetobj_sub, features = genelist3,angle = 0,hjust = 0.5,label = FALSE,
                group.colors = c("#DC050C","#FB8072","#1965B0","#7BAFDE","#FF7F00","#FDB462","#E7298A","#33A02C","#B2DF8A","#55A1B1","#8DD3C7","#A6761D","#7570B3","#BEAED4","#6a3d9a"))&
  scale_fill_gradientn(colors = c("white", "#fee0d2", "#de2d26"))&theme(axis.text.y = element_blank())

miceobj_sub <- ScaleData(miceobj_sub, features = genelist1)
miceobj_sub$celltype <- factor(miceobj_sub$celltype,levels = c("Astrocyte","Bergmann","Choroid","Endothelial","Fibroblast","Golgi","Granule","Microglia","MLI1","MLI2","ODC","OPC","PLI","Purkinje","UBC"))
Idents(miceobj_sub) <- "celltype"

p3 <- DoHeatmap(miceobj_sub, features = genelist1,angle = 0,hjust = 0.5,label = FALSE,
                group.colors = c("#DC050C","#FB8072","#1965B0","#7BAFDE","#FF7F00","#FDB462","#E7298A","#33A02C","#B2DF8A","#55A1B1","#8DD3C7","#A6761D","#7570B3","#BEAED4","#6a3d9a"))&
  scale_fill_gradientn(colors = c("white", "#fff7bc", "#d95f0e"))&theme(axis.text.y = element_blank())
