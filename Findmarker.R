library(Seurat)
library(ggplot2)
args<-commandArgs(T)

obj <- readRDS(args[1])

gene_overlap <- read.table(args[2],header=TRUE)

gene <- gene_overlap[,1]
obj2 <- obj[gene,]


Idents(obj2) <- args[3]
marker <- FindMarkers(obj2, ident.1 = args[4], only.pos = TRUE, min.pct = 0, logfc.threshold = 0)

write.csv(marker,args[5],quote = FALSE)
