#gemc output cell_colname region_colname
args<-commandArgs(T)
print("script start")
mat2 <- data.table::fread(args[1])
prefix <- args[2]
#bin <- args[3]
#bin <- as.numeric(bin)
#mat2$x <- round(mat2$x/bin,0)
#mat2$y <- round(mat2$y/bin,0)
mat2 <- as.data.frame(mat2)

#qianqian
mat2$geneID <- mat2$gene
mat2$MIDCount <- mat2$umi_count

mat2$cluster <- mat2[,args[3]]
mat2 <- mat2[mat2$cluster!=0,]
region <- aggregate(mat2[,args[4]],by=list("cellID"=paste0("cell_", mat2$cluster)),FUN=function(x){names(sort(table(x),decreasing = T))[1]})
rownames(region) <- region$cellID


cell <- 1:length(unique(mat2$cluster))
names(cell) <- paste0("cell_", unique(mat2$cluster))
gene <- 1:length(unique(mat2$geneID))
names(gene) <- unique(mat2$geneID)
library(Matrix)
mat <- sparseMatrix(j=cell[paste0("cell_", mat2$cluster)], i=gene[mat2$geneID], x=mat2$MIDCount)
rownames(mat) <- names(gene)
colnames(mat) <- names(cell)
library(Seurat)
obj <- CreateSeuratObject(counts = mat, min.features = 0, min.cells = 0)
coord <- aggregate(mat2[, c("x", "y")], by=list("cluster"=mat2$cluster), mean)
rownames(coord) <- paste0("cell_", coord$cluster)
obj@meta.data[,c("coor_x","coor_y")] <- coord[rownames(obj@meta.data), c("x", "y")]
obj$region <- region[rownames(obj@meta.data),2]
if (!dir.exists("STrds")){
    print("create STrds")
    dir.create("STrds")
}
saveRDS(obj,file=paste0("STrds/",prefix,".rds"))
print(paste0("file save in STrds/",prefix,".rds"))
if (!dir.exists("STh5ad")){
    print("create STh5ad")
    dir.create("STh5ad")
}

SaveH5Seurat(obj, filename = paste0("STh5ad/",prefix,".h5Seurat"))
Convert(paste0("STh5ad/",prefix,".h5Seurat"), dest = "h5ad")
unlink(paste0("STh5ad/",prefix,".h5Seurat"))

print("mission completed")

