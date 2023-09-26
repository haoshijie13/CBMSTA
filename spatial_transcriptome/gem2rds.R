args<-commandArgs(T)
print("script start")
mat2=data.table::fread(args[1])
prefix=args[2]
bin=args[3]
bin=as.numeric(bin)
mat2$x <- round(mat2$x/bin,0)
mat2$y <- round(mat2$y/bin,0)
mat2$cluster <- paste0(mat2$x,mat2$y)

cell=1:length(unique(mat2$cluster))
names(cell)=paste0("cell_",unique(mat2$cluster))
gene=1:length(unique(mat2$geneID))
names(gene)=unique(mat2$geneID)
library(Matrix)
mat=sparseMatrix(j=cell[paste0("cell_",mat2$cluster)],i=gene[mat2$geneID],x=mat2$MIDCount)
rownames(mat)=names(gene)
colnames(mat)=names(cell)
library(Seurat)
obj=CreateSeuratObject(counts = mat,min.features = 0,min.cells = 0)
coord=aggregate(mat2[,c("x","y")],by=list("cluster"=mat2$cluster),mean)
rownames(coord)=paste0("cell_",coord$cluster)
obj@meta.data[,c("coor_x","coor_y")]=coord[rownames(obj@meta.data),c("x","y")]
if (!dir.exists("STrds")){
    print("create STrds")
    dir.create("STrds")
}
saveRDS(obj,file=paste0("STrds/",prefix,".rds"))
print(paste0("file save in STrds/",prefix,".rds")
print("mission completed")
