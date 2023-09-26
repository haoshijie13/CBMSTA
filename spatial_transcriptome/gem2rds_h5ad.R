library(argparse)

parser <- ArgumentParser()
parser$add_argument("-g", "--gem", type="character", default='',
        help="input gem (gem/gem.gz)",
        metavar="xxx.gem")
parser$add_argument("-f", "--func", type="character", default='all',
        help="function to create h5ad, rds or all, default all",
        metavar="h5ad/rds/all")
parser$add_argument("-o", "--output", type="character", default='',
        help="output file name",
        metavar="output")
parser$add_argument("-C", "--columns", type="character", default=c('geneID', 'x', 'y', 'MIDCount', 'cell', 'region'), nargs='+',
        help="set the gem columns names, must include(geneID x y MIDCount) default: geneID x y MIDCount cell region",
        metavar="geneID x y MIDCount cell region")
parser$add_argument("-b", "--bin", type="character", default='50',
        help="input cell or bin(10/20/50/100/200/...), if you choose cellbin, your columns must include cell, default 50",
        metavar="(10/20/50/100/200/...)/cell")
parser$add_argument("-r", "--regions", type="character", default=c('region'), nargs='+',
        help="the columns name of region, default region",
        metavar="region1 region2 ...")
parser$add_argument("-c", "--coordinate", type="character", default=c('x', 'y'), nargs='+',
        help="the columns name of coordinate,must include(x y), default: x y",
        metavar="x y rx ry ...")
args <- parser$parse_args()

print("script start")
mat2 <- data.table::fread(args$gem)
prefix <- args$output
mat2 <- as.data.frame(mat2)

if (length(args$columns)==length(colnames(mat2))){
colnames(mat2) <- args$columns
}else{
    print("columns number must equal with -C")
    q()
}

if (args$bin!='cell'){
    bin <- as.numeric(args$bin)
    mat2$x <- round(mat2$x/bin,0)
    mat2$y <- round(mat2$y/bin,0)
    mat2$cell <- paste0(mat2$x,mat2$y)
    mat2$x <- mat2$x*bin
    mat2$y <- mat2$y*bin
} else if (args$bin=='cell'){
    if (!('cell'%in%args$columns)){
        print("Error, if you choose cellbin, your columns must include cell")
        q()
    }
    mat2 <- mat2[mat2$cell!=0,]
}

cell <- 1:length(unique(mat2$cell))
names(cell) <- paste0("cell_", unique(mat2$cell))
gene <- 1:length(unique(mat2$geneID))
names(gene) <- unique(mat2$geneID)
library(Matrix)
mat <- sparseMatrix(j=cell[paste0("cell_", mat2$cell)], i=gene[mat2$geneID], x=mat2$MIDCount)
rownames(mat) <- names(gene)
colnames(mat) <- names(cell)
library(Seurat)
obj <- CreateSeuratObject(counts = mat, min.features = 0, min.cells = 0)

#add coordinate to obj
if(length(args$coordinate)%%2!=0){
    print("Error, coordinate parse must include x y pairs")
    q()
}

if (length(args$coordinate[args$coordinate%in%args$columns])==length(args$coordinate)){
    for (i in 1:as.integer(length(args$coordinate)/2)){
        x_column <- args$coordinate[i*2-1]
        y_column <- args$coordinate[i*2]
        coord <- aggregate(mat2[, c(x_column, y_column)], by=list("cell"=mat2$cell), mean)
        coord[, c(x_column, y_column)] <- coord[, c(x_column, y_column)]/bin
        rownames(coord) <- paste0("cell_", coord$cell)
        obj@meta.data[,c(paste0("coor_",x_column),paste0("coor_",y_column))] <- coord[rownames(obj@meta.data), c(x_column, y_column)]
    }
}else{print("Error, coordinate parse must all in columns")
      q()
}

#add region to obj
if (args$regions%in%args$columns){
    if (length(args$regions[args$regions%in%args$columns])==length(args$regions)){
        for (i in args$regions){
            region <- aggregate(mat2[,i],by=list("cellID"=paste0("cell_", mat2$cell)),FUN=function(x){names(sort(table(x),decreasing = T))[1]})
            rownames(region) <- region$cellID
            obj@meta.data[,i] <- region[rownames(obj@meta.data),2]
        }
    }else{print("Error, region parse must all in columns")
          q()
    }
}
#save obj to rds
if (args$func=="rds"||args$func=="all"){
    if (!dir.exists("STrds")){
        print("create STrds")
        dir.create("STrds")
    }
    saveRDS(obj,file=paste0("STrds/",prefix,".",args$bin,".rds"))
    print(paste0("file save in STrds/",prefix,".",args$bin,".rds"))
}
#save obj to h5ad
if (args$func=="h5ad"||args$func=="all"){
    if (!dir.exists("STh5ad")){
        print("create STh5ad")
        dir.create("STh5ad")
    }
    library(SeuratDisk)
    SaveH5Seurat(obj, filename = paste0("STh5ad/",prefix,".",args$bin,".h5Seurat"))
    Convert(paste0("STh5ad/",prefix,".",args$bin,".h5Seurat"), dest = "h5ad")
    unlink(paste0("STh5ad/",prefix,".",args$bin,".h5Seurat"))
    print(paste0("file save in STh5ad/",prefix,".",args$bin,".h5ad"))
}

print("mission completed")

