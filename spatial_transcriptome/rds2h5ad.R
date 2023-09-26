library(Seurat)
library(SeuratData)
library(SeuratDisk)

args <- commandArgs(T)
SeuObj <- readRDS(args[1])
SaveH5Seurat(SeuObj, filename = paste0(args[2],".h5Seurat"))
Convert(paste0(args[2],".h5Seurat"), dest = "h5ad")
unlink(paste0(args[2],".h5Seurat"))
