library(Seurat)
library(parallel)
args<-commandArgs(T)

obj2 <- readRDS(args[1])


A <- c('LobuleI','LobuleII','LobuleIII','LobuleIV_V','LobuleIV','LobuleV')
P <- c('LobuleVI','Simplex lobule','CrusI','LobuleVII','LobuleVIII','CrusII','LobuleIX','Copula','Parmedian lobule')

AH <- c('LobuleI','LobuleII','LobuleIII','LobuleIV_V','LobuleVI','Simplex lobule','CrusI','LobuleIV','LobuleV')
PH <- c('LobuleVII','LobuleVIII','CrusII','LobuleIX','Copula','Parmedian lobule')

V <- c('LobuleI','LobuleII','LobuleIII','LobuleIV_V','LobuleVI','LobuleVII','LobuleVIII','LobuleIX','LobuleIV','LobuleV')
H <- c('CrusII','Copula','Parmedian lobule','Simplex lobule','CrusI')

X <- c('LobuleX')
PFF <- c("Paraflocculus","flocculus")

obj2@meta.data$ident1 <- obj2@meta.data$lobulename2
obj2@meta.data[obj2@meta.data$lobulename2%in%X,'ident1'] <- 'X'
obj2@meta.data[obj2@meta.data$lobulename2%in%PFF,'ident1'] <- 'PFF'
Idents(obj2) <- 'ident1'


obj2@active.assay <- "RNA"

#mouse2 specific
if(args[2]=='mouse2'){
obj2$test <- '1'
obj2@meta.data[obj2@meta.data$nCount_RNA<50,'test'] <- '2'
obj2 <- subset(obj2,subset=test=='1')
}

obj2 <- NormalizeData(obj2)
obj2 <- ScaleData(obj2,features=rownames(obj2))
FindMarker <- function(ident1){

marker <- FindMarkers(obj2,logfc.threshold=0.05,only.pos = TRUE,ident.1=ident1,slot="scale.data")
#marker <- FindMarkers(obj2,logfc.threshold=0.05,only.pos = TRUE,ident.1=strsplit(ident1[1],"/")[[1]],ident.2=strsplit(ident1[2],"/")[[1]],slot="data")
marker$gene <- rownames(marker)
marker$ident1 <- ident1
return(marker)
}

#ident1s <- unique(obj2$ident1)
ident1s <- c('X','PFF')
#ident1s <- list(c('A','P','ident1'),c('P','A','ident1'),
#                c('AH','PH','ident2'),c('PH','AH','ident2'),
#                c('V','H','ident3'),c('H','V','ident3'),
#                c('PFF','rest/X','ident4'),c('X','rest/PFF','ident4'),
#                c('rest/X','PFF','ident4'),c('rest/PFF','X','ident4'))

cl<-makeCluster(2, type="FORK")
marker.list <- parLapply(cl, ident1s,FindMarker)
stopCluster(cl)

marker_raw <- marker.list[[1]]
for (i in 2:length(marker.list)){
    marker_raw <- rbind(marker_raw,marker.list[[i]])
    }

write.csv(marker_raw,paste0(args[2],".scale.lobule.marker.csv"),quote = FALSE)
