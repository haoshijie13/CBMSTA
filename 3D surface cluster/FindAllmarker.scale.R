library(Seurat)
library(parallel)
args<-commandArgs(T)

obj2 <- readRDS(args[1])

#AP-Primary
A <- c('LobuleI','LobuleII','LobuleI_II','LobuleIII','LobuleIV_V','LobuleIV','LobuleV')
P <- c('LobuleVI','Simplex lobule','CrusI','LobuleVII','LobuleVIII','CrusII','LobuleIX','Copula','Parmedian lobule')

#AP-Horizontal
AH <- c('LobuleI','LobuleII','LobuleI_II','LobuleIII','LobuleIV_V','LobuleVI','Simplex lobule','CrusI','LobuleIV','LobuleV')
PH <- c('LobuleVII','LobuleVIII','CrusII','LobuleIX','Copula','Parmedian lobule')

#AP-Prepyramidal
AP <- c('LobuleI','LobuleII','LobuleI_II','LobuleIII','LobuleIV_V','LobuleVI','Simplex lobule','CrusI','LobuleIV','LobuleV','LobuleVII','CrusII')
PP <- c('LobuleVIII','LobuleIX','Copula','Parmedian lobule')

#AP-Secondary
AS <- c('LobuleI','LobuleII','LobuleI_II','LobuleIII','LobuleIV_V','LobuleVI','Simplex lobule','CrusI','LobuleIV','LobuleV','LobuleVII','CrusII','LobuleVIII','Parmedian lobule')
PS <- c('LobuleIX','Copula')

V <- c('LobuleI','LobuleII','LobuleI_II','LobuleIII','LobuleIV_V','LobuleVI','LobuleVII','LobuleVIII','LobuleIX','LobuleIV','LobuleV')
H <- c('CrusII','Copula','Parmedian lobule','Simplex lobule','CrusI')

X <- c('LobuleX')
PFF <- c("Paraflocculus","flocculus")


obj2@meta.data$ident1 <- 'initial'
obj2@meta.data[obj2@meta.data$lobulename%in%A,'ident1'] <- 'A'
obj2@meta.data[obj2@meta.data$lobulename%in%P,'ident1'] <- 'P'
obj2@meta.data[obj2@meta.data$lobulename%in%X,'ident1'] <- 'X'
obj2@meta.data[obj2@meta.data$lobulename%in%PFF,'ident1'] <- 'PFF'

obj2@meta.data$ident2 <- 'initial'
obj2@meta.data[obj2@meta.data$lobulename%in%AH,'ident2'] <- 'AH'
obj2@meta.data[obj2@meta.data$lobulename%in%PH,'ident2'] <- 'PH'
obj2@meta.data[obj2@meta.data$lobulename%in%X,'ident2'] <- 'X'
obj2@meta.data[obj2@meta.data$lobulename%in%PFF,'ident2'] <- 'PFF'

obj2@meta.data$ident3 <- 'initial'
obj2@meta.data[obj2@meta.data$lobulename%in%V,'ident3'] <- 'V'
obj2@meta.data[obj2@meta.data$lobulename%in%H,'ident3'] <- 'H'
obj2@meta.data[obj2@meta.data$lobulename%in%X,'ident3'] <- 'X'
obj2@meta.data[obj2@meta.data$lobulename%in%PFF,'ident3'] <- 'PFF'

obj2@meta.data$ident4 <- 'rest'
obj2@meta.data[obj2@meta.data$lobulename%in%X,'ident4'] <- 'X'
obj2@meta.data[obj2@meta.data$lobulename%in%PFF,'ident4'] <- 'PFF'

obj2@meta.data$ident5 <- 'initial'
obj2@meta.data[obj2@meta.data$lobulename%in%AP,'ident5'] <- 'AP'
obj2@meta.data[obj2@meta.data$lobulename%in%PP,'ident5'] <- 'PP'
obj2@meta.data[obj2@meta.data$lobulename%in%X,'ident5'] <- 'X'
obj2@meta.data[obj2@meta.data$lobulename%in%PFF,'ident5'] <- 'PFF'

obj2@meta.data$ident6 <- 'initial'
obj2@meta.data[obj2@meta.data$lobulename%in%AS,'ident6'] <- 'AS'
obj2@meta.data[obj2@meta.data$lobulename%in%PS,'ident6'] <- 'PS'
obj2@meta.data[obj2@meta.data$lobulename%in%X,'ident6'] <- 'X'
obj2@meta.data[obj2@meta.data$lobulename%in%PFF,'ident6'] <- 'PFF'


if(args[2]=='mouse2'){
obj2$test <- '1'
obj2@meta.data[obj2@meta.data$nCount_RNA<50,'test'] <- '2'
obj2 <- subset(obj2,subset=test=='1')
}

obj2@active.assay <- "RNA"
obj2 <- NormalizeData(obj2)
obj2 <- FindVariableFeatures(obj2,nfeatures = 5000)
obj2 <- ScaleData(obj2,features=rownames(obj2))
FindMarker <- function(ident1){

Idents(obj2) <- ident1[3]
marker <- FindMarkers(obj2,logfc.threshold=0.05,only.pos = TRUE,ident.1=strsplit(ident1[1],"/")[[1]],ident.2=strsplit(ident1[2],"/")[[1]],slot="scale.data")
#marker <- FindMarkers(obj2,logfc.threshold=0.05,only.pos = TRUE,ident.1=strsplit(ident1[1],"/")[[1]],ident.2=strsplit(ident1[2],"/")[[1]],slot="data")
marker$gene <- rownames(marker)
marker$ident1 <- ident1[1]
marker$ident2 <- ident1[2]
return(marker)
}

#ident1s <- unique(obj2$ident1)
ident1s <- list(c('A','P','ident1'),c('P','A','ident1'),
                c('AH','PH','ident2'),c('PH','AH','ident2'),
                c('V','H','ident3'),c('H','V','ident3'),
                c('AP','PP','ident5'),c('PP','AP','ident5'),
                c('AS','PS','ident6'),c('PS','AS','ident6'))

cl<-makeCluster(5, type="FORK")
marker.list <- parLapply(cl, ident1s,FindMarker)
stopCluster(cl)

marker_raw <- marker.list[[1]]
for (i in 2:length(marker.list)){
    marker_raw <- rbind(marker_raw,marker.list[[i]])
    }

write.csv(marker_raw,paste0(args[2],".scale.marker.csv"),quote = FALSE)
