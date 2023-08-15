args<-commandArgs(T)
homolist <- read.table(args[1],sep = "\t",header = TRUE,)
homolist <- homolist[,c("macaqueGene","marmosetGene","mouseGene")]

sp <- "macaque"
macaquefiles <- list.files(paste0(args[2],sp,"/"),pattern = ".csv")
df_a <- as.data.frame(matrix(nrow=0,ncol=1))
names(df_a) <- c("mouseGene")
for (mf in macaquefiles) {
  celltype <- sub("_positive.csv","",sub("macaque_cluster_marker_","",mf))
  print(celltype)
  df <- read.csv(paste0(args[2],mf),header = TRUE,row.names = 1)
  df$gene <- rownames(df)
  df <- df[df$p_val<0.01,c("gene","avg_log2FC","pct.1")]
  #设置在细胞类型的表达FC的阈值
  df <- df[df$avg_log2FC>0,]
  df <- df[df$pct.1>0.25,]
  
  df <- head(df[order(df$avg_log2FC,decreasing = TRUE),],200)
  df$avg_log2FC <- 2
  
  df <- df[,c("gene","avg_log2FC")]
  
  names(df) <- c("gene",paste0(celltype,"_avg_log2FC_",sp))
  df <- merge(homolist,df,by.x = "macaqueGene",by.y = "gene")
  df <- df[,c("mouseGene",paste0(celltype,"_avg_log2FC_",sp))]
  df_a <- merge(df_a,df,by = "mouseGene",all = TRUE)
}

sp <- "marmoset"
marmosetfiles <- list.files(paste0(args[3],sp,"_0802","/"),pattern = ".csv")
for (mf in marmosetfiles) {
  celltype <- sub("_positive.csv","",sub("marmoset_cluster_marker_","",mf))
  print(celltype)
  df <- read.csv(paste0(args[3],sp,"/",mf),header = TRUE,row.names = 1)
  df$gene <- rownames(df)
  df <- df[df$p_val<0.01,c("gene","avg_log2FC","pct.1")]
  #设置在细胞类型的表达FC的阈值
  df <- df[df$avg_log2FC>0,]
  df <- df[df$pct.1>0.25,]
  
  df <- head(df[order(df$avg_log2FC,decreasing = TRUE),],200)
  df$avg_log2FC <- 2
  
  df <- df[,c("gene","avg_log2FC")]
  
  names(df) <- c("gene",paste0(celltype,"_avg_log2FC_",sp))
  df <- merge(homolist,df,by.x = "marmosetGene",by.y = "gene")
  df <- df[,c("mouseGene",paste0(celltype,"_avg_log2FC_",sp))]
  df_a <- merge(df_a,df,by = "mouseGene",all = TRUE)
}

sp <- "Mice"
micefiles <- list.files(paste0(args[4],sp,"/"),pattern = ".csv")
for (mf in micefiles) {
  celltype <- sub("_positive.csv","",sub("mice_Annotated_Cerebellum_marker_","",mf))
  print(celltype)
  df <- read.csv(paste0(args[4],sp,"/",mf),header = TRUE,row.names = 1)
  df$gene <- rownames(df)
  df <- df[df$p_val<0.01,c("gene","avg_log2FC","pct.1")]
  #设置在细胞类型的表达FC的阈值
  df <- df[df$avg_log2FC>0,]
  df <- df[df$pct.1>0.25,]
  
  df <- head(df[order(df$avg_log2FC,decreasing = TRUE),],200)
  df$avg_log2FC <- 2
  
  df <- df[,c("gene","avg_log2FC")]
  
  names(df) <- c("gene",paste0(celltype,"_avg_log2FC_",sp))
  df <- merge(homolist,df,by.x = "mouseGene",by.y = "gene")
  df <- df[,c("mouseGene",paste0(celltype,"_avg_log2FC_",sp))]
  df_a <- merge(df_a,df,by = "mouseGene",all = TRUE)
}

df_a <- merge(homolist,df_a,by = "mouseGene")
df_b <- df_a[,c("mouseGene","macaqueGene","marmosetGene","Astrocyte_avg_log2FC_macaque","Astrocyte_avg_log2FC_marmoset","Astrocyte_avg_log2FC_Mice","Bergmann_avg_log2FC_macaque","Bergmann_avg_log2FC_marmoset","Bergmann_avg_log2FC_Mice","Choroid_avg_log2FC_macaque","Choroid_avg_log2FC_marmoset","Choroid_avg_log2FC_Mice","Endo_stalk_avg_log2FC_macaque","Endo_stalk_avg_log2FC_marmoset","Endo_stalk_avg_log2FC_Mice","Fibroblast_avg_log2FC_macaque","Fibroblast_avg_log2FC_marmoset","Fibroblast_avg_log2FC_Mice","Golgi_avg_log2FC_macaque","Golgi_avg_log2FC_marmoset","Golgi_avg_log2FC_Mice","Granule_avg_log2FC_macaque","Granule_avg_log2FC_marmoset","Granule_avg_log2FC_Mice","Microglia_avg_log2FC_macaque","Microglia_avg_log2FC_marmoset","Microglia_avg_log2FC_Mice","MLI1_avg_log2FC_macaque","MLI1_avg_log2FC_marmoset","MLI1_avg_log2FC_Mice","MLI2_avg_log2FC_macaque","MLI2_avg_log2FC_marmoset","MLI2_avg_log2FC_Mice","ODC_avg_log2FC_macaque","ODC_avg_log2FC_marmoset","ODC_avg_log2FC_Mice","OPC_avg_log2FC_macaque","OPC_avg_log2FC_marmoset","OPC_avg_log2FC_Mice","PLI_avg_log2FC_macaque","PLI_avg_log2FC_marmoset","PLI_avg_log2FC_Mice","Purkinje_avg_log2FC_macaque","Purkinje_avg_log2FC_marmoset","Purkinje_avg_log2FC_Mice","UBC_avg_log2FC_macaque","UBC_avg_log2FC_marmoset","UBC_avg_log2FC_Mice")]

result <- as.data.frame(matrix(nrow=0,ncol=3))
names(result) <- c("gene","celltype","type")
for (n in seq(4, 46, by=3)) {
  df_b[is.na(df_b)] <- as.integer(0)
  rownames(df_b) <- df_b$mouseGene
  tmp <- df_b[,(n):(n+2)]
  tmp$ans<-apply(tmp,1,function(x) sum(x>1))
  tmp_uniq <- tmp[tmp$ans==1,]
  tmp_comm <- tmp[tmp$ans==3,]
  
  
  macaque_uniq <- as.data.frame(rownames(tmp_uniq)[tmp_uniq[,1]!=0])
  if(length(macaque_uniq[,1])>0){
    ct1 <- strsplit(colnames(tmp_uniq)[1],split = "_")[[1]][1]
    print(ct1)
    names(macaque_uniq) <- "gene"
    macaque_uniq$celltype <- ct1
    macaque_uniq$type <- "macaque_uniq"
  }else{
    macaque_uniq <- as.data.frame(matrix(nrow=0,ncol=3))
    names(macaque_uniq) <- c("gene","celltype","type")
  }
  
  
  marmoset_uniq <- as.data.frame(rownames(tmp_uniq)[tmp_uniq[,2]!=0])
  if(length(marmoset_uniq[,1])>0){
    ct1 <- strsplit(colnames(tmp_uniq)[2],split = "_")[[1]][1]
    print(ct1)
    names(marmoset_uniq) <- "gene"
    marmoset_uniq$celltype <- ct1
    marmoset_uniq$type <- "marmoset_uniq"
  }else{
    marmoset_uniq <- as.data.frame(matrix(nrow=0,ncol=3))
    names(marmoset_uniq) <- c("gene","celltype","type")
  }
  
  mice_uniq <- as.data.frame(rownames(tmp_uniq)[tmp_uniq[,3]!=0])
  if (length(mice_uniq[,1])>0) {
    ct1 <- strsplit(colnames(tmp_uniq)[3],split = "_")[[1]][1]
    print(ct1)
    names(mice_uniq) <- "gene"
    mice_uniq$celltype <- ct1
    mice_uniq$type <- "mice_uniq"
  }else{
    mice_uniq <- as.data.frame(matrix(nrow=0,ncol=3))
    names(mice_uniq) <- c("gene","celltype","type")
  }
  
  
  common <- as.data.frame(rownames(tmp_comm))
  if(length(common[,1])>0) {
    ct1 <- strsplit(colnames(tmp_comm)[3],split = "_")[[1]][1]
    print(ct1)
    names(common) <- "gene"
    common$celltype <- ct1
    common$type <- "common"
  }else{
    common <- as.data.frame(matrix(nrow=0,ncol=3))
    names(common) <- c("gene","celltype","type")
  }
  result <- rbind(result,macaque_uniq,marmoset_uniq,mice_uniq,common)
}

genefra <- as.data.frame(table(result$gene))
head(genefra)

result_sub <- result[result$gene%in%genefra$Var1[genefra$Freq==1],]

write.csv(result_sub,args[5],quote = FALSE)
write.csv(result,args[6],quote = FALSE)

#qugene <- as.character(genefra$Var1[genefra$Freq!=1])
library(VennDiagram)
setwd(args[7])
for (n in seq(4, 46, by=3)) {
  df_b[is.na(df_b)] <- as.integer(0)
  rownames(df_b) <- df_b$mouseGene
  tmp <- df_b[,(n):(n+2)]
  #print(head(tmp))
  ct1 <- strsplit(colnames(tmp)[1],split = "_")[[1]][1]
  print(ct1)
  macaque_gene <- rownames(tmp)[tmp[,1]==2]
  marmoset_gene <- rownames(tmp)[tmp[,2]==2]
  mice_gene <- rownames(tmp)[tmp[,3]==2]
  #macaque_gene_uniq <- macaque_gene[!(macaque_gene%in%qugene)]
  #marmoset_gene_uniq <- marmoset_gene[!(marmoset_gene%in%qugene)]
  #mice_gene_uniq <- mice_gene[!(mice_gene%in%qugene)]
  #print(length(macaque_gene_uniq))
  #print(length(marmoset_gene_uniq))
  #print(length(mice_gene_uniq))
  venn.diagram(list(Mice=mice_gene,Macaque=macaque_gene,Marmoset=marmoset_gene),
               main=ct1,
               filename = paste0("venn0811/",ct1,"_VennDiagram.png"))
}