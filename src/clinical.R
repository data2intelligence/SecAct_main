
items <- c(
  "TGFB_GSE174686","TGFB3_GSE174686",
  "IL1B_GSE80060","IL1B_GSE57253","IL1B_GSE70019",
  "IFNG_GSE100093","IFNG_GSE78193",
  "TNFSF12_GSE42049","TNFSF12_GSE42048", # tumor
  "IL6_GSE62941","IL6_GSE45867","IL6_GSE61201", # tumor x x
  "TNF_GSE48498","TNF_GSE11903",
  "IL1A_GSE57253","IL1A_GSE70019",
  "NTN1_GSE225691", # tumor
  "IL11_Widjaja2024"
)

for(item in items)
{
  gene <- strsplit(item,"_") [[1]][1]
  dataset <- strsplit(item,"_") [[1]][2]

  if(item%in%c("TGFB_GSE174686","TGFB3_GSE174686"))
  {
    dataPath <- "../../SpaCE/data/"
    f3 <- as.matrix(read.csv(gzfile(paste0(dataPath,"aTGFB_XOMA.diff.human.gz")),as.is=T,sep="\t",check.names=F))

    if(item%in%c("TGFB_GSE174686"))
    {
      f3 <- f3[,"mouse.Pan",drop=F]
    }
    if(item=="TGFB3_GSE174686")
    {
      f3 <- f3[,"Pan-1+2",drop=F]
    }
  }else if(item%in%c("VEGFA_GSE72951","VEGFA_E-MTAB-3267")){
    dataPath <- "../../SpaCE/code/CytoSig_prediction-master/data/bulk/tumor/"

    if(item=="VEGFA_GSE72951")
    {
      f3 <- as.matrix(read.csv(gzfile(paste0(dataPath,"GSE72951.self_subtract.gz")),as.is=T,sep="\t"))
    }else if(item=="VEGFA_E-MTAB-3267"){
      f3 <- as.matrix(read.csv(gzfile(paste0(dataPath,"E-MTAB-3267.norm_subtract.gz")),as.is=T,sep="\t",check.names=F))
    }

  }else if(item%in%c("NTN1_GSE225691")){

    NTN_bulk <- as.matrix(read.csv("/data/rub2/project/Secretome/data/GSE225691_NTN1/GSE225687_Patient_Raw_Count.csv",sep=";",row.names=1,header=T))
    rownames(NTN_bulk) <- substr(rownames(NTN_bulk),1,15)
    NTN_bulk <- rm_duplicates(NTN_bulk)

    gene_anno <- as.matrix(read.csv(gzfile(paste0("../_raw/BRCA_10x_Datasets/Version1.0.0_Breast.Cancer_rep1/filtered_feature_bc_matrix/features.tsv.gz")),as.is=T,header=F,sep="\t"))

    olp <- intersect(rownames(NTN_bulk),gene_anno[,1])

    NTN_bulk <- NTN_bulk[olp,]
    rownames(NTN_bulk) <- gene_anno[match(rownames(NTN_bulk),gene_anno[,1]),2]
    rownames(NTN_bulk) <- transferSymbol(rownames(NTN_bulk))
    NTN_bulk <- rm_duplicates(NTN_bulk)

    expr <- NTN_bulk

    expr.scaled <- t(t(expr)*1e6/colSums(expr))
    expr.scaled.log <- log2(expr.scaled + 1 )

    f3 <- expr.scaled.log[,(1:12)*2] - expr.scaled.log[,(1:12)*2 - 1]

  }else if(item%in%c("IL11_Widjaja2024")){

    library(rio)
    xlsx <- import_list(paste0(validationPath,"/IL11/IL11_Widjaja2024_TableS5.xlsx"))

    s0 <- as.data.frame(xlsx[[1]])
    s0 <- s0[s0[,2]>0,]
    s0[,3] <- as.numeric(s0[,3])
    s0 <- s0[order(s0[,2],decreasing=T),]
    s0 <- s0[!duplicated(s0[,8]),]
    rownames(s0) <- s0[,8]
    s1 <- s0[,3,drop=F]

    s0 <- as.data.frame(xlsx[[2]])
    s0 <- s0[s0[,2]>0,]
    s0[,3] <- as.numeric(s0[,3])
    s0 <- s0[order(s0[,2],decreasing=T),]
    s0 <- s0[!duplicated(s0[,8]),]
    rownames(s0) <- s0[,8]
    s2 <- s0[,3,drop=F]

    s0 <- as.data.frame(xlsx[[3]])
    s0 <- s0[s0[,2]>0,]
    s0[,3] <- as.numeric(s0[,3])
    s0 <- s0[order(s0[,2],decreasing=T),]
    s0 <- s0[!duplicated(s0[,8]),]
    rownames(s0) <- s0[,8]
    s3 <- s0[,3,drop=F]

    olp <- intersect(rownames(s1), intersect(rownames(s2),rownames(s3)))

    f3 <- cbind(vWat=s1[olp,], cbind(liver=s2[olp,], gastro=s3[olp,]) )
    rownames(f3) <- olp

    f3 <- cbind(f3, meanLogFC=apply(f3,1,mean))

    f3 <- transferMouseToHuman(f3)

    f3 <- f3[,4,drop=F]

  }else{
    dataPath <- "../../SpaCE/code/CytoSig_prediction-master/data/bulk/inflam/"
    allFiles <- list.files(dataPath)

    flag1 <- grepl(dataset,allFiles)
    flag2 <- grepl(".diff.1.cntmap",allFiles)

    fileName_full <- allFiles[flag1&flag2]
    fileName_full <- gsub(".diff.1.cntmap","",fileName_full)

    f3 <- as.matrix(read.csv(gzfile(paste0(dataPath,fileName_full,".diff.1")),as.is=T,sep="\t"))

  }

  rownames(f3) <- transferSymbol(rownames(f3))
  f3 <- rm_duplicates(f3)

  library(SecAct)

  for(sig in c(0,112,36,618))
  {
    if(sig==0)
    {
      ref <- paste0("/data/rub2/project/Secretome/code/CytoSig-master/CytoSig/signature.centroid")
    }else if(sig==618){
      ref <- paste0("/data/rub2/project/Secretome/results/signatureComb/AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3.tsv")
    }else if(sig==6182){
      ref <- paste0("/data/rub2/project/Secretome/results/signatureComb/AllSigFilteredBy_MoranI_TCGA_0.25_ds3_ICGC_highConfidence.tsv")
    }else{
      ref <- paste0("/data/rub2/project/Secretome/code/CytoSig-master/CytoSig/signature.centroid.",sig)
    }

    aMat <- read.csv(paste0(QCSummaryPath,"lambda_summary/lambda_compare_",sig,".csv"))
    a <- aMat[aMat[,2]==max(aMat[,2]),1]

    res <- SecAct.inference(f3, SigMat=ref, lambda=a, nrand=1000)
    save(res, file = paste0(validationPath,"clinical/",item,"_",sig,"_",a,".RData"))


    print(paste0(a,"_",sig))

  }
}







compSigNames <- c(
  "CytoSig","NicheNet",
  "ImmuneDic","SecAct"
)
CompSigAlphas <- c(
  1e+07,1e+07,
  50000,5e+05
)


for(compSig in compSigs)
{
  clinical <- data.frame()

  for(item in items)
  {
    gene <- strsplit(item,"_") [[1]][1]
    dataset <- strsplit(item,"_") [[1]][2]

    #for(a in c(10000,50000,1e+05,5e+05,1e+06,5e+06,1e+07,5e+07) )
    #{
    #if(!file.exists(paste0(validationPath,"clinical/",item,"_",compSig,"_",a,".RData"))) next
    a <- CompSigAlphas[compSigs==compSig]

    load(paste0(validationPath,"clinical/",item,"_",compSig,"_",a,".RData"))
    f3 <- as.matrix(res$zscore)

    if(compSig==0)
    {
      rownames(f3) <- transferSymbol(rownames(f3))
      f3 <- rm_duplicates(f3)
    }

    f3 <- rbind(f3,TGFB=colMeans(f3[rownames(f3)%in%c("TGFB1","TGFB2","TGFB3"),,drop=F]))

    f3 <- f3[,!grepl("Placebo",colnames(f3),fixed=T) & !is.null(colnames(f3)),drop=F]
    f3 <- f3[,!grepl("non.responder",colnames(f3),fixed=T) & !is.null(colnames(f3)),drop=F]

    if(item=="IFN1_GSE72754")
    {
      f3 <- f3[,"IFNA.high.whole.blood",drop=F]
    }
    if(item=="IL1B_GSE80060")
    {
      f3 <- f3[,c("IL1B.Day3.Canakinumab_100.0","IL1B.Day3.Canakinumab_90.0","IL1B.Day3.Canakinumab_70.0"),drop=F]
    }
    if(item=="TNFSF12_GSE42048")
    {
      f3 <- f3[,c("TWEAK.ACHN.xenograft.72h","TWEAK.ACHN.xenograft.24h"),drop=F]
    }
    if(item=="TNF_GSE48498")
    {
      f3 <- f3[,"TNFA.whole.blood.cells",drop=F]
    }

    f3 <- apply(f3,1,function(x) mean(x,na.rm=T))

    clinical[paste0(item,"_",compSig,"_",a),"item"] <- item
    clinical[paste0(item,"_",compSig,"_",a),"compSig"] <- paste0(compSig,"_",a)
    clinical[paste0(item,"_",compSig,"_",a),"value"] <- f3[gene]
    #} #a
  }


  mat = reshape2::dcast( clinical , compSig~item )
  rownames(mat) <- mat[,1]
  mat <- t(mat[,-1])

  mat <- mat[order(apply(mat,1,function(x) median(x,na.rm=T))),,drop=F]
  mat <- mat[,order(apply(mat,2,function(x) median(x,na.rm=T))),drop=F]


  library(ComplexHeatmap)

  png(paste0(validationPath,"validation_clinical_blockade_",compSigNames[which(compSigs==compSig)],".png"), width = 20, height = 20, res=200, units = "cm")

  row_ha <- rowAnnotation(
    z = anno_boxplot(as.matrix(mat), height = unit(3, "cm") )
  )
  column_ha <- columnAnnotation(
    z = anno_boxplot(as.matrix(mat), height = unit(3, "cm") )
  )

  ht <- Heatmap(as.matrix(mat),
                column_title=compSigNames[which(compSigs==compSig)],
                name = "Risk z",
                col = circlize::colorRamp2(c(-3, 0,3), c("green", "white", "red")),
                row_names_max_width = max_text_width(
                  rownames(mat),
                  gp = gpar(fontsize = 12)
                ),
                column_names_max_height = max_text_width(
                  colnames(mat),
                  gp = gpar(fontsize = 12)
                ),
                top_annotation = column_ha,
                right_annotation = row_ha,
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                cell_fun = function(j, i, x, y, width, height, fill) {grid.text(sprintf("%.3f", mat[i, j]), x, y, gp = gpar(fontsize = 10))}
  )

  draw(ht)
  dev.off()

  #clinicalBest <- clinical[clinical[,2]==colnames(mat)[1],,drop=F]
  clinicalBest <- clinical

  if(compSig==0)
  {
    clinicalComb <- clinicalBest
  }else{
    clinicalComb <- rbind(clinicalComb,clinicalBest)
  }
}


for(compSig in compSigs)
{
  clinicalComb[clinicalComb[,2]==paste0(compSig,"_",CompSigAlphas[compSigs==compSig]),2] <- compSigNames[compSigs==compSig]
}


source("createLRdb.R")
LRdb <- rbind(LRdb,c("NTN1","DCC"))
LRdb <- rbind(LRdb,c("NTN1","UNC5C"))

for(item in items)
{
  gene <- strsplit(item,"_") [[1]][1]
  dataset <- strsplit(item,"_") [[1]][2]

  if(item%in%c("TGFB_GSE174686","TGFB3_GSE174686"))
  {
    f3 <- as.matrix(read.csv(gzfile("../../SpaCE/data/aTGFB_XOMA.diff.human.gz"),as.is=T,sep="\t",check.names=F))
    rownames(f3) <- transferSymbol(rownames(f3))
    f3 <- rm_duplicates(f3)

    if(item%in%c("TGFB_GSE174686"))
    {
      f3 <- f3[,"mouse.Pan",drop=F]
      f3 <- rbind(f3,TGFB=colMeans(f3[rownames(f3)%in%c("TGFB1","TGFB2","TGFB3"),,drop=F]))
    }
    if(item=="TGFB3_GSE174686")
    {
      f3 <- f3[,"Pan-1+2",drop=F]
    }

  }else if(item%in%c("NTN1_GSE225691")){
    gene_anno <- as.matrix(read.csv(gzfile(paste0("../_raw/BRCA_10x_Datasets/Version1.0.0_Breast.Cancer_rep1/filtered_feature_bc_matrix/features.tsv.gz")),as.is=T,header=F,sep="\t"))

    NTN_bulk <- as.matrix(read.csv("/data/rub2/project/Secretome/data/GSE225691_NTN1/GSE225687_Patient_Raw_Count.csv",sep=";",row.names=1,header=T))
    rownames(NTN_bulk) <- substr(rownames(NTN_bulk),1,15)
    NTN_bulk <- rm_duplicates(NTN_bulk)

    olp <- intersect(rownames(NTN_bulk),gene_anno[,1])

    NTN_bulk <- NTN_bulk[olp,]
    rownames(NTN_bulk) <- gene_anno[match(rownames(NTN_bulk),gene_anno[,1]),2]
    rownames(NTN_bulk) <- transferSymbol(rownames(NTN_bulk))
    NTN_bulk <- rm_duplicates(NTN_bulk)

    expr <- NTN_bulk

    expr.scaled <- t(t(expr)*1e6/colSums(expr))
    expr.scaled.log <- log2(expr.scaled + 1 )


    f3 <- expr.scaled.log[,(1:12)*2] - expr.scaled.log[,(1:12)*2 - 1]

  }else if(item%in%c("IL11_Widjaja2024")){

    library(rio)
    xlsx <- import_list(paste0(validationPath,"/IL11/IL11_Widjaja2024_TableS5.xlsx"))

    s0 <- as.data.frame(xlsx[[1]])
    s0 <- s0[s0[,2]>0,]
    s0[,3] <- as.numeric(s0[,3])
    s0 <- s0[order(s0[,2],decreasing=T),]
    s0 <- s0[!duplicated(s0[,8]),]
    rownames(s0) <- s0[,8]
    s1 <- s0[,3,drop=F]

    s0 <- as.data.frame(xlsx[[2]])
    s0 <- s0[s0[,2]>0,]
    s0[,3] <- as.numeric(s0[,3])
    s0 <- s0[order(s0[,2],decreasing=T),]
    s0 <- s0[!duplicated(s0[,8]),]
    rownames(s0) <- s0[,8]
    s2 <- s0[,3,drop=F]

    s0 <- as.data.frame(xlsx[[3]])
    s0 <- s0[s0[,2]>0,]
    s0[,3] <- as.numeric(s0[,3])
    s0 <- s0[order(s0[,2],decreasing=T),]
    s0 <- s0[!duplicated(s0[,8]),]
    rownames(s0) <- s0[,8]
    s3 <- s0[,3,drop=F]

    olp <- intersect(rownames(s1), intersect(rownames(s2),rownames(s3)))

    f3 <- cbind(s1[olp,], cbind(s2[olp,], s3[olp,]) )
    rownames(f3) <- olp

    f3 <- cbind(f3, meanLogFC=apply(f3,1,mean))

    f3 <- transferMouseToHuman(f3)

    f3 <- f3[,4,drop=F]

    rownames(f3) <- transferSymbol(rownames(f3))
    f3 <- rm_duplicates(f3)

  }else{
    dataPath <- "../../SpaCE/code/CytoSig_prediction-master/data/bulk/inflam/"
    allFiles <- list.files(dataPath)

    flag1 <- grepl(dataset,allFiles)
    flag2 <- grepl(".diff.1.cntmap",allFiles)

    fileName_full <- allFiles[flag1&flag2]
    fileName_full <- gsub(".diff.1.cntmap","",fileName_full)

    f3 <- as.matrix(read.csv(gzfile(paste0(dataPath,fileName_full,".diff.1")),as.is=T,sep="\t"))
    rownames(f3) <- transferSymbol(rownames(f3))
    f3 <- rm_duplicates(f3)
  }


  f3 <- f3[,!grepl("Placebo",colnames(f3),fixed=T),drop=F]
  f3 <- f3[,!grepl("non.responder",colnames(f3),fixed=T),drop=F]

  if(item=="IFN1_GSE72754")
  {
    f3 <- f3[,"IFNA.high.whole.blood",drop=F]
  }
  if(item=="IL1B_GSE80060")
  {
    f3 <- f3[,c("IL1B.Day3.Canakinumab_100.0","IL1B.Day3.Canakinumab_90.0","IL1B.Day3.Canakinumab_70.0"),drop=F]
  }
  if(item=="TNFSF12_GSE42048")
  {
    f3 <- f3[,c("TWEAK.ACHN.xenograft.4h","TWEAK.ACHN.xenograft.8h","TWEAK.ACHN.xenograft.72h","TWEAK.ACHN.xenograft.24h"),drop=F]
  }
  if(item=="TNF_GSE48498")
  {
    f3 <- f3[,"TNFA.whole.blood.cells",drop=F]
  }

  f3 <- apply(f3,1,function(x) mean(x,na.rm=T))
  f3 <- scale(f3)[,1]

  clinicalComb[paste0(item,"LigandExp"),"item"] <- item
  clinicalComb[paste0(item,"LigandExp"),"compSig"] <- "LigandExp"
  clinicalComb[paste0(item,"LigandExp"),"value"] <- f3[gene]


  if(gene=="TGFB") gene <- c("TGFB1","TGFB2","TGFB3")
  clinicalComb[paste0(item,"ReceptorExp"),"item"] <- item
  clinicalComb[paste0(item,"ReceptorExp"),"compSig"] <- "ReceptorExp"
  clinicalComb[paste0(item,"ReceptorExp"),"value"] <- mean(f3[unique(LRdb[LRdb[,1]%in%gene,2])], na.rm=TRUE)

  clinicalComb[paste0(item,"LRsumExp"),"item"] <- item
  clinicalComb[paste0(item,"LRsumExp"),"compSig"] <- "LRsumExp"
  clinicalComb[paste0(item,"LRsumExp"),"value"] <- mean(f3[c(gene,unique(LRdb[LRdb[,1]%in%gene,2]))], na.rm=TRUE)
}







mat = reshape2::dcast( clinicalComb , compSig~item )
rownames(mat) <- mat[,1]
mat <- t(mat[,-1])
mat <- mat[!is.na(mat[,"SecAct"]),]

mat <- mat[order(apply(mat,1,function(x) median(x,na.rm=T))),]
mat <- mat[,order(apply(mat,2,function(x) median(x,na.rm=T)))]


TNFSF12_GSE42048 <- read.csv("/data/rub2/project/Secretome/data/GSE42048.HG-U133_Plus_2.rma",sep="\t",row.names=1)
TNFSF12_GSE42049 <- read.csv("/data/rub2/project/Secretome/data/GSE42049.HG-U133_Plus_2.rma",sep="\t",row.names=1)

TNFSF12_GSE42048_logFC <- rowMeans(TNFSF12_GSE42048[,22:41])-rowMeans(TNFSF12_GSE42048[,42:62])
TNFSF12_GSE42049_logFC <- rowMeans(TNFSF12_GSE42049[,1:4])-rowMeans(TNFSF12_GSE42049[,5:8])

TNFSF12_GSE42048_logFC <- scale(TNFSF12_GSE42048_logFC)
TNFSF12_GSE42049_logFC <- scale(TNFSF12_GSE42049_logFC)

mat["TNFSF12_GSE42048","LigandExp"] <- TNFSF12_GSE42048_logFC["205611_at@TNFSF12",1]
mat["TNFSF12_GSE42048","ReceptorExp"]<- TNFSF12_GSE42048_logFC["218368_s_at@TNFRSF12A",1]
mat["TNFSF12_GSE42048","LRsumExp"] <- mean(TNFSF12_GSE42048_logFC[c("205611_at@TNFSF12","218368_s_at@TNFRSF12A"),1])

mat["TNFSF12_GSE42049","LigandExp"] <- TNFSF12_GSE42049_logFC["205611_at@TNFSF12",1]
mat["TNFSF12_GSE42049","ReceptorExp"] <- TNFSF12_GSE42049_logFC["205611_at@TNFSF12",1]
mat["TNFSF12_GSE42049","LRsumExp"] <- mean(TNFSF12_GSE42049_logFC[c("205611_at@TNFSF12","218368_s_at@TNFRSF12A"),1])



write.csv(mat,paste0(validationPath,"validation_clinical_final.csv"),quote=F)

meta <- read.csv(paste0(validationPath,"cytokine_disease_GEOID.csv"),as.is=T,row.names=1)
meta_sorted <- meta[rownames(mat),]
write.csv(meta_sorted,paste0(validationPath,"cytokine_disease_GEOID_sorted.csv"))


rownames(mat)[rownames(mat)=="TGFB_GSE174686"] <- "TGFB_GSE174686^"
rownames(mat)[rownames(mat)=="TGFB3_GSE174686"] <- "TGFB3_GSE174686^"
rownames(mat)[rownames(mat)=="TNFSF12_GSE42048"] <- "TNFSF12_GSE42048^"
rownames(mat)[rownames(mat)=="TNFSF12_GSE42049"] <- "TNFSF12_GSE42049^"
rownames(mat)[rownames(mat)=="IL6_GSE62941"] <- "IL6_GSE62941^"
rownames(mat)[rownames(mat)=="IL11_Widjaja2024"] <- "IL11_Widjaja2024^"




png(paste0(validationPath,"validation_clinical_final.png"), width = 17.5, height = 19, res=200, units = "cm")

library(ComplexHeatmap)
disease_vec <- meta_sorted[,"Disease"]
dataset_vec <- rownames(mat)

#row_ha <- rowAnnotation(
#	z = anno_boxplot(as.matrix(mat), height = unit(3, "cm") )
#)
row_ha <- rowAnnotation(
  Disease = disease_vec,
  col = list(Disease = c( "Cancer" = "#88cdbc", "Inflammatory" = "#e0d579", "Aging"="#eca680")),
  labels = anno_text(dataset_vec, which = "row", gp = gpar(fontsize = 10) ),
  width = max(grobWidth(textGrob(dataset_vec)))
)

column_ha <- columnAnnotation(
  z = anno_boxplot(as.matrix(mat), height = unit(3, "cm"), gp = gpar(fill = sigColors[colnames(mat)]) )
  #z = anno_boxplot(as.matrix(mat), height = unit(3, "cm") )
)

ht <- Heatmap(as.matrix(mat),
              name = "Activity Diff",
              rect_gp = gpar(col = "white", lwd = 2),
              col = circlize::colorRamp2(c(-5, 0,5), c("#91bfdb", "white", "#fc8d59")),
              row_names_max_width = max_text_width(
                rownames(mat),
                gp = gpar(fontsize = 12)
              ),
              column_names_max_height = max_text_width(
                colnames(mat),
                gp = gpar(fontsize = 12)
              ),
              show_row_names = FALSE,
              show_column_names = TRUE,
              column_names_rot = 38,
              top_annotation = column_ha,
              right_annotation = row_ha,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              cell_fun = function(j, i, x, y, width, height, fill) {grid.text(sprintf("%.3f", mat[i, j]), x, y, gp = gpar(fontsize = 10))}
)
draw(ht)

dev.off()
