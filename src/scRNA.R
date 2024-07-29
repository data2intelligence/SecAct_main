

for(cancer in cancers)
{
  inputPath <- paste0("/data/rub2/data/single_cell_TF_activity/",cancer,".pickle.gz")
  outputPath <- paste0(validationPath,"singlecell/",cancer,".csv.gz")

  transformFormat <- paste0("python Secretome_s5_11_singlecell.py ",inputPath," ",outputPath)
  system(transformFormat)


  cdata_T <- as.matrix(read.csv(gzfile(paste0(validationPath,"singlecell/",cancer,".csv.gz")),row.names=1,check.names=F))
  rownames(cdata_T) <- transferSymbol(rownames(cdata_T))
  cdata_T <- rm_duplicates(cdata_T)


  cdata_T_minusBG <- cdata_T-rowMeans(cdata_T, na.rm=T)
  cdata_T_minusBG <- round(cdata_T_minusBG,3)


  library(SecAct)
  sigs <- c(0,112,36,618)
  sigNames <- c("CytoSig","NicheNet","ImmuneDic","SecAct")

  sigs <- c(618)
  sigNames <- c("SecAct")

  for(sig in sigs)
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

    res <- SecAct.inference(cdata_T_minusBG, SigMat=ref, lambda=a, nrand=1000)
    save(res, file = paste0(validationPath,"singlecell/",cancer,"_",sigNames[which(sigs==sig)],".RData"))
  }


}




for(cancer in cancers)
{

  sp_tf <- list(
    IFN.STAT1=list(c("IFN1","IFNG","IFNL","IFNL1","IL21","IL27") , "STAT1"),
    IL6and10fam.STAT3=list(c("IL6","LIF","OSM","IL10","IL22","IL9","IL21","IL27","IL17A","IL23A") , "STAT3"),
    IL12fam.STAT4=list(c("IL12","IL23A","IL27") , "STAT4"),
    IL2fam.STAT5=list(c("IL2","IL3","IL7","IL9","IL15","IL27") , "STAT5"),
    IL1andTNF.RELA=list(c("IL1A","IL1B","TNF") , "RELA"),
    TNFSF.RELB=list(c("LTA","CD40LG","TNFSF11","TNFSF12","TNFSF13B","TNFSF14") , "RELB"),
    BMP.SMAD1=list(c("BMP2","BMP3","BMP4","BMP5","BMP6","BMP7","BMP8A") , "SMAD1"),
    TGFB.SMAD2=list(c("INHBA","GDF11","TGFB1","TGFB2","TGFB3") , "SMAD2"),
    TGFB.SMAD3=list(c("INHBA","GDF11","TGFB1","TGFB2","TGFB3") , "SMAD3"),
    BMPandTGFB.SMAD4=list(c("BMP2","BMP3","BMP4","BMP5","BMP6","BMP7","BMP8A","INHBA","GDF11","TGFB1","TGFB2","TGFB3") , "SMAD4")
  )


  for(sigName in c("CytoSig","NicheNet","ImmuneDic","SecAct"))
  {

    summaryTable <- data.frame()

    load(paste0(validationPath,"singlecell/",cancer,"_",sigName,".RData"))
    sp <- as.matrix(res$zscore)

    if(sigName=="CytoSig")
    {
      rownames(sp) <- transferSymbol(rownames(sp))
      rownames(sp)[1] <- "INHBA"
    }

    tf <- as.matrix(read.table(paste0("/data/rub2/data/single_cell_TF_activity/",cancer,".Rabit.Cistrome.t.feature.gz"),sep="\t",check.names=F)	)
    tf <- t(tf)

    olp <- intersect(colnames(sp),colnames(tf))
    sp_olp <- sp[,olp]
    tf_olp <- tf[,olp]

    tfs <- sapply(strsplit(rownames(tf_olp),".",fixed=T),function(x) return(x[2]))
    celltypes <- sapply(strsplit(colnames(sp_olp),".",fixed=T),function(x) return(x[1]))
    celltypes <- gsub(",","",celltypes)

    rownames(tf_olp) <- tfs

    for(celltype in unique(celltypes))
    {
      sp_olp_sub <- sp_olp[,celltypes==celltype]
      tf_olp_sub <- tf_olp[,celltypes==celltype,drop=F]
      if(ncol(tf_olp_sub)<2) next

      tf_olp_sub <- tf_olp_sub[rowSums(tf_olp_sub)!=0,]

      for(i in names(sp_tf))
      {
        if(sum(rownames(sp_olp_sub)%in%sp_tf[[i]][[1]])==0) next

        sp_act <- colMeans(sp_olp_sub[rownames(sp_olp_sub)%in%sp_tf[[i]][[1]],,drop=F])

        if(!sp_tf[[i]][[2]]%in%rownames(tf_olp_sub)) next
        tf_act <- colMeans(tf_olp_sub[sp_tf[[i]][[2]],,drop=F])

        if(sum(tf_act>0)==0 | sum(tf_act<=0)==0) next

        fg.df <- data.frame(sp_act=sp_act,tf_act=tf_act)
        fg.df <- fg.df[fg.df[,2]!=0,]
        fg.df[,2] <- fg.df[,2]>0
        if(length(unique(fg.df[,2]))==1) next


        library(ROCR)
        predM <- prediction(fg.df$sp_act, fg.df$tf_act)
        roc = performance(predM, measure = "auc")

        summaryTable[paste0(cancer,celltype,i),"cancer"] <- cancer
        summaryTable[paste0(cancer,celltype,i),"celltype"] <- celltype
        summaryTable[paste0(cancer,celltype,i),"SP_TF"] <- i
        summaryTable[paste0(cancer,celltype,i),"AUC"] <- roc@ y.values
      }
    }

    write.csv(summaryTable,paste0(validationPath,"singlecell/",cancer,"_",sigName,"_SP_TF.csv"),quote=F)

  } # sigName

}



summaryTable.median.median <- data.frame() # sig ~ cancer
summaryTable.median.median2 <- data.frame() # sig ~ SPTF

for(sigName in c("CytoSig","NicheNet","ImmuneDic","SecAct"))
{
  cancers <- list.files(paste0(validationPath,"singlecell/"))
  cancers <- cancers[grepl(paste0(sigName,"_SP_TF.csv"),cancers)]
  cancers <- gsub(paste0("_",sigName,"_SP_TF.csv"),"",cancers)

  summaryTable <- data.frame()
  for(cancer in cancers)
  {
    temp <- read.csv(paste0(validationPath,"singlecell/",cancer,"_",sigName,"_SP_TF.csv"),row.names=1,header=T)
    summaryTable <- rbind(summaryTable,temp)
  }

  summaryTable.median <- data.frame()
  for(cancer in cancers)
  {
    summaryTableCancer <- summaryTable[summaryTable[,"cancer"]==cancer,]
    for(i in unique(summaryTableCancer[,"SP_TF"]))
    {
      summaryTable_sub <- summaryTableCancer[summaryTableCancer[,"SP_TF"]==i,]

      summaryTable.median[paste0(cancer,i),"cancer"] <- cancer
      summaryTable.median[paste0(cancer,i),"SP_TF"] <- i
      summaryTable.median[paste0(cancer,i),"AUC.median"] <- median(summaryTable_sub[,4])
    }
  }

  for(cancer in cancers)
  {
    summaryTable.median_sub <- summaryTable.median[summaryTable.median[,"cancer"]==cancer,]

    summaryTable.median.median[paste0(sigName,cancer),"sig"] <- sigName
    summaryTable.median.median[paste0(sigName,cancer),"cancer"] <- cancer
    summaryTable.median.median[paste0(sigName,cancer),"AUC.median.median"] <- median(summaryTable.median_sub[,"AUC.median"])
  }

  for(SPTF in unique(summaryTable.median[,"SP_TF"]))
  {
    summaryTable.median_sub <- summaryTable.median[summaryTable.median[,"SP_TF"]==SPTF,]

    summaryTable.median.median2[paste0(sigName,SPTF),"sig"] <- sigName
    summaryTable.median.median2[paste0(sigName,SPTF),"SPTF"] <- SPTF
    summaryTable.median.median2[paste0(sigName,SPTF),"AUC.median.median"] <- median(summaryTable.median_sub[,"AUC.median"])
  }

  #write.csv(summaryTable.median,paste0(validationPath,"validation_singlecell.csv"),quote=F)




  sigOrder <- data.frame()
  for(i in unique(summaryTable.median[,"cancer"]))
  {
    sigOrder[i,"r_mean"] <- median(summaryTable.median[summaryTable.median[,"cancer"]==i,"AUC.median"],na.rm=T)
  }
  sigOrder <- sigOrder[order(sigOrder[,"r_mean"],decreasing=T),,drop=F]

  summaryTable.median[,"cancer"] <- factor(summaryTable.median[,"cancer"],levels=rownames(sigOrder))

  library(ggplot2)
  p1 <- ggplot(summaryTable.median, aes(x = cancer, y = AUC.median))+
    geom_hline(yintercept=0.5, color = "grey1", linewidth=0.8, linetype = 'dotted')+
    geom_boxplot(width=0.6, outlier.shape=NA, fill="#fde0dd", alpha=0.5)+
    geom_jitter(aes(colour=SP_TF),size=0.8)+
    #scale_colour_manual(name="SP.TF", values = c(RColorBrewer::brewer.pal(9,"Set1"),"black","green","cyan","blueviolet") )+
    scale_colour_manual(name="SP.TF", values = c(RColorBrewer::brewer.pal(10,"Paired")) )+
    #ggtitle("SP-TF relation in scRNA-seq data")+
    guides(colour = guide_legend(override.aes = list(size=2)))+
    coord_cartesian(ylim = c(0.4, 1))+
    ylab("AUC")+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(size=10,colour = "black"),
      axis.title = element_text(size=10,colour = "black"),
      axis.text.x = element_text(size=8, angle = 90, hjust = 1, vjust = 0.5),
      axis.title.x = element_blank(),
      legend.title = element_blank(),
      legend.position="bottom"
    )

  ggsave(paste0(validationPath,"validation_singlecell_",sigName,".png"), p1, width = 16, height = 13.5, dpi=200, units = "cm",limitsize = FALSE)



  mat = reshape2::dcast( summaryTable.median , cancer~SP_TF )
  rownames(mat) <- mat[,1]
  mat <- t(mat[,-1])

  mat <- mat[order(apply(mat,1,function(x) median(x,na.rm=T)),decreasing=T),]
  #mat <- mat[,order(apply(mat,2,function(x) median(x,na.rm=T)))]
  mat <- mat[,sort(colnames(mat))]

  mat <- t(mat)
  colnames(mat) <- gsub("and","&", colnames(mat))

  widthValue <- ifelse(sigName=="ImmuneDic",19.9,21)

  png(paste0(validationPath,"validation_singlecell_",sigName,"_heatmap.png"), width = widthValue, height = 20, res=200, units = "cm")

  library(ComplexHeatmap)

  technique_vec <- sapply(strsplit(rownames(mat),".",fixed=T),function(x) return(x[3]))
  dataset_vec <- sapply(strsplit(rownames(mat),".",fixed=T),function(x) return(paste0(x[1],".",x[2],"          ")))

  #row_ha <- rowAnnotation(
  #	AUC = anno_boxplot(as.matrix(mat), height = unit(3, "cm") )
  #)
  row_ha <- rowAnnotation(
    Seq = technique_vec,
    col = list(Seq = c("InDrop" = "#B395BD", "10x" = "#7DAEE0", "SS2" = "#EA8379")),
    labels = anno_text(dataset_vec, which = "row", gp = gpar(fontsize = 10) ),
    width = max(grobWidth(textGrob(dataset_vec)))
  )
  column_ha <- columnAnnotation(
    AUC = anno_boxplot(as.matrix(mat), height = unit(2.8, "cm"), gp = gpar(fill = "skyblue"))
  )

  ht <- Heatmap(as.matrix(mat),
                name = "AUC",
                rect_gp = gpar(col = "grey10", lwd = 2),
                row_title = " \n \n \n \n ",
                column_title = sigName,
                col = circlize::colorRamp2(c(0,0.5,1), c("#aad962", "white", "#ef6a32")),
                row_names_max_width = max_text_width(
                  rownames(mat),
                  gp = gpar(fontsize = 11)
                ),
                column_names_max_height = max_text_width(
                  colnames(mat),
                  gp = gpar(fontsize = 12)
                ),
                show_row_names = FALSE,
                show_column_names = TRUE,
                top_annotation = column_ha,
                right_annotation = row_ha,
                column_names_rot = 38,
                cluster_rows = FALSE,
                cluster_columns = FALSE

  )

  draw(ht)

  # y = 0.5
  decorate_annotation("AUC", {
    grid.lines(c(.5, 10.5), c(0.5, 0.5), gp = gpar(lty = 2, col = "grey20"),
               default.units = "native")
  })

  dev.off()

}




sigOrder <- data.frame()
for(i in unique(summaryTable.median.median[,"sig"]))
{
  sigOrder[i,"r_mean"] <- median(summaryTable.median.median[summaryTable.median.median[,"sig"]==i,"AUC.median.median"],na.rm=T)
}
sigOrder <- sigOrder[order(sigOrder[,"r_mean"],decreasing=T),,drop=F]

summaryTable.median.median[,1] <- factor(summaryTable.median.median[,1], levels=rownames(sigOrder))

p2 <- ggplot(summaryTable.median.median, aes(x = sig, y = AUC.median.median, fill=sig))+
  geom_hline(yintercept=0.5, color = "grey1", linewidth=0.8, linetype = 'dotted')+
  geom_boxplot(width=0.7,outlier.shape=NA)+
  scale_fill_manual( values=sigColors[rownames(sigOrder)] )+
  geom_jitter(size=0.1,alpha=0.2)+
  ggtitle(" ")+
  xlab(" ")+
  ylab("AUC")+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size=10,colour = "black"),
    axis.title = element_text(size=10,colour = "black"),
    axis.text.x = element_text(size=10,angle = 45, hjust = 1),
    legend.position="none"
  )

ggsave(paste0(validationPath,"validation_singlecell_comb_dataset.png"), p2, width = 6.2, height = 8, dpi=200, units = "cm",limitsize = FALSE)


sigOrder <- data.frame()
for(i in unique(summaryTable.median.median2[,"sig"]))
{
  sigOrder[i,"r_mean"] <- median(summaryTable.median.median2[summaryTable.median.median2[,"sig"]==i,"AUC.median.median"],na.rm=T)
}
sigOrder <- sigOrder[order(sigOrder[,"r_mean"],decreasing=T),,drop=F]

summaryTable.median.median2[,1] <- factor(summaryTable.median.median2[,1], levels=rownames(sigOrder))

p2 <- ggplot(summaryTable.median.median2, aes(x = sig, y = AUC.median.median, fill=sig))+
  geom_hline(yintercept=0.5, color = "grey1", linewidth=0.8, linetype = 'dotted')+
  geom_boxplot(width=0.7,outlier.shape=NA)+
  scale_fill_manual( values=sigColors[rownames(sigOrder)] )+
  geom_jitter(size=0.5,alpha=0.5,color="red",width=0.1)+
  ggtitle(" ")+
  xlab(" ")+
  ylab("AUC")+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size=10,colour = "black"),
    axis.title = element_text(size=10,colour = "black"),
    axis.text.x = element_text(size=10,angle = 45, hjust = 1),
    legend.position="none"
  )

ggsave(paste0(validationPath,"validation_singlecell_comb_SPTF.png"), p2, width = 6.2, height = 8, dpi=400, units = "cm",limitsize = FALSE)




mat = reshape2::dcast( summaryTable.median.median2 , sig ~ SPTF )
rownames(mat) <- mat[,1]
mat <- t(mat[,-1])

mat <- mat[order(apply(mat,1,function(x) median(x,na.rm=T))),]
mat <- mat[,order(apply(mat,2,function(x) median(x,na.rm=T)))]


png(paste0(validationPath,"validation_singlecell_comb_SPTF_heatmap.png"), width = 18, height = 20, res=200, units = "cm")

library(ComplexHeatmap)

row_ha <- rowAnnotation(
  z = anno_boxplot(as.matrix(mat), height = unit(3, "cm") )
)
column_ha <- columnAnnotation(
  z = anno_boxplot(as.matrix(mat), height = unit(3, "cm"), gp = gpar(fill = sigColors[colnames(mat)]) )
)

ht <- Heatmap(as.matrix(mat),
              name = "AUC",
              rect_gp = gpar(col = "white", lwd = 2),
              col = circlize::colorRamp2(c(-1, 0,1), c("#91bfdb", "white", "#fc8d59")),
              row_names_max_width = max_text_width(
                rownames(mat),
                gp = gpar(fontsize = 12)
              ),
              column_names_max_height = max_text_width(
                colnames(mat),
                gp = gpar(fontsize = 12)
              ),
              column_names_rot = 38,
              top_annotation = column_ha,
              #right_annotation = row_ha,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              cell_fun = function(j, i, x, y, width, height, fill) {grid.text(sprintf("%.3f", mat[i, j]), x, y, gp = gpar(fontsize = 10))}
)
draw(ht)

dev.off()
