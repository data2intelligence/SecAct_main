

# CPTAC
dataPath <- "/data/rub2/data/CPTAC/"
cancers <- c("BRCA", "COAD", "GBM", "HNSC", "KIRC", "LUAD", "LUSC", "OV", "PDAC", "UCEC")

geneAnno <- read.csv(paste0(dataPath,"Proteome_BCM_GENCODE_v34_harmonized_v1/README/Gene_annotation_and_representable_isoform_mapping_table.txt"),sep="\t",check.names=F)
geneAnno[geneAnno[,1]=="ENSG00000284024.2",4] <- "MSANTD7"

geneAnno <- geneAnno[,c(1,4)]
geneAnno <- geneAnno[!duplicated(geneAnno[,1]),]
rownames(geneAnno) <- geneAnno[,1]


for(cancer in cancers)
{
  cdata_T <- as.matrix(read.csv(paste0(dataPath,"RNA_BCM_v1/",cancer,"_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt"),sep="\t",row.names=1,check.names=F))

  rownames(cdata_T) <- geneAnno[rownames(cdata_T),2]
  rownames(cdata_T) <- transferSymbol(rownames(cdata_T))
  cdata_T <- rm_duplicates(cdata_T)

  if(file.exists(paste0(dataPath,"RNA_BCM_v1/",cancer,"_RNAseq_gene_RSEM_coding_UQ_1500_log2_Normal.txt")))
  {
    cdata_N <- as.matrix(read.csv(paste0(dataPath,"RNA_BCM_v1/",cancer,"_RNAseq_gene_RSEM_coding_UQ_1500_log2_Normal.txt"),sep="\t",row.names=1,check.names=F))

    rownames(cdata_N) <- geneAnno[rownames(cdata_N),2]
    rownames(cdata_N) <- transferSymbol(rownames(cdata_N))
    cdata_N <- rm_duplicates(cdata_N)

    olp <- intersect(rownames(cdata_T),rownames(cdata_N))
    cdata_T_minusBG <- cdata_T[olp,]-rowMeans(cdata_N[olp,])

  }else{
    cdata_T_minusBG <- cdata_T-rowMeans(cdata_T)
  }


  library(SecAct)
  sigs <- c(0,112,36,618)
  sigNames <- c("CytoSig","NicheNet","ImmuneDic","SecAct")

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
    save(res, file = paste0(validationPath,"CPTAC/",cancer,"_",sigNames[which(sigs==sig)],".RData"))
  }

}


dataPath <- "/data/rub2/data/CPTAC/"
cancers <- c("BRCA", "COAD", "GBM", "HNSC", "KIRC", "LUAD", "LUSC", "OV", "PDAC", "UCEC")

geneAnno <- read.csv(paste0(dataPath,"Proteome_BCM_GENCODE_v34_harmonized_v1/README/Gene_annotation_and_representable_isoform_mapping_table.txt"),sep="\t",check.names=F)
geneAnno[geneAnno[,1]=="ENSG00000284024.2",4] <- "MSANTD7"

geneAnno <- geneAnno[,c(1,4)]
geneAnno <- geneAnno[!duplicated(geneAnno[,1]),]
rownames(geneAnno) <- geneAnno[,1]

fg.df <- data.frame()

sigNames <- c("CytoSig","NicheNet","ImmuneDic","SecAct")
for(sigName in sigNames)
{
  for(cancer in cancers)
  {
    act_abu <- data.frame()
    act_mut <- data.frame()
    abu_mut <- data.frame()

    load(paste0(validationPath,"CPTAC/",cancer,"_",sigName,".RData"))
    Act <- as.matrix(res$zscore)

    if(sigName=="CytoSig")
    {
      rownames(Act) <- transferSymbol(rownames(Act))
      rownames(Act)[1] <- "INHBA"
    }

    Protein <- as.matrix(read.csv(paste0(dataPath,"Proteome_BCM_GENCODE_v34_harmonized_v1/",cancer,"_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt"),sep="\t",row.names=1,check.names=F))
    rownames(Protein) <- geneAnno[rownames(Protein),2]
    rownames(Protein) <- transferSymbol(rownames(Protein))
    Protein <- rm_duplicates(Protein)
    Protein <- Protein[rowSums(!is.na(Protein))>ncol(Protein)*0.2,,drop=F]

    olp_c <- intersect(colnames(Act),colnames(Protein))
    olp_r <- intersect(rownames(Act),rownames(Protein))

    Act_olp <- Act[olp_r,olp_c]
    Protein_olp <- Protein[olp_r,olp_c]

    genes <- rownames(Act_olp)
    corVec <- c()
    aucVec <- c()
    for(gene in genes)
    {
      x <- unlist(Act_olp[gene,])
      y <- unlist(Protein_olp[gene,])

      x <- x[!is.na(y)]
      y <- y[!is.na(y)]

      rv <- cor(x,y,use="pairwise.complete.obs")
      corVec <- c(corVec,rv)

      y_alt <- as.numeric(y >= quantile(y)[4]) ################ =?

      library(ROCR)
      predM <- prediction(x, y_alt)
      roc <-  performance(predM, measure = "auc")
      av <- unlist(roc@ y.values)

      aucVec <- c(aucVec,av)
    }

    fg.df[paste0(sigName,cancer,genes),"sig"] <- sigName
    fg.df[paste0(sigName,cancer,genes),"cancer"] <- cancer
    fg.df[paste0(sigName,cancer,genes),"gene"] <- genes
    fg.df[paste0(sigName,cancer,genes),"cor"] <- corVec
    fg.df[paste0(sigName,cancer,genes),"auc"] <- aucVec

  }
}




fg.df <- fg.df[!is.na(fg.df[,"auc"]),]

library(ggplot2)
library(patchwork)
library(reshape2)

fg.df_secact <- fg.df[fg.df[,"sig"]=="SecAct",]

p1 <- ggplot(fg.df_secact, aes(x = cancer, y = auc))+
  geom_hline(yintercept=0.5, color = "grey1", linewidth=0.8, linetype = 'dotted')+
  geom_violin(aes(group=cancer,fill=cancer), trim=FALSE, alpha=0.3)+
  geom_boxplot(width=0.15,outlier.shape=NA)+
  annotate("text", x = 6.1, y=0.08, label = paste0("n = 1043"))+
  ggtitle(" ")+
  xlab(" ")+
  ylab("AUC")+
  theme_classic()+
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position="none"
  )

ggsave(paste0(validationPath,"validation_CPTAC.png"), p1, width = 11, height = 7, dpi=200, units = "cm", limitsize = FALSE)







fg.df_gene <- data.frame()
for(gene in unique(fg.df[,"gene"]) )
{
  for(sig in unique(fg.df[,"sig"]) )
  {
    fg.df_sub <- fg.df[fg.df[,"gene"]==gene&fg.df[,"sig"]==sig,]

    if(nrow(fg.df_sub)==0) next

    fg.df_gene[paste0(gene,sig),"gene"] <- gene
    fg.df_gene[paste0(gene,sig),"sig"] <- sig
    fg.df_gene[paste0(gene,sig),"auc"] <- median(fg.df_sub[,"auc"],na.rm=T)
  }
}

sigOrder <- data.frame()
for(i in unique(fg.df_gene[,"sig"]))
{
  sigOrder[i,"r_mean"] <- mean(fg.df_gene[fg.df_gene[,"sig"]==i,"auc"])
}
sigOrder <- sigOrder[order(sigOrder[,"r_mean"],decreasing=T),,drop=F]

fg.df_gene[,2] <- factor(fg.df_gene[,2], levels=rownames(sigOrder))

sig_stat <- table(fg.df_gene[,"sig"])

p2 <- ggplot(fg.df_gene, aes(x = sig, y = auc))+
  geom_hline(yintercept=0.5, color = "grey1", linewidth=0.8, linetype = 'dotted')+
  geom_violin(aes(group=sig, fill=sig),width=0.8,trim=FALSE)+
  geom_boxplot(width=0.15,outlier.shape=NA)+
  scale_fill_manual( values=sigColors[rownames(sigOrder)] )+
  annotate("text", x = 1, y=0.1, label = paste0("n = "))+
  annotate("text", x = 1, y=0.02, label = paste0(sig_stat[rownames(sigOrder)[1]] ))+
  annotate("text", x = 2, y=0.02, label = paste0(sig_stat[rownames(sigOrder)[2]] ))+
  annotate("text", x = 3, y=0.02, label = paste0(sig_stat[rownames(sigOrder)[3]] ))+
  annotate("text", x = 4, y=0.02, label = paste0(sig_stat[rownames(sigOrder)[4]] ))+
  ylim(0,1)+
  ggtitle("All")+
  xlab(" ")+
  ylab("AUC")+
  theme_classic()+
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position="none"
  )

ggsave(paste0(validationPath,"validation_CPTAC_compare_all.png"), p2, width = 6.5, height = 7.8, dpi=200, units = "cm", limitsize = FALSE)



gene_stat <- table(fg.df_gene[,1])
olpGenes <- names(gene_stat)[gene_stat==length(sigNames)]


fg.df_gene <- data.frame()
for(gene in olpGenes)
{
  for(sig in unique(fg.df[,"sig"]) )
  {
    fg.df_sub <- fg.df[fg.df[,"gene"]==gene&fg.df[,"sig"]==sig,]

    if(nrow(fg.df_sub)==0) next

    fg.df_gene[paste0(gene,sig),"gene"] <- gene
    fg.df_gene[paste0(gene,sig),"sig"] <- sig
    fg.df_gene[paste0(gene,sig),"auc"] <- median(fg.df_sub[,"auc"],na.rm=T)
  }
}

sigOrder <- data.frame()
for(i in unique(fg.df_gene[,"sig"]))
{
  sigOrder[i,"r_mean"] <- mean(fg.df_gene[fg.df_gene[,"sig"]==i,"auc"])
}
sigOrder <- sigOrder[order(sigOrder[,"r_mean"],decreasing=T),,drop=F]

fg.df_gene[,2] <- factor(fg.df_gene[,2], levels=rownames(sigOrder))

p2 <- ggplot(fg.df_gene, aes(x = sig, y = auc))+
  geom_hline(yintercept=0.5, color = "grey1", linewidth=0.8, linetype = 'dotted')+
  geom_violin(aes(group=sig, fill=sig),width=0.8,trim=FALSE)+
  geom_boxplot(width=0.15,outlier.shape=NA)+
  scale_fill_manual( values=sigColors[rownames(sigOrder)] )+
  annotate("text", x = 2.5, y=0.05, label = paste0("n = ",length(olpGenes)))+
  ylim(0,1)+
  ggtitle("Shared")+
  xlab(" ")+
  ylab("AUC")+
  theme_classic()+
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position="none"
  )

ggsave(paste0(validationPath,"validation_CPTAC_compare_olp.png"), p2, width = 6.5, height = 7.8, dpi=200, units = "cm", limitsize = FALSE)
