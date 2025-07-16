source("Secretome_s0_path.R")
source("../../SpaCE/code/ifun.R")

#SWARM -t 2 -g 20 --time 03:00:00
args = commandArgs(trailingOnly=TRUE)
cancer <- args[1]

# CPTAC
dataPath <- "/data/rub2/data/CPTAC/"
cancers <- c("BRCA", "COAD", "GBM", "HNSC", "KIRC", "LUAD", "LUSC", "OV", "PDAC", "UCEC")

geneAnno <- read.csv(paste0(dataPath,"Proteome_BCM_GENCODE_v34_harmonized_v1/README/Gene_annotation_and_representable_isoform_mapping_table.txt"),sep="\t",check.names=F)
geneAnno[geneAnno[,1]=="ENSG00000284024.2",4] <- "MSANTD7"

geneAnno <- geneAnno[,c(1,4)]
geneAnno <- geneAnno[!duplicated(geneAnno[,1]),]
rownames(geneAnno) <- geneAnno[,1]


# CPTAC RNA
#for(cancer in cancers)
#{
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
	
	
	#write.table(cdata_T_minusBG,paste0(validationPath,"CPTAC/",cancer,".diff"),sep="\t",quote=F)	
	#
	#runCytoSig <- 
	#	paste0(
	#	"CytoSig_run.py -a 1000000 -s 15 -i ",
	#	paste0(validationPath,"CPTAC/",cancer,".diff"),
	#	" -o ",
	#	paste0(validationPath,"CPTAC/",cancer)
	#	)
	#
	#system("module load python") # run it in node first
	#system(runCytoSig)
	
	
	library(SecAct)
	sigs <- c(0,11,112,36,618)
	sigNames <- c("CytoSig","NicheNet.v1","NicheNet.v2","ImmuneDic","SecAct")
	
	sigs <- c(61811,61822,61833,61844,61855,61866,61877)
	sigNames <- c("SecAct.w1","SecAct.w2","SecAct.v1","SecAct.v2","SecAct.v2_0.9","SecAct.v2_0.8")
	
	sigs <- c(61888)
	sigNames <- c("SecAct.1k_0.9")
	
	sigs <- c(61899)
	sigNames <- c("SecAct.1k_residual")

	sigs <- c(100011,100012,100021,100022,100031,100032,100041,100042,100051,100052)
	sigNames <- c("SecAct.1k_11","SecAct.1k_12","SecAct.1k_21","SecAct.1k_22","SecAct.1k_31","SecAct.1k_32","SecAct.1k_41","SecAct.1k_42","SecAct.1k_51","SecAct.1k_52")
	
	
	for(sig in sigs)
	{
		if(sig==0)
		{
			ref <- paste0("/data/rub2/project/Secretome/code/CytoSig-master/CytoSig/signature.centroid")
		}else if(sig==618){
			ref <- paste0("/data/rub2/project/Secretome/results/signatureComb/AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3.tsv")
		}else if(sig==6182){
			ref <- paste0("/data/rub2/project/Secretome/results/signatureComb/AllSigFilteredBy_MoranI_TCGA_0.25_ds3_ICGC_highConfidence.tsv")
		}else if(sig==61811){
			ref <- paste0("/data/rub2/project/Secretome/results/signatureComb/AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3_weighted1.tsv")
		}else if(sig==61822){
			ref <- paste0("/data/rub2/project/Secretome/results/signatureComb/AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3_weighted2.tsv")
		}else if(sig==61833){
			ref <- paste0("/data/rub2/project/Secretome/results/signatureComb/SecAct.signature.adjusted_v1")
		}else if(sig==61844){
			ref <- paste0("/data/rub2/project/Secretome/results/signatureComb/SecAct.signature.adjusted_v2")
		}else if(sig==61855){
			ref <- paste0("/data/rub2/project/Secretome/results/signatureComb/SecAct.signature.adjusted_v2_0.1")
		}else if(sig==61866){
			ref <- paste0("/data/rub2/project/Secretome/results/signatureComb/SecAct.signature.adjusted_v2_0.2")
		}else if(sig==61877){
			ref <- paste0("/data/rub2/project/Secretome/results/signatureComb/AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3_weighted2_0.1")
		}else if(sig==61888){
			ref <- paste0("/data/rub2/project/Secretome/results/signatureComb/AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3_1k_0.1")
		}else if(sig==61899){
			ref <- paste0("/data/rub2/project/Secretome/results/signatureComb/AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3_1k_residual_0.1.tsv")
		}else if(sig==100011){
			ref <- paste0(signatureCombPath,"AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3_1k_vst.tsv")
		}else if(sig==100012){
			ref <- paste0(signatureCombPath,"AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3_1k_vst_0.9.tsv")
		}else if(sig==100021){
			ref <- paste0(signatureCombPath,"AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3_1k_weighted.tsv")
		}else if(sig==100022){
			ref <- paste0(signatureCombPath,"AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3_1k_weighted_0.9.tsv")
		}else if(sig==100031){
			ref <- paste0(signatureCombPath,"AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3_1k_weighted2.tsv")
		}else if(sig==100032){
			ref <- paste0(signatureCombPath,"AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3_1k_weighted2_0.9.tsv")
		}else if(sig==100041){
			ref <- paste0(signatureCombPath,"AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3_1k_vst_condition_logUMI_cellType.tsv")
		}else if(sig==100042){
			ref <- paste0(signatureCombPath,"AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3_1k_vst_condition_logUMI_cellType_0.9.tsv")
		}else if(sig==100051){
			ref <- paste0(signatureCombPath,"AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3_1k_vst_condition_cellType.tsv")
		}else if(sig==100052){
			ref <- paste0(signatureCombPath,"AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3_1k_vst_condition_cellType_0.9.tsv")
		}else{
			ref <- paste0("/data/rub2/project/Secretome/code/CytoSig-master/CytoSig/signature.centroid.",sig)
		}
		
		#aMat <- read.csv(paste0(QCSummaryPath,"lambda_summary/lambda_compare_",sig,".csv"))
		#a <- aMat[aMat[,2]==max(aMat[,2]),1]
		a <- 5e+05
		
		res <- SecAct.inference(Y=cdata_T_minusBG, SigMat=ref, lambda=a, nrand=1000)
		save(res, file = paste0(validationPath,"CPTAC/",cancer,"_",sigNames[which(sigs==sig)],".RData"))
	}

#}


if(FALSE)
{

dataPath <- "/data/rub2/data/CPTAC/"
cancers <- c("BRCA", "COAD", "GBM", "HNSC", "KIRC", "LUAD", "LUSC", "OV", "PDAC", "UCEC")

geneAnno <- read.csv(paste0(dataPath,"Proteome_BCM_GENCODE_v34_harmonized_v1/README/Gene_annotation_and_representable_isoform_mapping_table.txt"),sep="\t",check.names=F)
geneAnno[geneAnno[,1]=="ENSG00000284024.2",4] <- "MSANTD7"

geneAnno <- geneAnno[,c(1,4)]
geneAnno <- geneAnno[!duplicated(geneAnno[,1]),]
rownames(geneAnno) <- geneAnno[,1]

fg.df <- data.frame()

sigNames <- c("CytoSig","NicheNet.v1","NicheNet.v2","ImmuneDic","SecAct")
sigNames <- c("CytoSig","NicheNet.v1","NicheNet.v2","ImmuneDic","SecAct","SecAct.w1","SecAct.w2","SecAct.w2_0.9","SecAct.v1","SecAct.v2","SecAct.v2_0.9","SecAct.v2_0.8","SecAct.1k_0.9","SecAct.1k_residual","SecAct.1k_11","SecAct.1k_12","SecAct.1k_21","SecAct.1k_22","SecAct.1k_31","SecAct.1k_32","SecAct.1k_41","SecAct.1k_42","SecAct.1k_51","SecAct.1k_52")
sigNames <- c("CytoSig","NicheNet.v1","NicheNet.v2","ImmuneDic","SecAct.1k_42")
for(sigName in sigNames)
{
for(cancer in cancers)
{
	act_abu <- data.frame()
	act_mut <- data.frame()
	abu_mut <- data.frame()
	
	#Act <- as.matrix(read.table(paste0(validationPath,"CPTAC/",cancer,".Zscore"),sep="\t",check.names=F)	)
	load(paste0(validationPath,"CPTAC/",cancer,"_",sigName,".RData"))
	Act <- as.matrix(res$zscore)
	Act <- expand_rows(Act)
	
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
		
		y_alt <- as.numeric(y >= quantile(y)[4] ) ################ =?
		
		library(ROCR)
		predM <- prediction(x, y_alt)
		roc <-  performance(predM, measure = "auc")
		av <- unlist(roc@ y.values)
		
		aucVec <- c(aucVec,av)
		
		
		if(sigName=="SecAct.1k_42"&cancer=="BRCA"&gene=="GZMA")
		{
			roc.df <- data.frame(expr=as.character(y_alt),act=x)
			
			library(ggplot2)
			p2 <- ggplot(roc.df,aes(x=expr,y=act)) + 
				geom_jitter(aes(fill=expr),colour="black",width=0.2,shape = 21,alpha=0.6)+
				scale_color_manual(values=c("#00BFC4","#F8766D"))+
				scale_x_discrete(labels= c("Low","High"))+
				ylab("GZMA SecAct Activity")+
				xlab("GZMA Proteomics Abundance")+
				theme_classic()+ 
				theme(
					panel.grid = element_blank(),
			  		panel.background = element_blank(),
			  		plot.title = element_text(hjust = 0.5),
    				axis.title = element_text(colour = "black"),
					axis.text = element_text(colour = "black"),
					legend.position="none"
				)
			ggsave(paste0(validationPath,"validation_CPTAC_single_jitter.png"), p2, width = 6.8, height = 6.5, dpi=200, units = "cm",limitsize = FALSE)
			
			
			library(pROC)
			library(ggplot2)
			rocobj <- roc(expr ~ act, data = roc.df)
			auc <- round(auc(expr ~ act, data = roc.df),2)

			# plot on a single plot with AUC in labels
			p2 <- ggroc(rocobj,legacy.axes=TRUE, color="purple")+ 
				ylab("True positive rate")+
				xlab("False positive rate")+
				annotate("text", x = 0.65, y=0.3, label = paste0("AUC = ",auc))+
				geom_abline(slope=1, intercept=0, linetype="dashed")+
				theme_classic()+ 
				theme(
					panel.grid = element_blank(),
			  		panel.background = element_blank(),
			  		plot.title = element_text(hjust = 0.5),
    				axis.title = element_text(colour = "black"),
					axis.text = element_text(colour = "black"),
					legend.position="none"
				)
			ggsave(paste0(validationPath,"validation_CPTAC_single_ROC.png"), p2, width = 6.8, height = 6.5, dpi=200, units = "cm",limitsize = FALSE)

			
		}
		
	}
	
	fg.df[paste0(sigName,cancer,genes),"sig"] <- sigName
	fg.df[paste0(sigName,cancer,genes),"cancer"] <- cancer
	fg.df[paste0(sigName,cancer,genes),"gene"] <- genes
	fg.df[paste0(sigName,cancer,genes),"cor"] <- corVec
	fg.df[paste0(sigName,cancer,genes),"auc"] <- aucVec

	
}
}


fg.df <- fg.df[!is.na(fg.df[,"auc"]),]
fg.df[fg.df[,"sig"]=="SecAct.1k_42","sig"] <- "SecAct"

gn <- length(unique(fg.df[,"gene"]))


library(ggplot2)
library(patchwork)
library(reshape2)

fg.df_secact <- fg.df[fg.df[,"sig"]=="SecAct",]


p1 <- ggplot(fg.df_secact, aes(x = cancer, y = auc))+
	geom_hline(yintercept=0.5, color = "grey1", linewidth=0.75, linetype = 'dotted')+
	geom_violin(aes(group=cancer,fill=cancer), trim=FALSE, alpha=0.3, width=0.8)+
	geom_boxplot(width=0.15,outlier.shape=NA)+
	annotate("text", x = 6.1, y=0.08, label = paste0("n = 952"))+
	ggtitle("SecAct")+
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

ggsave(paste0(validationPath,"validation_CPTAC.png"), p1, width = 12.1, height = 7, dpi=200, units = "cm", limitsize = FALSE)


cancer <- "BRCA"
gene <- "GZMA"
fg.df_secact[fg.df_secact[,"cancer"]==cancer,] -> fg.df_secact_brca
fg.df_secact_brca[order(fg.df_secact_brca[,5]),]



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
	sigOrder[i,"r_mean"] <- median(fg.df_gene[fg.df_gene[,"sig"]==i,"auc"])
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
	annotate("text", x = 5, y=0.02, label = paste0(sig_stat[rownames(sigOrder)[5]] ))+
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

ggsave(paste0(validationPath,"validation_CPTAC_compare_all_weighted.png"), p2, width = 17, height = 7.8, dpi=200, units = "cm", limitsize = FALSE)
ggsave(paste0(validationPath,"validation_CPTAC_compare_all.png"), p2, width = 7, height = 7.8, dpi=200, units = "cm", limitsize = FALSE)



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
	sigOrder[i,"r_mean"] <- median(fg.df_gene[fg.df_gene[,"sig"]==i,"auc"])
}
sigOrder <- sigOrder[order(sigOrder[,"r_mean"],decreasing=T),,drop=F]

fg.df_gene[,2] <- factor(fg.df_gene[,2], levels=rownames(sigOrder))

p2 <- ggplot(fg.df_gene, aes(x = sig, y = auc))+
	geom_hline(yintercept=0.5, color = "grey1", linewidth=0.8, linetype = 'dotted')+
	geom_violin(aes(group=sig, fill=sig),width=0.8,trim=FALSE)+
	geom_boxplot(width=0.15,outlier.shape=NA)+
	scale_fill_manual( values=sigColors[rownames(sigOrder)] )+
	annotate("text", x = 3, y=0.05, label = paste0("n = ",length(olpGenes)))+
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

ggsave(paste0(validationPath,"validation_CPTAC_compare_olp.png"), p2, width = 7, height = 7.8, dpi=200, units = "cm", limitsize = FALSE)






#Mutation <- as.matrix(read.csv(paste0(dataPath,"Mutation_BCM_v1/",cancer,"_somatic_mutation_gene_level_binary.txt"),sep="\t",row.names=1,check.names=F))
#rownames(Mutation) <- geneAnno[rownames(Mutation),2]
#rownames(Mutation) <- transferSymbol(rownames(Mutation))
#Mutation <- rm_duplicates(Mutation)
#
#
#Mutation <- Mutation[rownames(Mutation)%in%genes,]
#Mutation <- Mutation[,colnames(Mutation)%in%olp_c]
#
#olp_c2 <- intersect(colnames(Mutation),olp_c)
#
#
#highly_mutated_genes <- rownames(Mutation)[rowSums(Mutation)>= length(olp_c2)*0.05]
#if(length(highly_mutated_genes)==0) next
#
#
#for(gene in highly_mutated_genes)
#{
#	act_vec <- RNA_olp[gene, olp_c2 ]
#	abu_vec <- Protein_olp[gene, olp_c2 ]
#	mut <- as.character(Mutation[gene, olp_c2 ])
#			
#	pv1 <- wilcox.test(act_vec~mut)$p.value
#	
#	#pv2 <- wilcox.test(abu_vec~mut)$p.value
#	
#	if(pv1<0.05)
#	{
#		cor_res <- cor.test(Protein_olp[gene,],RNA_olp[gene,])
#		rv <- round(cor_res$estimate,3)
#		pv <- signif(cor_res$p.value,3)
#		act_abu[paste0(gene,olp_c),"gene"] <- paste0(gene,", r=",rv,", p=",pv)
#		act_abu[paste0(gene,olp_c),"abu"] <- Protein_olp[gene,]
#		act_abu[paste0(gene,olp_c),"act"] <- RNA_olp[gene,]
#		
#		act_mut[paste0(gene,olp_c2),"gene"] <- gene
#		act_mut[paste0(gene,olp_c2),"mut"] <- as.character(Mutation[gene, olp_c2 ])
#		act_mut[paste0(gene,olp_c2),"act"] <- RNA_olp[gene, olp_c2 ]
#		
#		abu_mut[paste0(gene,olp_c2),"gene"] <- gene
#		abu_mut[paste0(gene,olp_c2),"mut"] <- as.character(Mutation[gene, olp_c2 ])
#		abu_mut[paste0(gene,olp_c2),"abu"] <- Protein_olp[gene, olp_c2 ]
#	}
#}
#
#if(nrow(act_mut)==0) next
#highly_mutated_genes <- highly_mutated_genes[highly_mutated_genes%in%unique(act_mut[,"gene"])]
#if(length(highly_mutated_genes)==0) next
#
#
#library(ggplot2)
#library(patchwork)
#library(ggsignif)
#p1 <- ggplot(act_abu,aes(x = abu, y = act))+
#	geom_point()+
#	theme_bw()+ 
#	theme(
#	  panel.grid = element_blank(),
#	  panel.background = element_blank(),
#	  axis.text = element_text(size=10,colour = "black"),
#	  axis.title = element_text(size=10,colour = "black"),
#	  legend.position="none"
#	)+ facet_wrap(~ gene, ncol =length(highly_mutated_genes), scales = "free") 
#
#p2 <- ggplot(act_mut,aes(x = mut, y = act))+
#	geom_boxplot(outlier.shape=NA)+
#	geom_jitter()+
#	geom_signif( comparisons = list(c("0", "1"))  )+
#	theme_bw()+ 
#	theme(
#	  panel.grid = element_blank(),
#	  panel.background = element_blank(),
#	  axis.text = element_text(size=10,colour = "black"),
#	  axis.title = element_text(size=10,colour = "black"),
#	  legend.position="none"
#	)+ facet_wrap(~ gene, ncol =length(highly_mutated_genes), scales = "free") 
#
#p3 <- ggplot(abu_mut,aes(x = mut, y = abu))+
#	geom_boxplot(outlier.shape=NA)+
#	geom_jitter()+
#	geom_signif( comparisons = list(c("0", "1"))  )+
#	theme_bw()+ 
#	theme(
#	  panel.grid = element_blank(),
#	  panel.background = element_blank(),
#	  axis.text = element_text(size=10,colour = "black"),
#	  axis.title = element_text(size=10,colour = "black"),
#	  legend.position="none"
#	)+ facet_wrap(~ gene, ncol =length(highly_mutated_genes), scales = "free") 
#
#
#ggsave(paste0(validationPath,"/mutation/",cancer,".png"), p1/p2/p3, width = 11*length(highly_mutated_genes), height = 33, dpi=200, units = "cm",limitsize = FALSE)


}	
