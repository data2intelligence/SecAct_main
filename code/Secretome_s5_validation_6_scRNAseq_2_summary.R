source("Secretome_s0_path.R")


inputPath <- paste0(dataValPath,"scRNAseq/")
outputPath <- paste0(validationPath,"scRNAseq/")

cancers <- list.dirs(outputPath,full.names = FALSE,recursive = FALSE)
cancers <- setdiff(cancers,"") #includes the input path itself
cancers <- setdiff(cancers,c("Nasopharyngeal.GSE162025.10x","Breast.GSE176078.10x")) #includes the input path itself
cancers <- setdiff(cancers,c("Fallopian.GSE156728.10x","Cholangio.GSE156728.SS2")) #includes the input path itself

AUPRC_baseline_vec <- c()
summaryTable.median.median <- data.frame() # sig ~ cancer
summaryTable.median.median2 <- data.frame() # sig ~ SPTF


for(sigName in sigNames)
{
	summaryTable <- data.frame()
	for(cancer in cancers)
	{
		temp <- read.csv(paste0(outputPath,cancer,"/",sigName,"_SP_TF.csv"),row.names=1,header=TRUE)
		summaryTable <- rbind(summaryTable,temp)
	}
	AUPRC_baseline_vec[sigName] <- median(summaryTable[,"AUPRC_baseline"])
	
	summaryTable.median <- data.frame()
	for(cancer in cancers)
	{
		summaryTableCancer <- summaryTable[summaryTable[,"cancer"]==cancer,]
		for(i in unique(summaryTableCancer[,"SP_TF"]))
		{
			summaryTable_sub <- summaryTableCancer[summaryTableCancer[,"SP_TF"]==i,]
			
			summaryTable.median[paste0(cancer,i),"cancer"] <- cancer
			summaryTable.median[paste0(cancer,i),"SP_TF"] <- i
			summaryTable.median[paste0(cancer,i),"AUC.median"] <- median(summaryTable_sub[,"AUC"])
			summaryTable.median[paste0(cancer,i),"AUPRC.median"] <- median(summaryTable_sub[,"AUPRC"])
		}
	}
	
	for(cancer in cancers)
	{
		summaryTable.median_sub <- summaryTable.median[summaryTable.median[,"cancer"]==cancer,]
		
		summaryTable.median.median[paste0(sigName,cancer),"sig"] <- sigName
		summaryTable.median.median[paste0(sigName,cancer),"cancer"] <- cancer
		summaryTable.median.median[paste0(sigName,cancer),"AUC.median.median"] <- median(summaryTable.median_sub[,"AUC.median"])
		summaryTable.median.median[paste0(sigName,cancer),"AUPRC.median.median"] <- median(summaryTable.median_sub[,"AUPRC.median"])
	}
	
	for(SPTF in unique(summaryTable.median[,"SP_TF"]))
	{
		summaryTable.median_sub <- summaryTable.median[summaryTable.median[,"SP_TF"]==SPTF,]
		
		summaryTable.median.median2[paste0(sigName,SPTF),"sig"] <- sigName
		summaryTable.median.median2[paste0(sigName,SPTF),"SPTF"] <- SPTF
		summaryTable.median.median2[paste0(sigName,SPTF),"AUC.median.median"] <- median(summaryTable.median_sub[,"AUC.median"])
		summaryTable.median.median2[paste0(sigName,SPTF),"AUPRC.median.median"] <- median(summaryTable.median_sub[,"AUPRC.median"])
	}
	
	
	# draw heatmap summary
	mat = reshape2::dcast( summaryTable.median[,1:3], cancer~SP_TF )
	rownames(mat) <- mat[,1]
	mat <- t(mat[,-1])
	
	print(sigName)
	res <- wilcox.test(apply(mat,1,function(x) median(x,na.rm=T)), mu=0.5)
	print(res$p.value)
	
	mat <- mat[order(apply(mat,1,function(x) median(x,na.rm=T)),decreasing=T),]
	mat <- mat[,sort(colnames(mat))]
	
	mat <- t(mat)
	colnames(mat) <- gsub("and","&", colnames(mat))
	rownames(mat) <- gsub("_","", rownames(mat))
	
	widthValue <- ifelse(sigName=="ImmuneDic",19.9,21)
	
	write.csv(mat, paste0(outputPath,"validation_scRNAseq_",sigName,"_heatmap.csv"), quote=FALSE)	
	pdf(paste0(outputPath,"validation_scRNAseq_",sigName,"_heatmap.pdf"), width = 6.9, height = 7.8)
	
	library(ComplexHeatmap)
	
	technique_vec <- sapply(strsplit(rownames(mat),".",fixed=T),function(x) return(x[3])) 	
	dataset_vec <- sapply(strsplit(rownames(mat),".",fixed=T),function(x) return(paste0(x[1],"_",x[2],"          "))) 	
	
	row_ha <- rowAnnotation(
		Sequencing = technique_vec,
		col = list(Sequencing = c("InDrop" = "#B395BD", "10x" = "#7DAEE0", "SS2" = "#EA8379")),
		labels = anno_text(dataset_vec, which = "row", gp = gpar(fontsize = 10) ), 
		width = max(grobWidth(textGrob(dataset_vec)))
	)
	column_ha <- columnAnnotation(
		AUC = anno_boxplot(as.matrix(mat), height = unit(2.8, "cm"), gp = gpar(fill = "skyblue")),
		annotation_name_side = "left"
	)
			
	ht <- Heatmap(as.matrix(mat),
		name = "AUC",
		rect_gp = gpar(col = "grey10", lwd = 2),
		row_title = " \n \n \n \n ",
		column_title = sigName,
		col = circlize::colorRamp2(c(0,0.5,1), c("#aad962", "white", "#ef6a32")),
		row_names_max_width = max_text_width(
	        rownames(mat), 
	        gp = gpar(fontsize = 10)
	        ),
	    column_names_max_height = max_text_width(
	        colnames(mat), 
	        gp = gpar(fontsize = 11)
	        ),
	    show_row_names = FALSE,
	    show_column_names = TRUE,
	    top_annotation = column_ha,
	    right_annotation = row_ha,
	    column_names_rot = 45,
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


for(aa in c("AUC","AUPRC"))
{	
	sigOrder <- data.frame()
	for(i in unique(summaryTable.median.median2[,"sig"]))
	{
		sigOrder[i,"r_mean"] <- median(summaryTable.median.median2[summaryTable.median.median2[,"sig"]==i,paste0(aa,".median.median")],na.rm=TRUE)
	}
	sigOrder <- sigOrder[order(sigOrder[,"r_mean"],decreasing=TRUE),,drop=FALSE]
	
	summaryTable.median.median2_sub <- summaryTable.median.median2[,c("sig","SPTF",paste0(aa,".median.median"))]
	colnames(summaryTable.median.median2_sub)[3] <- "value"
	summaryTable.median.median2_sub[,1] <- factor(summaryTable.median.median2_sub[,1], levels=rownames(sigOrder))
	
	library(ggplot2)
	p2 <- ggplot(summaryTable.median.median2_sub, aes(x = sig, y = value, fill=sig))+
		geom_hline(yintercept=ifelse(aa=="AUC",0.5,AUPRC_baseline_vec["SecAct"]), color = "grey1", linewidth=0.5, linetype = 'dotted')+
		geom_boxplot(width=0.7,outlier.shape=NA)+
		scale_fill_manual( values=sigColors[rownames(sigOrder)] )+
		geom_jitter(size=0.8,alpha=0.5,color="red",width=0.01)+
		ggtitle(" ")+
		xlab(" ")+
		ylab(aa)+
		theme_classic()+ 
		theme(
		  panel.grid = element_blank(),
		  panel.background = element_blank(),
		  plot.title = element_text(hjust = 0.5),
		  axis.text = element_text(size=10,colour = "black"),
		  axis.title = element_text(size=10,colour = "black"),
		  axis.text.x = element_text(size=10,angle = 45, hjust = 1),
		  axis.title.x = element_blank(),
		  legend.position="none"
		)
	if(aa=="AUC") p2 <- p2+scale_y_continuous(breaks = c(0, 0.5, 1))
	
	ggsave(paste0(outputPath,"validation_scRNAseq_comb_SP_TF_",aa,".png"), p2, width = 6.2, height = 7.7, dpi=400, units = "cm",limitsize = FALSE)
	write.csv(summaryTable.median.median2_sub, paste0(outputPath,"validation_scRNAseq_comb_SP_TF_",aa,".csv"), quote=FALSE)

	
	#sigOrder <- data.frame()
	#for(i in unique(summaryTable.median.median[,"sig"]))
	#{
	#	sigOrder[i,"r_mean"] <- median(summaryTable.median.median[summaryTable.median.median[,"sig"]==i,paste0(aa,".median.median")],na.rm=TRUE)
	#}
	#sigOrder <- sigOrder[order(sigOrder[,"r_mean"],decreasing=TRUE),,drop=FALSE]
	#
	#
	#summaryTable.median.median_sub <- summaryTable.median.median[,c("sig","cancer",paste0(aa,".median.median"))]
	#colnames(summaryTable.median.median_sub)[3] <- "value"
	#summaryTable.median.median_sub[,1] <- factor(summaryTable.median.median_sub[,1], levels=rownames(sigOrder))
	#
	#library(ggplot2)
	#p2 <- ggplot(summaryTable.median.median_sub, aes(x = sig, y = value, fill=sig))+
	#	geom_hline(yintercept=0.5, color = "grey1", linewidth=0.8, linetype = 'dotted')+
	#	geom_boxplot(width=0.7,outlier.shape=NA)+
	#	scale_fill_manual( values=sigColors[rownames(sigOrder)] )+
	#	geom_jitter(size=0.1,alpha=0.2)+
	#	ggtitle(" ")+
	#	xlab(" ")+
	#	ylab(aa)+
	#	theme_bw()+ 
	#	theme(
	#	  panel.grid = element_blank(),
	#	  panel.background = element_blank(),
	#	  plot.title = element_text(hjust = 0.5),
	#	  axis.text = element_text(size=10,colour = "black"),
	#	  axis.title = element_text(size=10,colour = "black"),
	#	  axis.text.x = element_text(size=10,angle = 45, hjust = 1),
	#	  legend.position="none"
	#	)
	#
	#ggsave(paste0(outputPath,"validation_scRNAseq_comb_cancer_",aa,".png"), p2, width = 6.2, height = 8, dpi=200, units = "cm",limitsize = FALSE)

}



##############
# fig 3h, i, j
##############

sigName <- "SecAct"
cancers <- c("Head_Neck.GSE103322.SS2")

sp_tf <- list(
	IFN.STAT1=list(c("IFN1","IFNG","IFNL","IFNL1","IL21","IL27") , "STAT1"),
	#IFN1.STAT2=list(c("IFNAR2") , "STAT2"),
	IL6and10fam.STAT3=list(c("IL6","LIF","OSM","IL10","IL22","IL9","IL21","IL27","IL17A","IL23A") , "STAT3"),
	IL12fam.STAT4=list(c("IL12","IL23A","IL27") , "STAT4"), 
	IL2fam.STAT5=list(c("IL2","IL3","IL7","IL9","IL15","IL27") , "STAT5"), 
	#IL4.STAT6=list(c("IL4R") , "STAT6"), #IL4, IL13
	IL1andTNF.RELA=list(c("IL1A","IL1B","TNF") , "RELA"),
	TNFSFfam.RELB=list(c("LTA","CD40LG","TNFSF11","TNFSF12","TNFSF13B","TNFSF14") , "RELB"),
	BMPfam.SMAD1=list(c("BMP2","BMP3","BMP4","BMP5","BMP6","BMP7","BMP8A") , "SMAD1"),
	TGFBfam.SMAD2=list(c("INHBA","GDF11","TGFB1","TGFB2","TGFB3") , "SMAD2"),
	TGFBfam.SMAD3=list(c("INHBA","GDF11","TGFB1","TGFB2","TGFB3") , "SMAD3"),
	BMPandTGFBfam.SMAD4=list(c("BMP2","BMP3","BMP4","BMP5","BMP6","BMP7","BMP8A","INHBA","GDF11","TGFB1","TGFB2","TGFB3") , "SMAD4")
	#TGFB.SMAD5=list(c("BMP2","BMP3","BMP4","BMP5","BMP6","BMP7","BMP8A") , "SMAD5") # few cells have SMAD5 activity
)

for(cancer in cancers)
{
	summaryTable <- data.frame()
	
	load(paste0(outputPath,cancer,"/",sigName,".RData"))
	sp <- res$zscore
	
	tf <- read.table(paste0(inputPath,cancer,".Rabit.Cistrome.t.feature.gz"),sep="\t",check.names=F)
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
		
		for(i in names(sp_tf)[7])
		{
			if(sum(rownames(sp_olp_sub)%in%sp_tf[[i]][[1]])==0) next
			
			sp_act <- colMeans(sp_olp_sub[rownames(sp_olp_sub)%in%sp_tf[[i]][[1]],,drop=F])
			
			if(!sp_tf[[i]][[2]]%in%rownames(tf_olp_sub)) next
			tf_act <- colMeans(tf_olp_sub[sp_tf[[i]][[2]],,drop=F])
			
			if(sum(tf_act>0)==0 | sum(tf_act<0)==0) next
			fg.df <- data.frame(sp_act=sp_act,tf_act=tf_act)
			fg.df <- fg.df[fg.df[,2]!=0,]
			
			if(sum(fg.df[,2]>0)<3) next
			if(sum(fg.df[,2]<0)<3) next
			if(nrow(fg.df)<10) next
			
			x1 <- fg.df[fg.df[,2]>0,1]
			x2 <- fg.df[fg.df[,2]<0,1]
			pv <- signif(wilcox.test(x1, x2)$p.value,2)
			
			
			library(ggplot2)
			p1 <- ggplot(fg.df,aes(x=tf_act, y=sp_act)) + 
				geom_point(color="#551177", alpha=0.4, size=0.6)+
				annotate("text", x = 2.2, y=-4, label = paste0("p = ",pv))+
				ggtitle(" ")+
				xlab("SMAD1 activity")+
				ylab("BMP family activity")+
				theme_classic()+ 
				theme(
					plot.background = element_blank(),
					panel.grid = element_blank(),
					plot.title = element_text(hjust = 0.5),
					axis.title = element_text(colour = "black"),
					axis.text = element_text(colour = "black"),
					legend.position="none"
				)			
			ggsave(paste0(outputPath,"validation_scRNAseq_",cancer,"_",i,"_",celltype,"_scatter.png"), p1, width = 5.9, height = 7, dpi=300, units = "cm")
			#ggsave(paste0(validationPath,"validation_singlecell_",cancer,"_",i,"_",celltype,"_scatter.pdf"), p1, width = 5.9, height = 6, dpi=400, units = "cm")
			if(celltype=="Tumor_Cancer") write.csv(fg.df, paste0(outputPath,"validation_scRNAseq_",cancer,"_",i,"_",celltype,"_scatter.csv"), quote=FALSE)
			
			
			fg.df[,2] <- fg.df[,2]>0
			if(length(unique(fg.df[,2]))==1) next
			if(nrow(fg.df)<10) next
			
			library(ROCR)
			predM <- prediction(fg.df$sp_act, fg.df$tf_act)
			roc = performance(predM, measure = "auc")
			roc_plot = performance(predM, measure = "tpr", x.measure = "fpr")
			
			summaryTable[paste0(cancer,celltype,i),"cancer"] <- cancer
			summaryTable[paste0(cancer,celltype,i),"celltype"] <- celltype
			summaryTable[paste0(cancer,celltype,i),"SP_TF"] <- i
			summaryTable[paste0(cancer,celltype,i),"AUC"] <- roc@ y.values
			
			
			library(pROC)
			library(ggplot2)
			rocobj <- roc(tf_act ~ sp_act, data = fg.df)
			auc <- round(auc(tf_act ~ sp_act, data = fg.df),2)
	
			p2 <- ggroc(rocobj,legacy.axes=TRUE, color="purple")+ 
				ylab("True positive rate")+
				xlab("False positive rate")+
				annotate("text", x = 0.65, y=0.3, label = paste0("AUC = ",auc))+
				geom_abline(slope=1, intercept=0, linetype="dashed")+
				scale_x_continuous(breaks = c(0, 0.5, 1))+
				scale_y_continuous(breaks = c(0, 0.5, 1))+
				theme_classic()+ 
				theme(
					panel.grid = element_blank(),
			  		panel.background = element_blank(),
			  		plot.title = element_text(hjust = 0.5),
	    			axis.title = element_text(colour = "black"),
					axis.text = element_text(colour = "black"),
					legend.position="none"
				)
			ggsave(paste0(outputPath,"validation_scRNAseq_",cancer,"_",i,"_",celltype,"_ROC.png"), p2, width = 6.8, height = 6.5, dpi=200, units = "cm",limitsize = FALSE)
	
		}
	}
	
	if(nrow(summaryTable)==0) next
	summaryTable <- summaryTable[order(summaryTable[,"AUC"],decreasing=T),]
	summaryTable[,"celltype"] <- factor(summaryTable[,"celltype"],levels=as.character(summaryTable[,"celltype"])	)
	
	library(ggplot2)
	p1 <- ggplot(summaryTable, aes(x=celltype, y=AUC)) +
		geom_hline(yintercept=0.5, color = "grey10", linewidth=0.8, linetype = 'dotted')+
	  	geom_bar(stat="identity", fill="#eca680", color="grey10", width=0.7, alpha=0.25) +
		#geom_hline(yintercept=AUCmedian, color = "red", linewidth=0.8, linetype = 'dotted')+
		#annotate("text", x = 3, y=0.85, color="red", label = paste0("Median = ",round(AUCmedian,2)))+
	  	ggtitle(cancer)+
	  	ylab("AUC")+
		coord_cartesian(ylim = c(0.4, 1))+
		theme_classic()+ 
		theme(
			panel.background = element_blank(),
			panel.grid = element_blank(),
			axis.title = element_text(colour = "black"),
			axis.text = element_text(colour = "black"),
			axis.text.x = element_text(angle = 38,hjust = 1),
			axis.title.x =element_blank(),
			legend.position = "none"
		)
	ggsave(paste0(outputPath,"validation_scRNAseq_",cancer,"_",i,"_",sigName,"_ROC_comb.png"), p1, width = 20, height = 7.5, dpi=500, units = "cm", limitsize =FALSE)


} # cancer




summaryTable <- summaryTable[summaryTable[,2]!="Unknown",]
summaryTable[,2] <- sapply(strsplit(as.character(summaryTable[,2]),"_",fixed=T),function(x) return(x[2]))

library(dplyr)
sumsumMean <- summaryTable %>% 
	group_by(celltype) %>% 
	summarise(AUC = mean(AUC,na.rm=TRUE))
	
summaryTable <- as.data.frame(sumsumMean)

summaryTable <- summaryTable[order(summaryTable[,"AUC"],decreasing=T),]
AUCmedian <- median(summaryTable[,"AUC"])
summaryTable[,"celltype"] <- factor(summaryTable[,"celltype"],levels=as.character(summaryTable[,"celltype"])	)

library(ggplot2)
p1 <- ggplot(summaryTable, aes(x=celltype, y=AUC)) +
	geom_hline(yintercept=0.5, color = "grey10", linewidth=0.8, linetype = 'dotted')+
  	geom_bar(stat="identity", fill="#eca680", color="grey10", width=0.7, alpha=0.25) +
	#geom_hline(yintercept=AUCmedian, color = "red", linewidth=0.8, linetype = 'dotted')+
	#annotate("text", x = 3, y=0.85, color="red", label = paste0("Median = ",round(AUCmedian,2)))+
  	ggtitle(" ")+
  	ylab("AUC")+
	#coord_cartesian(ylim = c(0, 0.9))+
	scale_y_continuous(breaks = c(0, 0.5, 0.8, 1))+
	theme_classic()+ 
	theme(
		panel.background = element_blank(),
		panel.grid = element_blank(),
		axis.title = element_text(colour = "black"),
		axis.text = element_text(colour = "black"),
		axis.text.x = element_text(angle = 38,hjust = 1),
		axis.title.x =element_blank(),
		legend.position = "none"
	)
ggsave(paste0(outputPath,"validation_scRNAseq_",cancer,"_",i,"_",sigName,"_ROC_comb.png"), p1, width = 6.5, height = 7.2, dpi=300, units = "cm", limitsize =FALSE)
write.csv(summaryTable,paste0(outputPath,"validation_scRNAseq_",cancer,"_",i,"_",sigName,"_ROC_comb.csv"), quote=FALSE)


	
