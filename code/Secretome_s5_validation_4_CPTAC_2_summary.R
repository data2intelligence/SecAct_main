source("Secretome_s0_path.R")


inputPath <- paste0(dataValPath,"CPTAC/")
outputPath <- paste0(validationPath,"CPTAC/")


# CPTAC
cancers <- c("BRCA", "COAD", "GBM", "HNSC", "KIRC", "LUAD", "LUSC", "OV", "PDAC", "UCEC")


# Ensemble ID -> gene symbol
geneAnno <- read.csv(paste0(inputPath,"Proteome_BCM_GENCODE_v34_harmonized_v1/README/Gene_annotation_and_representable_isoform_mapping_table.txt"),sep="\t",check.names=F)
geneAnno[geneAnno[,1]=="ENSG00000284024.2",4] <- "MSANTD7"
geneAnno <- geneAnno[,c(1,4)]
geneAnno <- geneAnno[!duplicated(geneAnno[,1]),]
rownames(geneAnno) <- geneAnno[,1]

secreted_pathway <- read.csv(paste0(dataValPath,"Reactome/secreted_to_pathway_relation.txt"),sep="\t",row.names=2)

# compare different signatures
matched <- FALSE # CPTAC matched tumor and normal samples
for(responseType in c("proteomics","pathway"))
{
	if(responseType=="pathway") sigNames <- setdiff(sigNames,c("NicheNet.v1","NicheNet.v2"))
	
	fg.df <- data.frame()
	for(sigName in sigNames)
	{
		for(cancer in cancers)
		{
			load(paste0(outputPath,"/",cancer,"_",sigName,".RData"))
			Act <- as.matrix(res$zscore)
			
			# prepare predictor ~ response
			if(responseType=="proteomics")
			{
				predictor <- Act
				
				# proteomics
				Protein <- as.matrix(read.csv(paste0(inputPath,"Proteome_BCM_GENCODE_v34_harmonized_v1/",cancer,"_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt"),sep="\t",row.names=1,check.names=F))
				
				rownames(Protein) <- geneAnno[rownames(Protein),2]
				rownames(Protein) <- transferSymbol(rownames(Protein))
				Protein <- rm_duplicates(Protein)
				Protein <- Protein[rowSums(!is.na(Protein))>ncol(Protein)*0.2,,drop=F]
				
				response <- Protein
			}else{
				
				library(dplyr)
				library(stringr)
				# Convert to list
				gene_sets <- secreted_pathway %>%
				  mutate(genes = str_split(geneID, "/")) %>%
				  split(.$row.names) %>%
				  lapply(function(x) x$genes[[1]])
				
				
				avg_Act <- sapply(gene_sets, function(genes) {
				  overlap <- intersect(genes, rownames(Act))  # keep only available genes
				  if(length(overlap) == 0) {
				    return(rep(NA, ncol(Act)))  # if no genes found
				  }
				  colMeans(Act[overlap, , drop=FALSE])  # average across genes
				})
				
				# Convert to data.frame
				avg_Act <- as.data.frame(avg_Act)
				rownames(avg_Act) <- colnames(Act)  # samples as rows
				
				avg_Act <- t(avg_Act)
				avg_Act <- avg_Act[!is.na(avg_Act[,1]),]
				
				predictor <- avg_Act
				

				# pathway
				load(paste0(outputPath,"/",cancer,"_Pathway.RData"))
				response <- gsva_scores
			}
			
			olp_c <- intersect(colnames(predictor),colnames(response))
			olp_r <- intersect(rownames(predictor),rownames(response))
			
			predictor_olp <- predictor[olp_r,olp_c,drop=F]
			response_olp <- response[olp_r,olp_c,drop=F]
			
			genes <- olp_r
			genes_alt <- c()
			corVec <- c()
			aucVec <- c()
			auprcVec <- c()
		
			for(gene in genes)
			{
				x <- unlist(predictor_olp[gene,])
				y <- unlist(response_olp[gene,])
				
				x <- x[!is.na(y)]
				y <- y[!is.na(y)]
				
				if(matched==TRUE)
				{
					if(! (sum(y>0) >= 5 & sum(y<0) >= 5) ) next
				}else{
					if(length(y) < 10 ) next
				}
				
				genes_alt <- c(genes_alt, gene)
				rv <- cor(x,y,use="pairwise.complete.obs")
				corVec <- c(corVec,rv)
				
				if(matched==TRUE)
				{
					y_alt <- as.numeric(y >= 0 ) 
				}else{
					y_alt <- as.numeric(y >= quantile(y)[4] ) ################ =?
				}
				
				library(ROCR)
				predM <- prediction(x, y_alt)
				roc <-  performance(predM, measure = "auc")
				roc2 = performance(predM, measure = "aucpr")
				
				av <- unlist(roc@ y.values)
				apv <- unlist(roc2@ y.values)
				
				aucVec <- c(aucVec,av)
				auprcVec <- c(auprcVec,apv)
				
				
				if(sigName=="SecAct"&cancer=="BRCA"&gene=="R-HSA-6783783")
				{
					roc.df <- data.frame(expr=as.character(y_alt),act=x)
					
					x1 <- roc.df[roc.df[,1]=="1",2]
					x2 <- roc.df[roc.df[,1]=="0",2]
					pv <- signif(wilcox.test(x1, x2)$p.value,2)
					
					library(ggplot2)
					p2 <- ggplot(roc.df,aes(x=expr,y=act)) + 
						geom_jitter(aes(fill=expr),colour="black",width=0.2,shape = 21,alpha=0.6)+
						scale_fill_manual(values=c("#00BFC4","#F8766D"))+
						scale_x_discrete(labels= c("Low","High"))+
						annotate("text", x = 2, y=-15, label = paste0("p = ",pv))+
						ylab("IL10 SecAct Activity")+
						xlab("IL10 Pathway Activity")+
						theme_classic()+ 
						theme(
							panel.grid = element_blank(),
					  		panel.background = element_blank(),
					  		plot.title = element_text(hjust = 0.5),
    						axis.title = element_text(colour = "black"),
							axis.text = element_text(colour = "black"),
							legend.position="none"
						)
					ggsave(paste0(outputPath,"validation_CPTAC_single_jitter.png"), p2, width = 6.8, height = 6.5, dpi=200, units = "cm",limitsize = FALSE)
					write.csv(roc.df,paste0(outputPath,"validation_CPTAC_single_jitter.csv"),quote=F)

					
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
					ggsave(paste0(outputPath,"validation_CPTAC_single_ROC.png"), p2, width = 6.8, height = 6.5, dpi=200, units = "cm",limitsize = FALSE)
		
				}
				
			}
			
			fg.df[paste0(sigName,cancer,genes_alt),"sig"] <- sigName
			fg.df[paste0(sigName,cancer,genes_alt),"cancer"] <- cancer
			fg.df[paste0(sigName,cancer,genes_alt),"gene"] <- genes_alt
			fg.df[paste0(sigName,cancer,genes_alt),"cor"] <- corVec
			fg.df[paste0(sigName,cancer,genes_alt),"AUC"] <- aucVec
			fg.df[paste0(sigName,cancer,genes_alt),"AUPRC"] <- auprcVec
		
		} #cancer
		
		
		
		if(responseType=="pathway"&sigName=="SecAct")
		{
			fg.df_secact <- fg.df[fg.df[,"sig"]==sigName,]
			fg.df_secact[,"cancer"] <- CPTAC_rename[match(fg.df_secact[,"cancer"],CPTAC_rename[,1]),2]
			
			p1 <- ggplot(fg.df_secact, aes(x = cancer, y = AUC))+
				geom_hline(yintercept=0.5, color = "grey1", linewidth=0.75, linetype = 'dotted')+
				geom_violin(aes(group=cancer), fill =sigColors[sigName], trim=FALSE, alpha=0.8, width=0.8)+
				geom_boxplot(width=0.15,outlier.shape=NA)+
				annotate("text", x = 5.5, y=0.08, label = paste0("n = ", length(unique(fg.df_secact[,"gene"])) ))+
				ylim(0,1)+
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
			
			ggsave(paste0(outputPath,"validation_CPTAC.png"), p1, width = 14, height = 8.8, dpi=200, units = "cm", limitsize = FALSE)
			write.csv(fg.df_secact,paste0(outputPath,"validation_CPTAC.csv"),quote=F)

		}

		
	}#sigName

	
	
	
	
	
	
	
	
	for(aa in c("AUC","AUPRC"))
	{
		fg.df_gene <- data.frame()
		for(gene in unique(fg.df[,"gene"]) )
		{
			for(sig in unique(fg.df[,"sig"]) )
			{
				fg.df_sub <- fg.df[fg.df[,"gene"]==gene&fg.df[,"sig"]==sig,]
				
				if(nrow(fg.df_sub)==0) next
				
				fg.df_gene[paste0(gene,sig),"gene"] <- gene
				fg.df_gene[paste0(gene,sig),"sig"] <- sig
				fg.df_gene[paste0(gene,sig),"value"] <- median(fg.df_sub[,aa],na.rm=TRUE)
			}
		}
		
		if(responseType=="pathway") fg.df_gene[,"gene"] <- secreted_pathway[fg.df_gene[,"gene"],"Description"]
		
		fg.df_gene.d = reshape2::dcast( fg.df_gene , sig~gene )
		rownames(fg.df_gene.d) <- fg.df_gene.d[,1]
		fg.df_gene.d <- t(fg.df_gene.d[,-1])

		
		sigOrder <- data.frame()
		for(i in unique(fg.df_gene[,"sig"]))
		{
			sigOrder[i,"r_mean"] <- median(fg.df_gene[fg.df_gene[,"sig"]==i,"value"])
		}
		sigOrder <- sigOrder[order(sigOrder[,"r_mean"],decreasing=T),,drop=F]
		
		fg.df_gene[,2] <- factor(fg.df_gene[,2], levels=rownames(sigOrder))
		
		sig_stat <- table(fg.df_gene[,"sig"])
		
		library(ggplot2)
		p2 <- ggplot(fg.df_gene, aes(x = sig, y = value))+ 
			geom_hline(yintercept=ifelse(aa=="AUC",0.5,0.25), color = "grey1", linewidth=0.8, linetype = 'dotted')+
			geom_violin(aes(group=sig, fill=sig),width=0.8,trim=FALSE)+
			geom_boxplot(width=0.15,outlier.shape=NA)+
			scale_fill_manual( values=sigColors[rownames(sigOrder)] )+
			annotate("text", x = 1, y=0.1, label = paste0("n = "))+
			annotate("text", x = 1, y=0.02, label = paste0(sig_stat[rownames(sigOrder)[1]] ))+
			annotate("text", x = 2, y=0.02, label = paste0(sig_stat[rownames(sigOrder)[2]] ))+
			annotate("text", x = 3, y=0.02, label = paste0(sig_stat[rownames(sigOrder)[3]] ))+
			ylim(0,1)+
			ggtitle(" ")+
			xlab(" ")+
			ylab(aa)+
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
		
		
		if(responseType=="pathway")
		{
			ggsave(paste0(outputPath,"validation_CPTAC_compare_all_",responseType,"_",aa,".png"), p2, width = 5, height = 7.8, dpi=200, units = "cm", limitsize = FALSE)
		}else{
			p2 <- p2+
			annotate("text", x = 4, y=0.02, label = paste0(sig_stat[rownames(sigOrder)[4]] ))+
			annotate("text", x = 5, y=0.02, label = paste0(sig_stat[rownames(sigOrder)[5]] ))
			
			ggsave(paste0(outputPath,"validation_CPTAC_compare_all_",responseType,"_",aa,".png"), p2, width = 7, height = 7.8, dpi=200, units = "cm", limitsize = FALSE)
		}
		
		write.csv(fg.df_gene.d, paste0(outputPath,"validation_CPTAC_compare_all_",responseType,"_",aa,".csv"), quote=FALSE)
	
		#x1 <- fg.df_gene.d[,"SecAct"]
		#x2 <- fg.df_gene.d[,"CytoSig"]
		#x1 <- x1[!is.na(x1)]
		#x2 <- x2[!is.na(x2)]
		#print(wilcox.test(x1,x2))
	}#aa
	
}#responseType




####################################
# compare Leave-One-Out signatures #
####################################
stat <- read.csv(paste0(preprocessDataSummaryPath,"visium_stat_fullname.csv"))

fg.df <- data.frame()
for(cancer in CPTAC_rename[,1])
{
	fg.df[cancer,"cancer"] <- cancer
	
	if(cancer=="COAD")
	{
		fg.df[cancer,"no_LOO"] <- stat[grepl("CRC",stat[,1]),2]
	}else{
		fg.df[cancer,"no_LOO"] <- stat[grepl(cancer,stat[,1]),2]
	}
}

fg.df[,2] <- sum(stat[,2]) - fg.df[,2]

fg.df[,"cancer"] <- CPTAC_rename[match(fg.df[,"cancer"],CPTAC_rename[,1]),2]

fg.df <- fg.df[order(fg.df[,1]),]

fg.df <- rbind(c("Pan-cancer", sum(stat[,2])), fg.df)
fg.df[,2] <- as.numeric(fg.df[,2])

fg.df[,"Signature"] <- fg.df[,1]
fg.df[fg.df[,"Signature"]!="Pan-cancer","Signature"] <- "Leave-one-cancer-type-out"


fg.df[,1] <- factor(fg.df[,1], levels=fg.df[,1])
fg.df[,"Signature"] <- factor(fg.df[,"Signature"], levels=c("Pan-cancer","Leave-one-cancer-type-out"))
	 

p1 <- ggplot(fg.df, aes(x = cancer, y = no_LOO, label = no_LOO, fill=Signature))+
  		geom_bar(stat="identity", color="white", width=0.8, alpha=0.8) +
		scale_fill_manual(values=c("#9163b6","#e2975d"))+
		geom_text(vjust=-0.2)+
		ylab("# ST-samples")+
		ylim(0,1320)+
		theme_classic()+ 
		theme(
		  panel.grid = element_blank(),
		  panel.background = element_blank(),
		  plot.title = element_blank(),
		  axis.text = element_text(colour = "black"),
		  axis.title = element_text(colour = "black"),
		  axis.text.x = element_text(angle = 30, hjust = 1),
		  axis.title.x = element_blank(),
		  legend.position="top"
		)
	
ggsave(paste0(outputPath,"validation_LOO_sample_number.png"), p1, width = 15, height = 7.5, dpi=300, units = "cm", limitsize = FALSE)
write.csv(fg.df, paste0(outputPath,"validation_LOO_sample_number.csv"), quote=FALSE)
	



sigNames <- c("SecAct","SecAct.LeaveOneOut")

matched <- FALSE # CPTAC matched tumor and normal samples
for(responseType in c("proteomics","pathway"))
{
	if(responseType=="pathway") sigNames <- setdiff(sigNames,c("NicheNet.v1","NicheNet.v2"))
	
	fg.df <- data.frame()
	for(sigName in sigNames)
	{
		for(cancer in cancers)
		{
			load(paste0(outputPath,"/",cancer,"_",sigName,".RData"))
			Act <- as.matrix(res$zscore)
			
			# prepare predictor ~ response
			if(responseType=="proteomics")
			{
				predictor <- Act
				
				# proteomics
				Protein <- as.matrix(read.csv(paste0(inputPath,"Proteome_BCM_GENCODE_v34_harmonized_v1/",cancer,"_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt"),sep="\t",row.names=1,check.names=F))
				
				rownames(Protein) <- geneAnno[rownames(Protein),2]
				rownames(Protein) <- transferSymbol(rownames(Protein))
				Protein <- rm_duplicates(Protein)
				Protein <- Protein[rowSums(!is.na(Protein))>ncol(Protein)*0.2,,drop=F]
				
				response <- Protein
			}else{
				
				library(dplyr)
				library(stringr)
				# Convert to list
				gene_sets <- secreted_pathway %>%
				  mutate(genes = str_split(geneID, "/")) %>%
				  split(.$row.names) %>%
				  lapply(function(x) x$genes[[1]])
				
				
				avg_Act <- sapply(gene_sets, function(genes) {
				  overlap <- intersect(genes, rownames(Act))  # keep only available genes
				  if(length(overlap) == 0) {
				    return(rep(NA, ncol(Act)))  # if no genes found
				  }
				  colMeans(Act[overlap, , drop=FALSE])  # average across genes
				})
				
				# Convert to data.frame
				avg_Act <- as.data.frame(avg_Act)
				rownames(avg_Act) <- colnames(Act)  # samples as rows
				
				avg_Act <- t(avg_Act)
				avg_Act <- avg_Act[!is.na(avg_Act[,1]),]
				
				predictor <- avg_Act
				

				# pathway
				load(paste0(outputPath,"/",cancer,"_Pathway.RData"))
				response <- gsva_scores
			}
			
			olp_c <- intersect(colnames(predictor),colnames(response))
			olp_r <- intersect(rownames(predictor),rownames(response))
			
			predictor_olp <- predictor[olp_r,olp_c,drop=F]
			response_olp <- response[olp_r,olp_c,drop=F]
			
			genes <- olp_r
			genes_alt <- c()
			corVec <- c()
			aucVec <- c()
			auprcVec <- c()
		
			for(gene in genes)
			{
				x <- unlist(predictor_olp[gene,])
				y <- unlist(response_olp[gene,])
				
				x <- x[!is.na(y)]
				y <- y[!is.na(y)]
				
				if(matched==TRUE)
				{
					if(! (sum(y>0) >= 5 & sum(y<0) >= 5) ) next
				}else{
					if(length(y) < 10 ) next
				}
				
				genes_alt <- c(genes_alt, gene)
				rv <- cor(x,y,use="pairwise.complete.obs")
				corVec <- c(corVec,rv)
				
				if(matched==TRUE)
				{
					y_alt <- as.numeric(y >= 0 ) 
				}else{
					y_alt <- as.numeric(y >= quantile(y)[4] ) ################ =?
				}
				
				library(ROCR)
				predM <- prediction(x, y_alt)
				roc <-  performance(predM, measure = "auc")
				roc2 = performance(predM, measure = "aucpr")
				
				av <- unlist(roc@ y.values)
				apv <- unlist(roc2@ y.values)
				
				aucVec <- c(aucVec,av)
				auprcVec <- c(auprcVec,apv)
				
			}
			
			fg.df[paste0(sigName,cancer,genes_alt),"sig"] <- sigName
			fg.df[paste0(sigName,cancer,genes_alt),"cancer"] <- cancer
			fg.df[paste0(sigName,cancer,genes_alt),"gene"] <- genes_alt
			fg.df[paste0(sigName,cancer,genes_alt),"cor"] <- corVec
			fg.df[paste0(sigName,cancer,genes_alt),"AUC"] <- aucVec
			fg.df[paste0(sigName,cancer,genes_alt),"AUPRC"] <- auprcVec
		
		} #cancer
		
	}#sigName

	
	for(aa in c("AUC","AUPRC"))
	{
		fg.df_gene <- cbind(fg.df[,1:3], value=fg.df[,aa])
		fg.df_gene[,"cancer"] <- CPTAC_rename[match(fg.df_gene[,"cancer"],CPTAC_rename[,1]),2]
		fg.df_gene[fg.df_gene[,"sig"]!="SecAct","sig"] <- "Leave-one-cancer-type-out"
		fg.df_gene[fg.df_gene[,"sig"]=="SecAct","sig"] <- "Pan-cancer"
		fg.df_gene[,"sig"] <- factor(fg.df_gene[,"sig"], levels=c("Pan-cancer", "Leave-one-cancer-type-out") )
		
		library(ggplot2)
		library(ggsignif)
		p2 <- ggplot(fg.df_gene, aes(x = sig, y = value))+ 
			geom_hline(yintercept=ifelse(aa=="AUC",0.5,0.25), color = "grey1", linewidth=0.8, linetype = 'dotted')+
			geom_violin(aes(group=sig, fill=sig),width=0.8,trim=FALSE,alpha=0.8)+
			scale_fill_manual(values=c("#9163b6","#e2975d"))+
			geom_boxplot(width=0.15,outlier.shape=NA)+
			geom_signif( comparisons = list(c("Pan-cancer","Leave-one-cancer-type-out")) ,y_position = c(0.85), test = "wilcox.test")+
			#scale_fill_manual( values=sigColors[rownames(sigOrder)] )+
			#annotate("text", x = 1, y=0.1, label = paste0("n = "))+
			#annotate("text", x = 1, y=0.02, label = paste0(sig_stat[rownames(sigOrder)[1]] ))+
			#annotate("text", x = 2, y=0.02, label = paste0(sig_stat[rownames(sigOrder)[2]] ))+
			#annotate("text", x = 3, y=0.02, label = paste0(sig_stat[rownames(sigOrder)[3]] ))+
			ylim(0,1)+
			ggtitle(responseType)+
			xlab(" ")+
			ylab(aa)+
			theme_bw()+ 
			theme(
			  panel.grid = element_blank(),
			  panel.background = element_blank(),
			  plot.title = element_text(hjust = 0.5),
			  axis.text = element_text(colour = "black"),
			  axis.title = element_text(colour = "black"),
			  axis.text.x = element_blank(),
			  axis.ticks.x = element_blank(),
			  strip.background = element_rect(fill = "white"),
			  legend.position="none"
			)+ facet_wrap(~ cancer, nrow =1)
		
		ggsave(paste0(outputPath,"validation_CPTAC_compare_all_",responseType,"_",aa,"_LeaveOneOut.png"), p2, width = 30, height = 5, dpi=200, units = "cm", limitsize = FALSE)
		
		write.csv(fg.df_gene, paste0(outputPath,"validation_CPTAC_compare_all_",responseType,"_",aa,"_LeaveOneOut.csv"), quote=FALSE)
		
	}#aa
	
}#responseType
