source("Secretome_s0_path.R")

cancer <- "ProCan-DepMapSanger"


inputPath <- paste0(dataValPath,cancer,"/")
outputPath <- paste0(validationPath,cancer,"/")


# read proteomics
Protein <- read.csv(paste0(inputPath,"Proteomics_20250211/Protein_matrix_averaged_20250211.tsv"),check.names=F,sep="\t",comment.char = "#")
Protein_anno <- Protein[,1:2]
Protein <- Protein[,3:ncol(Protein)]
rownames(Protein) <- Protein_anno[,1]
Protein <- t(Protein)

Protein <- Protein[rowSums(!is.na(Protein))>ncol(Protein)*0.2,,drop=F]

rownames(Protein) <- transferSymbol(rownames(Protein))
Protein <- rm_duplicates(Protein)

# read pathway
load(paste0(outputPath,"/",cancer,"_Pathway.RData"))

secreted_pathway <- read.csv(paste0(dataValPath,"Reactome/secreted_to_pathway_relation.txt"),sep="\t",row.names=2)


# draw AUC
matched <- FALSE # CPTAC matched tumor and normal samples
for(responseType in c("proteomics","pathway"))
{
	if(responseType=="pathway") sigNames <- setdiff(sigNames,c("NicheNet.v1","NicheNet.v2"))
	
	fg.df <- data.frame()
	for(sigName in sigNames)
	{
		load(paste0(outputPath,cancer,"_",sigName,".RData"))
		Act <- as.matrix(res$zscore)
		
		if(responseType=="proteomics")
		{
			predictor <- Act
			
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
	} # sigName
		
	
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
				fg.df_gene[paste0(gene,sig),"value"] <- median(fg.df_sub[,aa],na.rm=T)
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
			ggsave(paste0(outputPath,"validation_ProCan_compare_all_",responseType,"_",aa,".png"), p2, width = 5, height = 7.8, dpi=200, units = "cm", limitsize = FALSE)
		}else{
			p2 <- p2+
			annotate("text", x = 4, y=0.02, label = paste0(sig_stat[rownames(sigOrder)[4]] ))+
			annotate("text", x = 5, y=0.02, label = paste0(sig_stat[rownames(sigOrder)[5]] ))
			
			ggsave(paste0(outputPath,"validation_ProCan_compare_all_",responseType,"_",aa,".png"), p2, width = 7, height = 7.8, dpi=200, units = "cm", limitsize = FALSE)
		}
		
		write.csv(fg.df_gene.d, paste0(outputPath,"validation_ProCan_compare_all_",responseType,"_",aa,".csv"), quote=FALSE)
	
		#x1 <- fg.df_gene.d[,"SecAct"]
		#x2 <- fg.df_gene.d[,"CytoSig"]
		#x1 <- x1[!is.na(x1)]
		#x2 <- x2[!is.na(x2)]
		#print(wilcox.test(x1,x2))
	} # aa

} # responseType
