source("Secretome_s0_path.R")

dataset_cancers <- list.files(lambdaPath)
dataset_cancers <- sapply(strsplit(dataset_cancers,"_",fixed=T),function(x) return(paste0(x[1],"_",x[2])))
dataset_cancers <- unique(dataset_cancers)

#for(sig in sigs)
runCompare <- function(sig)
{
	smysmyMat <- data.frame()
	for(dataset_cancer in dataset_cancers)
	{
		dataset <- strsplit(dataset_cancer,"_")[[1]][1]
		cancer <- strsplit(dataset_cancer,"_")[[1]][2]
		
		if(!grepl("-",cancer,fixed=T))
		{
			cdata <- read.Xena(cancer)
			
			rownames(cdata) <- transferSymbol(rownames(cdata))
			cdata <- rm_duplicates(cdata)
			
			cdata <- filter.counts(cdata)
			
			cdata_T <- extract.samples(cdata,"T")
			cdata_N <- extract.samples(cdata,"N")
			
			if(ncol(cdata_N)>5)
			{
				cdata_T_minusBG <- cdata_T-rowMeans(cdata_N)
			}else{
				cdata_T_minusBG <- cdata_T-rowMeans(cdata_T)
			}
			
			cdata_T_minusBG <- as.matrix(cdata_T_minusBG)
		}else{
			cdata_T <- as.matrix(read.csv(paste0("../data/1_Creation/ICGC/",cancer,".seq.expression.gz"),sep="\t",row.names=1))	
						
			rownames(cdata_T) <- transferSymbol(rownames(cdata_T))
			cdata_T <- rm_duplicates(cdata_T)
			
			cdata_T <- filter.counts(cdata_T)
			
			cdata_T_minusBG <- cdata_T-rowMeans(cdata_T)
			cdata_T_minusBG <- as.matrix(cdata_T_minusBG)
		}	
	
	
		smyMat <- data.frame()
		for(lambda in lambdas)
		{
			load(paste0(lambdaPath,dataset_cancer,"_",sig,"_",lambda,".RData"))
			
			X <- as.matrix(res$zscore)				
			Y <- as.matrix(cdata_T_minusBG)
			
			olp <- intersect(rownames(X),rownames(Y))
			X_olp <- X[olp,,drop=F]
			Y_olp <- Y[olp,,drop=F]
		
			for(i in olp)
			{
				smyMat[paste0(lambda,"_",i),"lambda"] <- as.character(lambda)
				smyMat[paste0(lambda,"_",i),"r"] <- cor(X_olp[i,],Y_olp[i,])
			}
			
			smysmyMat[paste0(cancer,"_",lambda),"cancer"] <- cancer
			smysmyMat[paste0(cancer,"_",lambda),"lambda"] <- as.character(lambda)
			smysmyMat[paste0(cancer,"_",lambda),"r"] <- median(smyMat[smyMat[,"lambda"]==lambda,"r"])
		}
		
		
		smyMat[,"lambda"] <- factor(smyMat[,"lambda"],levels=(as.character(lambdas)))
		
		#library(ggplot2)
		#p2 <- ggplot(smyMat,aes(x = a, y = r, fill=a))+
		#	geom_hline(yintercept=0, color = "grey", linewidth=0.6, linetype="dashed")+
		#	geom_boxplot(outlier.shape=NA)+
		#	ggtitle(cancer)+
		#	xlab("Alpha")+
		#	ylab("Pearson r value")+
		#	theme_bw()+ 
		#	theme(
		#	  panel.grid = element_blank(),
		#	  panel.background = element_blank(),
		#	  axis.text = element_text(size=10,colour = "black"),
		#	  axis.title = element_text(size=10,colour = "black")
		#	) 
		#ggsave(paste0(QCSummaryPath,"lambda_summary/",cancer,"_new.png"), p2, width = 15, height = 10, dpi=200, units = "cm", limitsize = FALSE)

	}

	smysmyMat[,"lambda"] <- factor(smysmyMat[,"lambda"],levels=(as.character(lambdas)))
	write.csv(smysmyMat,paste0(lambdaSummaryPath,"lambda_compare_",sig,"_cancerType.csv"),quote=F)
	
	#smysmyMat <- read.csv(paste0(lambdaSummaryPath,"lambda_compare_",sig,"_cancerType.csv"),row.names=1)
	#smysmyMat[,"lambda"] <- factor(smysmyMat[,"lambda"],levels=(as.character(lambdas)))
	
	smysmysmyMat <- data.frame()
	for(i in unique(smysmyMat[,"lambda"]))
	{
		smysmysmyMat[i,"r_mean"] <- mean(smysmyMat[smysmyMat[,"lambda"]==i,3])
		smysmysmyMat[i,"r_median"] <- median(smysmyMat[smysmyMat[,"lambda"]==i,3])
	}
	smysmysmyMat[,"Max_r_mean"] <- smysmysmyMat[,"r_mean"]==max(smysmysmyMat[,"r_mean"])
	bestLambda <- rownames(smysmysmyMat)[smysmysmyMat[,"Max_r_mean"]==TRUE]
	
	write.csv(smysmysmyMat,paste0(lambdaSummaryPath,"lambda_compare_",sig,".csv"),quote=F)
	writeLines(bestLambda, paste0(finalSignaturesPath,sig,"_lambda.txt"))
	
	
	library(ggplot2)
	p2 <- ggplot(smysmyMat,aes(x = lambda, y = r, fill=lambda))+
		geom_jitter(size=0.8, alpha=0.3, width=0.2)+
		geom_boxplot(outlier.shape=NA, alpha=0.3, width=0.5)+
		ggtitle(sig)+
		#xlab("Lambda")+
		xlab(" ")+
		ylab("Pearson r value")+
		theme_bw()+ 
		theme(
		  panel.grid = element_blank(),
		  panel.background = element_blank(),
		  plot.title = element_text(hjust = 0.5),
		  axis.text = element_text(size=11,colour = "black"),
		  axis.title = element_text(colour = "black"),
		  #axis.text.x = element_text(angle = 20,hjust = 1,color = ifelse(smysmysmyMat[,1]==max(smysmysmyMat[,1]), "red","black")),
		  axis.text.x = element_text(color = ifelse(smysmysmyMat[,"Max_r_mean"]==TRUE, "red","black")),
		  legend.position="none"
		) 
	ggsave(paste0(lambdaSummaryPath,"lambda_compare_",sig,".png"), p2, width = 13, height = 6, dpi=200, units = "cm", limitsize = FALSE)
	
}

parallel::mclapply(sigNames, runCompare, mc.cores=5) 

