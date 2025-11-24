source("Secretome_s0_path.R")

#SWARM -t 2 -g 80 --time 10:30:00
#only one SKCM need long time
args = commandArgs(trailingOnly=TRUE)
st <- args[1]
sampleName <- args[2]

sigType <- "Secreted"

sts <- unique(meta[,"Study"])
#for(st in sts)
#{
	signaturePath.st <- paste0(signaturePath,st,"/")
	dir.create(signaturePath.st)
	
	QCPath.st <- paste0(QCPath,st,"/")
	dir.create(QCPath.st)
	
	QCFilterPath.st <- paste0(QCFilterPath,st,"/")
	dir.create(QCFilterPath.st)
	
	sampleNames <- meta[meta[,"Study"]==st,"Sample_Name"]
	#for(sampleName in sampleNames)
	#{	
		signaturePath.st.sample <- paste0(signaturePath.st,sampleName,"/")
		dir.create(signaturePath.st.sample)
		QCPath.st.sample <- paste0(QCPath.st,sampleName,"/")
		dir.create(QCPath.st.sample)
		QCFilterPath.st.sample <- paste0(QCFilterPath.st,sampleName,"/")
		dir.create(QCFilterPath.st.sample)
		
		
		# gating
		for(version in c("vst_free","vst","vst_condition_logUMI_cellType"))
		{
			st.matrix.data.vst <- as.matrix(read.table(gzfile(paste0(signaturePath.st.sample, version, ".tsv.gz")), sep="\t",check.names=F))
			W <- calWeights(colnames(st.matrix.data.vst), radius=200, sigma=100, diagAsZero=TRUE)
			
			
			# autocorrelation
			st.matrix.data.vst.SPs <- st.matrix.data.vst[rownames(st.matrix.data.vst)%in%SPs,]			
			smry <- spatialAutoCorrelation(st.matrix.data.vst.SPs, W, permuteNo=1000)
			write.csv(smry,paste0(QCPath.st.sample,"AutoCorrelation_",version,".csv"),quote=F)
			
			
			# gating
			sigmat <- spatialCrossCorrelation(st.matrix.data.vst, W)
	
			# TCGA
			cancers <- list.files(TCGAPath)
			cancers <- sapply(strsplit(cancers,".",fixed=T),function(x) return(x[1]))
			cancers1 <- setdiff(unique(cancers), "LAML") # non solid
			
			for(cancer in cancers1)
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
				
				smyMat <- cor.act.exp(X=sigmat, Y=cdata_T_minusBG)
				
				write.csv(smyMat,paste0(QCPath.st.sample,"TCGA_",cancer,"_",version,".csv"),quote=F)
				
				# for fig1e-f and fig2a
				# random signature
				if(version=="vst"&st=="SKCM_2022_Sudmeier"&sampleName=="pt16")
				{
					set.seed(123)
					combPermute <- data.frame()
					for(i in 1:10)
					{
						sigmat_random <- sigmat
						rownames(sigmat_random) <- sample(rownames(sigmat))
						smyMat2 <- cor.act.exp(X=sigmat_random, Y=cdata_T_minusBG)
						
						combPermute[names(smyMat2),i] <- smyMat2
					}
					combPermute <- apply(combPermute,1,mean)
					write.csv(combPermute,paste0(QCPath.st.sample,"TCGA_",cancer,"_random.csv"),quote=F)
					
					if(cancer=="SKCM")
					{
						gene <- "TGFB1"
						
						X = sigmat
						Y = cdata_T_minusBG
						olp <- intersect(rownames(X),rownames(Y))
						X_olp <- X[olp,,drop=F]
						Y_olp <- Y[olp,,drop=F]
						
						# signature score
						cc_corr_r <- WGCNA::cor(X_olp,Y_olp,use="pairwise.complete.obs")	
						
						fg.df <- data.frame(x=Y_olp[gene,],y=cc_corr_r[gene,])
						write.csv(fg.df,paste0(QCPath.st.sample,"TCGA_",cancer,"_",gene,"_expr_vs_signatureSoce.csv"),quote=F)
					}
				}
			}
			
	
			# ICGC
			allFiles <- list.files(ICGCPath)
			allFiles <- allFiles[grepl("seq.expression",allFiles)]
			cancers <- gsub(".seq.expression.gz","",allFiles)
			cancers <- setdiff(cancers,"CLLE-ES") # non solid
			cancers <- setdiff(cancers,c("BPLL-FR","PAEN-AU","PRAD-FR")) # sample <40
			cancers2 <- cancers
	
			for(cancer in cancers2)
			{
				cdata_T <- as.matrix(read.csv(paste0(ICGCPath,cancer,".seq.expression.gz"),sep="\t",row.names=1))	
				#print(paste0(cancer," ",ncol(cdata_T) )) }
							
				rownames(cdata_T) <- transferSymbol(rownames(cdata_T))
				cdata_T <- rm_duplicates(cdata_T)
				
				cdata_T <- filter.counts(cdata_T)
				
				cdata_T_minusBG <- cdata_T-rowMeans(cdata_T)
				cdata_T_minusBG <- as.matrix(cdata_T_minusBG)
				
				smyMat <- cor.act.exp(X=sigmat, Y=cdata_T_minusBG)
				
				write.csv(smyMat,paste0(QCPath.st.sample,"ICGC_",cancer,"_",version,".csv"),quote=F)
				
				# for fig1e-f and fig2a
				# random signature
				if(version=="vst"&st=="SKCM_2022_Sudmeier"&sampleName=="pt16")
				{
					set.seed(123)
					combPermute <- data.frame()
					for(i in 1:10)
					{
						sigmat_random <- sigmat
						rownames(sigmat_random) <- sample(rownames(sigmat))
						smyMat2 <- cor.act.exp(X=sigmat_random, Y=cdata_T_minusBG)
						
						combPermute[names(smyMat2),i] <- smyMat2
					}
					combPermute <- apply(combPermute,1,mean)
					write.csv(combPermute,paste0(QCPath.st.sample,"ICGC_",cancer,"_random.csv"),quote=F)
				}
			}
		
			
			
			
			combStat <- list() # one shot
			combStat[["old"]] <- data.frame()
			combStat[["new_0.05"]] <- data.frame()
			combStat[["new_0.1"]] <- data.frame()
			combStat[["new_0.15"]] <- data.frame()
			combStat[["new_0.2"]] <- data.frame()
			combStat[["new_0.25"]] <- data.frame()
			
			comb <- data.frame() # gene * cancer (signature score)
			combTest <- data.frame() # wil.cox test secreted vs intra
			combComp <- data.frame() # cancer * gating strategy (filtered signature count)
			combComp_r <- data.frame() # cancer * gating strategy (r cutoff value)
			combFDR <- data.frame() # gene * cancer (sorted value)
			
			cancers <- c(paste0("TCGA_",cancers1), paste0("ICGC_",cancers2))
			
			for(cancer in cancers)
			{
				smyMat <- read.csv(paste0(QCPath.st.sample,cancer,"_",version,".csv"),row.names=1)
				smyMat <- smyMat[!is.na(smyMat[,1]),,drop=F]
				comb[rownames(smyMat),cancer] <- smyMat[,1]
				
				smyMat <- cbind(smyMat,Group="NA")
				smyMat[rownames(smyMat)%in%IPs,"Group"] <- "Intracellular"
				smyMat[rownames(smyMat)%in%MPs,"Group"] <- "Membrane"
				smyMat[rownames(smyMat)%in%SPs,"Group"] <- "Secreted"
				smyMat <- smyMat[smyMat[,"Group"]!="NA",]
				
				x1 <- smyMat[smyMat[,"Group"]==sigType,1]
				names(x1) <- rownames(smyMat)[smyMat[,"Group"]==sigType]
				x2 <- smyMat[smyMat[,"Group"]=="Intracellular",1]
				names(x2) <- rownames(smyMat)[smyMat[,"Group"]=="Intracellular"]
				
				pv <- signif(wilcox.test(x1, x2, alternative="greater")$p.value,2)
				
				combTest[cancer,"Var2"] <- cancer
				combTest[cancer,"bar"] <- quantile(x2,0.95)
				combTest[cancer,"Significant"] <- ifelse(pv<0.01,"Yes","No")
				
				
				# for old gating strategy
				#if(combTest[cancer,"Significant"]=="Yes")
				#{
					filterSPs <- names(x1)[x1>combTest[cancer,"bar"]]
					
					if(length(filterSPs)>0) combStat[["old"]][filterSPs,cancer] <- 1
					
					combComp[cancer,"old"] <- length(filterSPs)
					combComp_r[cancer,"old"] <- NA
				#}
				
				# for new gating strategy
				x2_sorted <- sort(x2,decreasing=T)
				tempMin <- 1
				for(i in length(x2_sorted):1 )
				{
					X <- sum( x2 >= x2_sorted[[i]] )/length(x2)
					Y <- sum( x1 >= x2_sorted[[i]] )/length(x1)
					combFDR[i,cancer] <- min(X/Y,tempMin)
					tempMin <- min(X/Y,tempMin)
				}
				
				
				for(ratio in c(0.05 * 1:5))
				{
					for(i in length(x2_sorted):1 )
					{
						if(combFDR[i,cancer] <= ratio)
						{
							combComp[cancer,paste0("new_",ratio)] <- sum(x1 > x2_sorted[i])
							combComp_r[cancer,paste0("new_",ratio)] <- x2_sorted[i]
	
							combStat[[paste0("new_",ratio)]][names(x1)[x1>x2_sorted[i]],cancer] <- 1
							break
						}
					}
				}
				
			} # cancer
			
			
			for(i in names(combStat) )
			{
				tempMat <- combStat[[i]]
				tempMat[is.na(tempMat)] <- 0
				write.csv(tempMat,paste0(QCFilterPath.st.sample,sigType,"_filterBar_",version,"_",i,".csv"), quote=F)
			}
			
			write.csv(combTest,paste0(QCFilterPath.st.sample,sigType,"_filterBar_",version,"_Test.csv"), quote=F)
			write.csv(combComp,paste0(QCFilterPath.st.sample,sigType,"_filterBar_",version,"_comp.csv"), quote=F)
			write.csv(combComp_r,paste0(QCFilterPath.st.sample,sigType,"_filterBar_",version,"_comp_r.csv"), quote=F)
			write.csv(combFDR,paste0(QCFilterPath.st.sample,sigType,"_filterBar_",version,"_FDR.csv"), quote=F)
			
			
			# draw figure gating figure
			comb.m <- reshape2::melt(as.matrix(comb))
			comb.m <- comb.m[!is.na(comb.m[,3]),]
			comb.m <- cbind(comb.m,Group="NA")
			
			comb.m[comb.m[,1]%in%IPs,"Group"] <- "Intracellular"
			comb.m[comb.m[,1]%in%MPs,"Group"] <- "Membrane"
			comb.m[comb.m[,1]%in%SPs,"Group"] <- "Secreted"
			
			#comb.m <- comb.m[!comb.m[,"Group"]%in%c("Membrane","NA"),]		
			comb.m <- comb.m[!comb.m[,"Group"]%in%c("NA"),]		
			
			combTest[,1] <- factor(combTest[,1],levels=unique(combTest[,1]))
			
			library(ggplot2)
			library(ggsignif)
			p2 <- ggplot(comb.m,aes(x = Group, y = value))+
				geom_hline(yintercept=0, color = "grey", linewidth=0.6, linetype="dashed")+
				geom_violin(aes(group=Group, fill=Group),trim=FALSE)+
				geom_boxplot(width=0.15,outlier.shape=NA)+
				geom_signif( comparisons = list(c(sigType, "Intracellular")) ,y_position = c(0.9), test = "wilcox.test", test.args ="greater")+
				geom_hline(data=combTest, aes(yintercept = bar, colour=Significant),linewidth=1.5, linetype="dashed")+
				xlab("Group")+
				ylab("Correlation between activity and expression")+
				theme_bw()+ 
				theme(
				  panel.grid = element_blank(),
				  panel.background = element_blank(),
				  axis.text = element_text(size=10,colour = "black"),
				  axis.text.x = element_blank(),
				  axis.title = element_text(size=10,colour = "black"),
				  legend.position="right",
				  strip.text = element_text(size = 10)
				) + facet_wrap(~ Var2, ncol = 11 ) #scales = "free"
				
			ggsave(paste0(QCFilterPath.st.sample,sigType,"_filterBar_",version,".png"), p2, width = 50, height = 22, dpi=200, units = "cm",limitsize = FALSE)
			
			
			
			comb.m <- reshape2::melt(as.matrix(combFDR))
			comb.m <- comb.m[!is.na(comb.m[,3]),]
			
			library(ggplot2)
			p2 <- ggplot(comb.m,aes(x = Var1, y = value))+
				geom_point(size=0.2)+
				geom_hline(yintercept=0.25, color = "grey", linewidth=0.6, linetype="dashed")+
				geom_hline(yintercept=0.20, color = "grey", linewidth=0.6, linetype="dashed")+
				geom_hline(yintercept=0.15, color = "grey", linewidth=0.6, linetype="dashed")+
				geom_hline(yintercept=0.10, color = "grey", linewidth=0.6, linetype="dashed")+
				geom_hline(yintercept=0.05, color = "grey", linewidth=0.6, linetype="dashed")+
				xlab("Rank")+
				ylab("Q")+
				theme_bw()+ 
				theme(
				  panel.grid = element_blank(),
				  panel.background = element_blank(),
				  axis.text = element_text(size=10,colour = "black"),
				  axis.title = element_text(size=10,colour = "black"),
				  legend.position="none",
				  strip.text = element_text(size = 10)
				) + facet_wrap(~ Var2, ncol = 11 ) #scales = "free"
				
			ggsave(paste0(QCFilterPath.st.sample,sigType,"_filterBar_",version,"_FDR.png"), p2, width = 50, height = 22, dpi=200, units = "cm",limitsize = FALSE)
		
		}# version
		
	#}
#}

