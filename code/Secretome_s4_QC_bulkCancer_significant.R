source("Secretome_s0_path.R")
source("../../SpaCE/code/ifun.R")
source("importHPA.R")

dataPath <- "/data/rub2/data/TCGA_firebrowse/mRNAseq_gene_exp_symbol"
cancers <- list.files(dataPath)
cancers <- setdiff(cancers, "LAML")
cancers1 <- cancers

dataPath <- "/data/rub2/data/ICGC/"
allFiles <- list.files(dataPath)
allFiles <- allFiles[grepl("seq.expression",allFiles)]
cancers <- gsub(".seq.expression.gz","",allFiles)
cancers <- setdiff(cancers,"CLLE-ES") 
cancers <- setdiff(cancers,c("BPLL-FR","PAEN-AU","PRAD-FR")) # sample <50
cancers2 <- cancers
#CLLE-ES	Chronic Lymphocytic Leukemia - ES	Spain
#MALY-DE	Malignant Lymphoma - DE	Germany


cancers <- c(cancers1,cancers2)


#SWARM -t 2 -g 5 --time 00:30:00
args = commandArgs(trailingOnly=TRUE)
st <- args[1]
sampleName <- args[2]

#st="SKCM_2022_Sudmeier"
#sampleName="pt16"

sigType <- "Secreted"

sts <- unique(meta[meta[,"Method"]=="Visium",1])
#for(st in sts)
#{
	QCPath.st <- paste0(QCPath,st,"/")
	dir.create(QCPath.st)
	QCFilterPath.st <- paste0(QCFilterPath,st,"/")
	dir.create(QCFilterPath.st)
	
	sampleNames <- meta[meta[,"Study"]==st&meta[,"Method"]=="Visium"&meta[,"Organism"]=="Human","Sample_Name"]
	#for(sampleName in sampleNames)
	#{	
		QCPath.st.sample <- paste0(QCPath.st,sampleName,"/")
		dir.create(QCPath.st.sample)
		QCFilterPath.st.sample <- paste0(QCFilterPath.st,sampleName,"/")
		dir.create(QCFilterPath.st.sample)
		
		for(version in c("vst","weighted","weighted2","vst_condition_logUMI_cellType","vst_condition_cellType"))
		{
			combStat <- list()
			combStat[["old"]] <- data.frame()
			combStat[["new_0.05"]] <- data.frame()
			combStat[["new_0.1"]] <- data.frame()
			combStat[["new_0.15"]] <- data.frame()
			combStat[["new_0.2"]] <- data.frame()
			combStat[["new_0.25"]] <- data.frame()
			
			comb <- data.frame()
			combTest <- data.frame()
			combFDR <- data.frame()
			combComp <- data.frame()
			combComp_r <- data.frame()
			for(cancer in cancers)
			{
				if(grepl("-",cancer))
				{
					#if(!file.exists(paste0(QCPath.st.sample,"ICGC_filter_",cancer,".csv"))) next
					smyMat <- read.csv(paste0(QCPath.st.sample,"ICGC_filter_",cancer,"_",version,".csv"),row.names=1)
					cancer <- paste0("ICGC_",cancer)
				}else{
					#if(!file.exists(paste0(QCPath.st.sample,"TCGA_Xena_filter_log_",cancer,".csv"))) next
					smyMat <- read.csv(paste0(QCPath.st.sample,"TCGA_Xena_filter_log_",cancer,"_",version,".csv"),row.names=1)
					cancer <- paste0("TCGA_",cancer)
				}
				
				smyMat <- smyMat[!is.na(smyMat[,1]),,drop=F]
				
				comb[rownames(smyMat),cancer] <- smyMat[,1]
				
				smyMat <- cbind(smyMat,Group="NA")
				smyMat[rownames(smyMat)%in%IPs,"Group"] <- "Intracellular"
				smyMat[rownames(smyMat)%in%MPs,"Group"] <- "Membrane"
				smyMat[rownames(smyMat)%in%SPs,"Group"] <- "Secreted"
				smyMat <- smyMat[smyMat[,"Group"]!="NA",]
				
				x1 <- smyMat[smyMat[,"Group"]==sigType,1]
				x2 <- smyMat[smyMat[,"Group"]=="Intracellular",1]
				pv <- signif(wilcox.test(x1, x2, alternative="greater")$p.value,2)
				#medianDif <- median(x1)-median(x2)
				
				combTest[cancer,"Var2"] <- cancer
				combTest[cancer,"bar"] <- quantile(x2,0.95)
				#combTest[cancer,"Significant"] <- ifelse(pv<0.01&medianDif>0,"Yes","No")
				combTest[cancer,"Significant"] <- ifelse(pv<0.01,"Yes","No")
				
				#if(combTest[cancer,"Significant"]=="Yes")
				#{
					names(x1) <- rownames(smyMat)[smyMat[,"Group"]==sigType]
					filterSPs <- names(x1)[x1>combTest[cancer,"bar"]]
					combStat[["old"]][filterSPs,cancer] <- 1
					
					combComp[cancer,"old"] <- length(filterSPs)
					combComp_r[cancer,"old"] <- NA
				#}
				
				
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
				
				
				if(TRUE)
				{
					if(sampleName=="pt16"&cancer=="TCGA_SKCM")
					{
						
						fg.df <- data.frame()
						
						x2_sorted <- sort(x2,decreasing=T)
						tempMin <- 1
						for(i in length(x2_sorted):1 )
						{
							X <- sum( x2 >= x2_sorted[[i]] )/length(x2)
							Y <- sum( x1 >= x2_sorted[[i]] )/length(x1)
							
							fg.df[i,"r"] <- x2_sorted[[i]]
							fg.df[i,"fraction"] <- (length(x2_sorted)-i+1)/length(x2_sorted)
							fg.df[i,"Intracellular"] <- X
							fg.df[i,"Secreted"] <- Y
							fg.df[i,"FDR"] <- min(X/Y,tempMin)
							
							tempMin <- min(X/Y,tempMin)
						}
						
						fg.df[,"FDR"] <- 1 - fg.df[,"FDR"]
						
						library(ggplot2)
						p2 <- ggplot(fg.df, aes(x = r))+
							geom_point(aes(y = Secreted), color="#ff8080", size=0.1)+
							geom_point(aes(y = Intracellular), color="skyblue", size=0.1)+
							geom_point(aes(y = FDR), color="orange", size=0.1)+
							#scale_colour_manual(values = c("orange","skyblue","#ff8080"))+
							geom_vline(xintercept = combComp_r["TCGA_SKCM","new_0.25"], colour="black", linewidth=0.6, linetype="dashed")+
							xlab("Correlation r")+
							ylab("Fraction")+
							scale_y_continuous(
  								  # Features of the first axis
  								  name = "Fraction",
  								  # Add a second axis and specify its features
  								  sec.axis = sec_axis(~1-., name="False Discovery Rate")
  								  #trans = "reverse"
  								)+
							theme_bw()+ 
							theme(
							  panel.grid = element_blank(),
							  panel.background = element_blank(),
							  panel.border = element_rect(colour = "black"),
							  axis.text = element_text(colour = "black"),
							  axis.title = element_text(colour = "black")
							) 
						ggsave("/data/rub2/project/Secretome/results/fig1_QC/QC2_pt16_SKCM_SPs_vs_IPs_FDR_new.png", p2, width = 7.0, height = 5.2, dpi=200, units = "cm",limitsize = FALSE)
						writeLines(as.character(combComp_r["TCGA_SKCM","new_0.25"]),"/data/rub2/project/Secretome/results/fig1_QC/QC2_pt16_SKCM_SPs_vs_IPs_FDR_new.txt")
					}
				}
				
				
			}
			
			write.csv(combTest,paste0(QCFilterPath.st.sample,"TCGA_Xena_filter_log_",sigType,"_filterBar_",version,"_Test.csv"), quote=F)
			write.csv(combFDR,paste0(QCFilterPath.st.sample,"TCGA_Xena_filter_log_",sigType,"_filterBar_",version,"_FDR.csv"), quote=F)
			write.csv(combComp,paste0(QCFilterPath.st.sample,"TCGA_Xena_filter_log_",sigType,"_filterBar_",version,"_comp.csv"), quote=F)
			
			for(i in names(combStat) )
			{
				tempMat <- combStat[[i]]
				tempMat[is.na(tempMat)] <- 0
				write.csv(tempMat,paste0(QCFilterPath.st.sample,"TCGA_Xena_filter_log_",sigType,"_filterBar_",version,"_",i,".csv"), quote=F)
			}
		
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
				
			ggsave(paste0(QCFilterPath.st.sample,"TCGA_Xena_filter_log_",sigType,"_filterBar_",version,".png"), p2, width = 50, height = 22, dpi=200, units = "cm",limitsize = FALSE)
			
			
			
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
				) + facet_wrap(~ Var2, ncol = 10 ) #scales = "free"
				
			ggsave(paste0(QCFilterPath.st.sample,"TCGA_Xena_filter_log_",sigType,"_filterBar_",version,"_FDR.png"), p2, width = 50, height = 22, dpi=200, units = "cm",limitsize = FALSE)
		
		}# version
		
#	}
#}




