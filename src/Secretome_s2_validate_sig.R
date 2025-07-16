source("Secretome_s0_path.R")
source("../../SpaCE/code/ifun.R")
source("../../SpaCE/code/SpaCE.R")
source("importHPA.R")

version <- "vst_condition_logUMI_cellType"

# validate TGFB1 and single target
if(FALSE)
{
	# TGFB1 to targets or partners
	genes <- c("TGFB1","TGFB2","TGFB3")
	genes <- c("TGFB1")
	pts  <- c("LTBP1","LTBP2","LTBP3","LTBP4","LRRC32","ITGAV","ITGB6","ITGB8","TGFBR1","TGFBR2","TGFBR3","SMAD2","SMAD3","SMAD4","CCN2","SERPINE1","TGFBI")
	pts  <- c("LTBP2","ITGAV","TGFBR2","SERPINE1")
	pts  <- c("LTBP1","LTBP2","LTBP3","LTBP4","LRRC32","ITGAV","ITGB6","ITGB8","TGFBR1","TGFBR2","TGFBR3","CCN2","SERPINE1","TGFBI")
	pts  <- c("SERPINE1","CCN2","TGFBI","LTBP2","TGFBR2","ITGAV")
	targets <- c("SERPINE1","CCN2","TGFBI")
	mediators <- c("LTBP2","TGFBR2","ITGAV")
	
	for(gene in genes)
	{
		comb <- read.table(paste0(signatureCombPath,"combSig_mean_SP_s0_r3/",gene,"_",version,".tsv"), sep="\t", check.names=F)
		comb <- comb[rowSums(is.na(comb))<ncol(comb)*0.3,,drop=F]
		
		fg.df <- data.frame()
		pvs <- c()
		for(pt in pts)
		{
			temp <- unlist(comb[pt,])
			temp <- temp[!is.na(temp)]
			if(length(temp)==0) next
			
			wilcox_res <- wilcox.test(temp, mu=0)
			pv <- signif(wilcox_res$p.value,2)
			pvs <- c(pvs, pv)
			
			fg.df[paste0(gene,pt,names(temp)),"x"] <- pt
			fg.df[paste0(gene,pt,names(temp)),"y"] <- temp
			fg.df[paste0(gene,pt,names(temp)),"z"] <- gene
		}
		names(pvs) <- pts
		
		fg.df[,1] <- factor(fg.df[,1],levels= c(targets, setdiff(names(sort(pvs)),targets) ) )
		fg.df <- cbind(fg.df,group="mediator")
		fg.df[fg.df[,1]%in%targets,"group"] <- "target"
		
		library(ggplot2)
		p1 <- ggplot(fg.df,aes(x=x,y=y,colour=group)) + 
			geom_hline(yintercept=0, color = "grey", linewidth=0.8, linetype="dashed")+
			geom_jitter(alpha=0.3, size=0.2, width=0.25,)+
			geom_boxplot(color="black",fill="white", alpha=0.5, width=0.5, outlier.shape = NA)+
			scale_colour_manual(values=c("#41b9c1","#8080FF"))+
			annotate("text", x = 0.75, y=0.2, label = paste0("p = "), size=2.8)+
			annotate("text", x = 1, y=0.16, label = paste0("",pvs[1]), size=2.8)+
			annotate("text", x = 2, y=0.16, label = paste0("",pvs[2]), size=2.8)+
			annotate("text", x = 3, y=0.16, label = paste0("",pvs[3]), size=2.8)+
			annotate("text", x = 4, y=0.16, label = paste0("",pvs[4]), size=2.8)+
			annotate("text", x = 5, y=0.16, label = paste0("",pvs[5]), size=2.8)+
			annotate("text", x = 6, y=0.16, label = paste0("",pvs[6]), size=2.8)+
			coord_cartesian(ylim = c(-0.1, 0.2))+
			ylab("Moran's I")+
			theme_classic()+ 
			theme(
				axis.text.x = element_text(angle = 38, hjust = 1),
				axis.title.x = element_blank(),
				legend.position="none"
			)
		ggsave(paste0("/data/rub2/project/Secretome/results/fig1_QC/moranI_example_single_",gene,".png"), p1, width = 10.8, height = 5.3, dpi=300, units = "cm",limitsize = FALSE)
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	# TGFB1-targets vs TGFB1-random
	gene <- "TGFB1"
	comb <- read.table(paste0(signatureCombPath,"combSig_mean_SP_s0_r3/",gene,"_",version,".tsv"), sep="\t", check.names=F)
	comb <- comb[rowSums(is.na(comb))<ncol(comb)*0.3,,drop=F]
		
	comb_vec <- apply(comb,1,function(x) median(x,na.rm=T))
	
	comb_vec["TGFBI"]
	comb_vec["TGFBR1"]
	comb_vec["SERPINE1"]
	
	#gmt <- read.gmt("/data/rub2/project/SpaCE/data/c2.all.v2023.2.Hs.symbols.gmt")
	TGFB1_targets <- c("ACVR1","APC","ARID4B","BCAR3","BMP2","BMPR1A","BMPR2","CDH1","CDK9","CDKN1C","CTNNB1","ENG","FKBP1A","FNTA","FURIN","HDAC1","HIPK2","ID1","ID2","ID3","IFNGR2","JUNB","KLF10","LEFTY2","LTBP2","MAP3K7","NCOR2","NOG","PMEPA1","PPM1A","PPP1CA","PPP1R15A","RAB31","RHOA","SERPINE1","SKI","SKIL","SLC20A1","SMAD1","SMAD3","SMAD6","SMAD7","SMURF1","SMURF2","SPTBN1","TGFB1","TGFBR1","TGIF1","THBS1","TJP1","TRIM33","UBE2D3","WWTR1","XIAP")
	#TGFB1_targets <- gmt[["KEGG_TGF_BETA_SIGNALING_PATHWAY"]] #KEGG_TGF_BETA_SIGNALING_PATHWAY REACTOME_SIGNALING_BY_TGFB_FAMILY_MEMBERS
	TGFB1_targets <- unique(transferSymbol(TGFB1_targets))
	TGFB1_targets <- intersect(names(comb_vec),TGFB1_targets)
	
	
	set.seed(123)
	TGFB1_nonTargets <- sample(setdiff(names(comb_vec),TGFB1_targets),length(TGFB1_targets))
	
	mean(comb_vec[TGFB1_targets])
	mean(comb_vec[TGFB1_nonTargets])
	wilcox_res <- wilcox.test(comb_vec[TGFB1_targets], comb_vec[TGFB1_nonTargets])
	pv <- signif(wilcox_res$p.value,3)
	
	fg.df <- data.frame(x=c(rep("Targets",length(TGFB1_targets)),rep("Random",length(TGFB1_nonTargets))),y=c(comb_vec[TGFB1_targets], comb_vec[TGFB1_nonTargets]))
	
	library(ggplot2)
	p2 <- ggplot(fg.df,aes(x=x, y=y, fill=x)) + 
		geom_hline(yintercept=0, color = "grey", linewidth=0.8, linetype="dashed")+
		geom_jitter(alpha=0.2,size=1)+
		geom_boxplot(colour="grey1", alpha=0.6, outlier.shape = NA)+
		scale_fill_manual(values=c("darkgrey","#8080FF"))+
		annotate("text", x = 1.5, y=0.028, label = paste0("p = ",pv))+
		ylab("Median Moran's I")+
		xlab("")+
		theme_bw()+ 
		theme(
			panel.background = element_blank(),
			panel.grid = element_blank(),
			axis.title = element_text(colour = "black"),
			axis.text = element_text(colour = "black"),
			legend.position="none"
		)
	
	library(patchwork)
	ggsave("/data/rub2/project/Secretome/results/fig1_QC/moranI_example_TGFB1.png", p2, width = 11, height = 11, dpi=300, units = "cm",limitsize = FALSE)
	
	



# visualize co-expression

temp_sorted <- sort(temp,decreasing=T)

#SWARM -t 2 -g 30 --time 01:00:00
args = commandArgs(trailingOnly=TRUE)
st <- args[1]
sampleName <- args[2]

sts <- unique(meta[meta[,"Method"]=="Visium",1])
#for(st in sts)
#{	
	sampleNames <- meta[meta[,"Study"]==st,"Sample_Name"]

	rawDataPath.st <- paste0(rawDataPath,st,"/")
	preprocessDataPath.st <- paste0(preprocessDataPath,st,"/")
	QCPath.st <- paste0(QCPath,st,"/")
	deconvResFigPath.st <- paste0(deconvResFigPath,st,"/")
	dir.create(deconvResFigPath.st)
	
	#for(sampleName in sampleNames)
	#{
		if(!paste0(st,"@",sampleName)%in%names(temp_sorted)) next
		
		i <- which(names(temp_sorted)==paste0(st,"@",sampleName))
		
		deconvResFigPath.st.sample <- paste0(deconvResFigPath.st,sampleName,"/")
		dir.create(deconvResFigPath.st.sample)
		QCPath.st.sample <- paste0(QCPath.st,sampleName,"/")
		dir.create(QCPath.st.sample)
		
		st.matrix.data <- as.matrix(read.table(paste0(preprocessDataPath.st,sampleName,"_counts.tsv"),check.names=FALSE))
		rownames(st.matrix.data) <- transferSymbol(rownames(st.matrix.data))
		st.matrix.data <- rm_duplicates(st.matrix.data)
		
		#st.matrix.data.scaled <- t(t(st.matrix.data)*1e5/colSums(st.matrix.data))
		#st.matrix.data.log <- log2(st.matrix.data.scaled + 1 )
		#res_deconv <- st.matrix.data.log[c("TGFB1","SERPINE1"),,drop=F]


		for(j in c("SERPINE1","LTBP2"))
		{
			res_deconv <- st.matrix.data[c("TGFB1",j),,drop=F]
			res_deconv[1,res_deconv[1,]>0] <- 1
			res_deconv[2,res_deconv[2,]>0] <- 2
			res_deconv <- rbind(res_deconv, comb=as.character(colSums(res_deconv)))
			res_deconv <- res_deconv[3,,drop=F]
			
			library(cowplot)
			icol <- 1
			irow <- nrow(res_deconv)/icol
	
			for(iround in 1:irow)
			{
				cellTypes <- rownames(res_deconv)[(icol*iround-icol+1):(icol*iround)]
			
				for(cellType in cellTypes)
				{
					Content <- res_deconv[cellType,]
						
					g <- visiualSpatial(
								visiualVector=Content,
								ggtitle=cellType,
								xlab="",
								ylab="",
								legendName="Expr", 
								legendExist="right", 
								scaleType="color-discrete",
								colorDiscreteOption="our",
								colors=c("grey","#FF8080","#8080FF","brown"),
								dotAlpha=1,
								HEanno=FALSE,
								bgImage=FALSE,
								label=FALSE
					)
				
					if(!exists("gg")){
						gg <- plot_grid(g + theme(legend.position="none"), nrow=1) 
					}else{
						gg <- plot_grid(gg, g + theme(legend.position="none"), nrow=1, rel_widths=c(which(cellTypes==cellType)-1, 1)) 
					}
				} # cellType
				
				if(!exists("ggg")){
					ggg <- plot_grid(gg, ncol=1) 
				}else{
					ggg <- plot_grid(ggg, gg, ncol=1, rel_heights=c(iround-1, 1)) 
				}
				
				rm(gg)
			}
			
			legend.obj <- get_legend(g)
	
			ggg <- plot_grid(
				ggg, legend.obj, nrow = 1, 
				rel_widths = c(1,0.05)
				)
			
			ggsave(paste0("/data/rub2/project/Secretome/results/fig1_QC/TGFB1/",i,"_",j,"_",st,"@",sampleName,"_comb.png"), ggg, width = icol*7, height = irow*7, dpi=300, units = "cm")
	
			rm(ggg)
		
		}
	#}
#}


	
}








# validate TGFB1 and signature score
if(TRUE)
{
	
st="SKCM_2022_Sudmeier"
sampleName="pt16"

#SWARM -t 2 -g 30 --time 00:20:00
args = commandArgs(trailingOnly=TRUE)
cancer <- args[1]

sts <- unique(meta[meta[,"Method"]=="Visium",1])
#for(st in sts)
#{
	signaturePath.st <- paste0(signaturePath,st,"/")
	dir.create(signaturePath.st)
	
	QCFilterPath.st <- paste0(QCFilterPath,st,"/")
	dir.create(QCFilterPath.st)
	
	deconvResPath.st <- paste0(deconvResPath,st,"/")
	dir.create(deconvResPath.st)
	
	sampleNames <- meta[meta[,"Study"]==st&meta[,"Method"]=="Visium"&meta[,"Organism"]=="Human","Sample_Name"]
	#for(sampleName in sampleNames)
	#{	
		signaturePath.st.sample <- paste0(signaturePath.st,sampleName,"/")
		dir.create(signaturePath.st.sample)
		
		QCFilterPath.st.sample <- paste0(QCFilterPath.st,sampleName,"/")
		dir.create(QCFilterPath.st.sample)
		
		deconvResPath.st.sample <- paste0(deconvResPath.st,sampleName,"/")
		dir.create(deconvResPath.st.sample)
	
		#x <- scan(paste0(signaturePath.st.sample,"spatial_sig_vst_inhouse_s0_r3.tsv"))
		#dims <- floor(sqrt(length(x) * 2))
		#m <- matrix(NA, dims, dims)
		#m[upper.tri(m, diag = TRUE)] <- x
		#m <- t(m)
		#m[upper.tri(m)] <- t(m)[upper.tri(m)]
		#
		#st.matrix.data.vst <- read.table(paste0(signaturePath.st.sample,"vst.tsv"), sep="\t", check.names=F)
		#rownames(m) <- rownames(st.matrix.data.vst)
		#colnames(m) <- rownames(st.matrix.data.vst)
		#
		#weight_vec <- read.table(paste0(deconvResPath.st.sample,"weighted.tsv"), sep="\t")
		#weight_vec <- weight_vec[,1]
		#
		#sigmat <- sweep(m, 1, weight_vec, "*")
		
		
		st.matrix.data.vst <- as.matrix(read.table(paste0(signaturePath.st.sample, version, ".tsv"), sep="\t",check.names=F))
		
		W <- calWeights(colnames(st.matrix.data.vst), r=3, diag0=TRUE)
		st.matrix.data.vst <- st.matrix.data.vst[,colnames(W)]
		
		m <- spatialCrossCorrelation(st.matrix.data.vst, W)
		
		sigmat <- m
		
		
		

		dataset <- strsplit(cancer,"_")[[1]][1]
		cancer <- strsplit(cancer,"_")[[1]][3]
		if(dataset=="TCGA")
		{
			cdata <- read.Xena(cancer)	
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
			dataPath <- "/data/rub2/data/ICGC/"
			cdata_T <- as.matrix(read.csv(paste0(dataPath,cancer,".seq.expression.gz"),sep="\t",row.names=1))	
			cdata_T <- filter.counts(cdata_T)
			
			rownames(cdata_T) <- transferSymbol(rownames(cdata_T))
			cdata_T <- rm_duplicates(cdata_T)
		
			cdata_T_minusBG <- cdata_T-rowMeans(cdata_T)
			cdata_T_minusBG <- as.matrix(cdata_T_minusBG)
		}
		
		
		X = sigmat
		Y = cdata_T_minusBG
		olp <- intersect(rownames(X),rownames(Y))
		X_olp <- X[olp,,drop=F]
		Y_olp <- Y[olp,,drop=F]
		
		
		# signature score
		cc_corr_r <- WGCNA::cor(X_olp,Y_olp,use="pairwise.complete.obs")	
		
		
		
		if(cancer=="SKCM")
		{
			for(i in c("TGFB1"))
			{
				fg.df <- data.frame(x=Y_olp[i,],y=cc_corr_r[i,])
				cor_res <- cor.test(fg.df[,1],fg.df[,2])
				rv <- round(cor_res$estimate,2)
				pv <- signif(cor_res$p.value,2)
				
				library(ggplot2)
				p1 <- ggplot(fg.df,aes(x=x, y=y)) + 
					geom_point(color="#ff8080", alpha=0.5, size=0.6)+
					annotate("text", x=-2.3, y=0.48, label=paste0("r = ",rv,"\np = ",pv) ,size=3)+
					xlab("TGFB1 Expression")+
					ylab("Signature Score")+
					theme_classic()+ 
					theme(
						panel.background = element_rect(fill='transparent'),
						plot.background = element_rect(fill='transparent', color=NA),
						axis.title = element_text(colour = "black"),
						axis.text = element_text(colour = "black")
					)#+geom_smooth(method = 'lm',color="brown")		
				ggsave(paste0("/data/rub2/project/Secretome/results/fig1_QC/QC0_scatter_",i,"_expr_sigScore.png"), p1, width = 5.6, height = 5.2, dpi=200, units = "cm")
			}
		}
		
		
		
	
		# correlation between signature score and expression

		# real
		smyMat <- cor.act.exp(X=sigmat, Y=cdata_T_minusBG)
		write.csv(smyMat,paste0("/data/rub2/project/Secretome/results/fig1_QC/Real_Random/QC0_pt16_",cancer,"_Real.csv"),quote=F)
		
		
		# random
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
		write.csv(combPermute,paste0("/data/rub2/project/Secretome/results/fig1_QC/Real_Random/QC0_pt16_",cancer,"_Random.csv"),quote=F)

#	} #sampleName
#} #st


}
	

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
if(FALSE)
{	
	# real vs random
	comb <- data.frame()
	for(group in c("Random","Real"))
	{
		smyMat <- read.csv(paste0("/data/rub2/project/Secretome/results/fig1_QC/Real_Random/QC0_pt16_SKCM_",group,".csv"),row.names=1)
		
		smyMat <- smyMat[!is.na(smyMat[,1]),,drop=F]
		
		comb[rownames(smyMat),group] <- smyMat[,1]
	}
	comb <- comb[rownames(comb)%in%SPs,]
	
	dim(comb)
	## [1] 928   2
	
	comb.m <- reshape2::melt(as.matrix(comb))
	comb.m <- comb.m[!is.na(comb.m[,3]),]
	
	comb.m[comb.m[,1]=="TGFB1",]
	## 	      Var1   Var2       value
	## 808  TGFB1 Random -0.01575136
	## 1736 TGFB1   Real  0.59244769

	x1 <- comb.m[comb.m[,"Var2"]=="Random",3]
	x2 <- comb.m[comb.m[,"Var2"]=="Real",3]
	pv <- signif(wilcox.test(x1, x2)$p.value,2)
		
	library(ggplot2)
	p2 <- ggplot(comb.m,aes(x = Var2, y = value, colour=Var2))+
		geom_hline(yintercept=0, color = "grey", linewidth=0.6, linetype="dashed")+
		geom_jitter(alpha=0.3, size=0.2, width=0.25,)+
		geom_boxplot(color="black",fill="white", alpha=0.5, width=0.5, outlier.shape = NA)+
		scale_color_manual(values = c("grey","#ff8080"))+
		annotate("text", x=1.2, y=0.85, label=paste0("p = ",pv) ,size=3)+
		xlab("Signature")+
		ylab("Correlation r")+
		theme_classic()+ 
		theme(
		  panel.grid = element_blank(),
		  panel.background = element_blank(),
		  axis.text = element_text(colour = "black"),
		  axis.title = element_text(colour = "black"),
		  legend.position="none"
		) 
		
	ggsave("/data/rub2/project/Secretome/results/fig1_QC/QC0_pt16_SKCM_SPs.png", p2, width = 5, height = 5, dpi=200, units = "cm",limitsize = FALSE)

	
	
	
	dataPath <- "/data/rub2/data/TCGA_firebrowse/mRNAseq_gene_exp_symbol"
	cancers <- list.files(dataPath)
	cancers <- setdiff(cancers, "LAML")
	cancers1 <- cancers
	
	dataPath <- "/data/rub2/data/ICGC/"
	allFiles <- list.files(dataPath)
	allFiles <- allFiles[grepl("seq.expression",allFiles)]
	cancers <- gsub(".seq.expression.gz","",allFiles)
	cancers <- setdiff(cancers,"CLLE-ES") 
	cancers <- setdiff(cancers,c("BPLL-FR","PAEN-AU","PRAD-FR")) # sample <40
	cancers2 <- cancers
	#CLLE-ES	Chronic Lymphocytic Leukemia - ES	Spain
	#MALY-DE	Malignant Lymphoma - DE	Germany
	
	cancers <- c(cancers1,cancers2)

	
	comb <- data.frame()
	for(cancer in cancers)
	{
		for(group in c("Random","Real"))
		{
			smyMat <- read.csv(paste0("/data/rub2/project/Secretome/results/fig1_QC/Real_Random/QC0_pt16_",cancer,"_",group,".csv"),row.names=1)
			
			smyMat <- smyMat[!is.na(smyMat[,1]),,drop=F]
			
			
			comb[paste0(cancer,group,rownames(smyMat)),"cancer"] <- ifelse(grepl("-",cancer),paste0("ICGC_",cancer),paste0("TCGA_",cancer))
			comb[paste0(cancer,group,rownames(smyMat)),"Group"] <- group
			comb[paste0(cancer,group,rownames(smyMat)),"gene"] <- rownames(smyMat)
			comb[paste0(cancer,group,rownames(smyMat)),"value"] <- smyMat[,1]

		}
	}
	comb.m <- comb[comb[,"gene"]%in%SPs,]
	
	
	comb.m[,"cancer"] <- factor(
		comb.m[,"cancer"],
		levels=sapply(cancers, function(x) ifelse(grepl("-",x),paste0("ICGC_",x),paste0("TCGA_",x)))
	)
	
	library(ggplot2)
	library(ggsignif)
	p2 <- ggplot(comb.m,aes(x = Group, y = value, colour=Group))+
		geom_hline(yintercept=0, color = "grey", linewidth=0.6, linetype="dashed")+
		geom_jitter(alpha=0.3, size=0.2, width=0.25,)+
		geom_boxplot(color="black",fill="white", alpha=0.5, width=0.5, outlier.shape = NA)+
		geom_signif( comparisons = list(c("Random","Real")) ,y_position = c(0.72), test = "wilcox.test", color="black")+
		scale_color_manual(values = c("grey","#ff8080"))+
		xlab("Signature")+
		ylab("Pearson r value")+
		guides(colour = guide_legend(override.aes = list(size=5)))+
		theme_bw()+ 
		theme(
		  panel.grid = element_blank(),
		  panel.background = element_blank(),
		  axis.text = element_text(size=15,colour = "black"),
		  axis.text.x = element_blank(),
		  axis.ticks.x = element_blank(),
		  axis.title.x = element_blank(),
		  axis.title = element_text(size=15,colour = "black"),
		  legend.position="bottom",
		  legend.text=element_text(size=15),
		  legend.title=element_text(size=15),
		  strip.text = element_text(size = 13)
		)+facet_wrap(~ cancer, ncol =11) #scales = "free"
	
	ggsave("/data/rub2/project/Secretome/results/fig1_QC/QC0_pt16_allCancerType_SPs.png", p2, width = 55, height = 23, dpi=200, units = "cm", limitsize = FALSE)

	
	


	
	# Secreted vs Intracellular
	
	smyMat <- read.csv(paste0("/data/rub2/project/Secretome/results/fig1_QC/Real_Random/QC0_pt16_SKCM_Real.csv"),row.names=1)
	smyMat <- smyMat[!is.na(smyMat[,1]),,drop=F]
	
	smyMat <- cbind(smyMat,Group="NA")
	smyMat[rownames(smyMat)%in%IPs,"Group"] <- "Intracellular"
	smyMat[rownames(smyMat)%in%MPs,"Group"] <- "Membrane"
	smyMat[rownames(smyMat)%in%SPs,"Group"] <- "Secreted"
	smyMat <- smyMat[smyMat[,"Group"]!="NA",]
	
	smyMat <- smyMat[smyMat[,"Group"]%in%c("Secreted","Intracellular"),]
	
	x1 <- smyMat[smyMat[,"Group"]=="Secreted",1]
	x2 <- smyMat[smyMat[,"Group"]=="Intracellular",1]
	pv <- signif(wilcox.test(x1, x2, alternative="greater")$p.value,2)
	
	
	library(ggplot2)
	library(ggsignif)
	p2 <- ggplot(smyMat,aes(x = Group, y = x))+
		#geom_hline(yintercept=0, color = "grey", linewidth=0.6, linetype="dashed")+
		geom_violin(aes(group=Group, fill=Group),trim=FALSE)+
		scale_fill_manual(values = c("skyblue","#ff8080"))+
		geom_boxplot(width=0.15,outlier.shape=NA)+
		geom_hline(yintercept = 0.45, colour="black", linewidth=0.6, linetype="dashed")+
		annotate("text", x=1.5, y=1, label=paste0("p = ",pv) ,size=3)+
		xlab("Group")+
		ylab("Correlation r")+
		theme_classic()+ 
		theme(
		  panel.grid = element_blank(),
		  panel.background = element_blank(),
		  axis.text = element_text(colour = "black"),
		  axis.title = element_text(colour = "black"),
		  axis.title.x = element_blank(),
		  legend.position="none"
		) 
	ggsave("/data/rub2/project/Secretome/results/fig1_QC/QC2_pt16_SKCM_SPs_vs_IPs.png", p2, width = 6, height = 5.2, dpi=200, units = "cm",limitsize = FALSE)

	
	
	
	
}







if(FALSE)
{

	comb <- c()
	gene <- "TGFB1"
	targ <- "IFNG"
	
	sts <- unique(meta[meta[,"Method"]=="Visium",1])
	for(st in sts)
	{
		signaturePath.st <- paste0(signaturePath,st,"/")
		 
		sampleNames <- meta[meta[,"Study"]==st&meta[,"Method"]=="Visium"&meta[,"Organism"]=="Human","Sample_Name"]
		for(sampleName in sampleNames)
		{	
			signaturePath.st.sample <- paste0(signaturePath.st,sampleName,"/")
			signaturePath.st.sample.singleSig <- paste0(signaturePath.st,sampleName,"/singleSig/")
			
			fileName <- paste0(signaturePath.st.sample.singleSig,gene,".tsv")
			if(file.exists(fileName))
			{
				m_sub <- read.table(fileName, sep="\t")
				if(targ%in%rownames(m_sub))
				{
					comb <- c(comb, m_sub[targ,1])
				}
			}
		}
	}

stat<-data.frame(comb)
library(ggplot2)
p0 <- ggplot(stat,aes(x=comb)) + 
	geom_histogram( colour="grey1", fill="gainsboro", alpha=0.5)+
	geom_density(color="grey3",linewidth=0.6)+
	ggtitle(paste0("MoranI_",gene,"_",targ, ": n=",length(comb),",mean=",mean(comb)))+
	ylab("")+
	xlab("")+
	theme_bw()+ 
	theme(
		panel.background = element_blank(),
		panel.grid = element_blank(),
		axis.title = element_text(size=14,colour = "black"),
		axis.text = element_text(size=13,colour = "black"),
		legend.position="none"
	)
ggsave(paste0("check_MoranI_",gene,"_",targ,".png"), p0, width = 12, height = 9, dpi=500, units = "cm", limitsize =FALSE)




	comb <- c()
	comb_names <- c()
	gene <- "TGFB1"
	targ <- "IFNG"
	
	sts <- unique(meta[meta[,"Method"]=="Visium",1])
	for(st in sts)
	{
		signaturePath.st <- paste0(signaturePath,st,"/")
		 
		sampleNames <- meta[meta[,"Study"]==st&meta[,"Method"]=="Visium"&meta[,"Organism"]=="Human","Sample_Name"]
		for(sampleName in sampleNames)
		{	
			signaturePath.st.sample <- paste0(signaturePath.st,sampleName,"/")
			signaturePath.st.sample.singleSig <- paste0(signaturePath.st,sampleName,"/singleSig/")
			
			fileName1 <- paste0(signaturePath.st.sample.singleSig,gene,".tsv")
			fileName2 <- paste0(signaturePath.st.sample.singleSig,targ,".tsv")
			if(file.exists(fileName1)&file.exists(fileName2))
			{
				m_sub1 <- read.table(fileName1, sep="\t")
				m_sub2 <- read.table(fileName2, sep="\t")
				olp <- intersect(rownames(m_sub1),rownames(m_sub2))
				comb <- c(comb,cor(m_sub1[olp,1],m_sub2[olp,1],use="pairwise.complete.obs"))
				comb_names <- c(comb_names,paste0(st,"@",sampleName))
			}
			
		}
	}

names(comb) <- comb_names
stat<-data.frame(comb)

write.csv(stat, paste0("check_corr_",gene,"_",targ,".csv"), quote=F)


library(ggplot2)
p0 <- ggplot(stat,aes(x=comb)) + 
	geom_histogram( colour="grey1", fill="gainsboro", alpha=0.5)+
	geom_density(color="grey3",linewidth=0.6)+
	ggtitle(paste0("sig_",gene,"_",targ, ": n=",length(comb),",mean=",mean(comb)))+
	ylab("")+
	xlab("")+
	theme_bw()+ 
	theme(
		panel.background = element_blank(),
		panel.grid = element_blank(),
		axis.title = element_text(size=14,colour = "black"),
		axis.text = element_text(size=13,colour = "black"),
		legend.position="none"
	)
ggsave(paste0("check_corr_",gene,"_",targ,".png"), p0, width = 12, height = 9, dpi=500, units = "cm", limitsize =FALSE)





	sample_neg <- names(comb)[comb<0]
	
	comb1 <- data.frame()
	comb2 <- data.frame()
	gene <- "TGFB1"
	targ <- "IFNG"
	
	sts <- unique(meta[meta[,"Method"]=="Visium",1])
	for(st in sts)
	{
		signaturePath.st <- paste0(signaturePath,st,"/")
		 
		sampleNames <- meta[meta[,"Study"]==st&meta[,"Method"]=="Visium"&meta[,"Organism"]=="Human","Sample_Name"]
		for(sampleName in sampleNames)
		{	
			signaturePath.st.sample <- paste0(signaturePath.st,sampleName,"/")
			signaturePath.st.sample.singleSig <- paste0(signaturePath.st,sampleName,"/singleSig/")
			
			fileName <- paste0(signaturePath.st.sample.singleSig,gene,".tsv")
			if(file.exists(fileName))
			{
				m_sub1 <- read.table(fileName, sep="\t")
				comb1[rownames(m_sub1),paste0(st,"@",sampleName)] <- m_sub1[,1]
			}
			fileName <- paste0(signaturePath.st.sample.singleSig,targ,".tsv")
			if(file.exists(fileName))
			{
				m_sub2 <- read.table(fileName, sep="\t")
				comb2[rownames(m_sub2),paste0(st,"@",sampleName)] <- m_sub2[,1]
			}
			
		}
	}


comb1 <- comb1[,sample_neg]
comb2 <- comb2[,sample_neg]

comb1 <- comb1[rowSums(is.na(comb1))<ncol(comb1)*0.2,,drop=F]
comb2 <- comb2[rowSums(is.na(comb2))<ncol(comb2)*0.2,,drop=F]

comb1.mean <- apply(comb1,1,function(x) mean(x,na.rm=T))
comb2.mean <- apply(comb2,1,function(x) mean(x,na.rm=T))

olp <- intersect(names(comb1.mean),names(comb2.mean))

length(olp)
12864

cor(comb1.mean[olp],comb2.mean[olp])
-0.2136392








	comb1 <- data.frame()
	comb2 <- data.frame()
	gene <- "TGFB1"
	targ <- "IFNG"
	
	sts <- unique(meta[meta[,"Method"]=="Visium",1])
	for(st in sts)
	{
		signaturePath.st <- paste0(signaturePath,st,"/")
		 
		sampleNames <- meta[meta[,"Study"]==st&meta[,"Method"]=="Visium"&meta[,"Organism"]=="Human","Sample_Name"]
		for(sampleName in sampleNames)
		{	
			signaturePath.st.sample <- paste0(signaturePath.st,sampleName,"/")
			signaturePath.st.sample.singleSig <- paste0(signaturePath.st,sampleName,"/singleSig/")
			
			fileName <- paste0(signaturePath.st.sample.singleSig,gene,".tsv")
			if(file.exists(fileName))
			{
				m_sub1 <- read.table(fileName, sep="\t")
				comb1[rownames(m_sub1),paste0(st,"@",sampleName)] <- m_sub1[,1]
			}
			fileName <- paste0(signaturePath.st.sample.singleSig,targ,".tsv")
			if(file.exists(fileName))
			{
				m_sub2 <- read.table(fileName, sep="\t")
				comb2[rownames(m_sub2),paste0(st,"@",sampleName)] <- m_sub2[,1]
			}
			
		}
	}


comb1 <- comb1[,c("CRC_2022_Qi@Patient6","CRC_2022_Wu@ST-liver2","BRCA_10x_Datasets@Version2.0.0_Breast.Cancer")]
comb2 <- comb2[,c("CRC_2022_Qi@Patient6","CRC_2022_Wu@ST-liver2","BRCA_10x_Datasets@Version2.0.0_Breast.Cancer")]

comb1 <- comb1[rowSums(is.na(comb1))<ncol(comb1)*0.2,,drop=F]
comb2 <- comb2[rowSums(is.na(comb2))<ncol(comb2)*0.2,,drop=F]

comb1.mean <- apply(comb1,1,function(x) mean(x,na.rm=T))
comb2.mean <- apply(comb2,1,function(x) mean(x,na.rm=T))

olp <- intersect(names(comb1.mean),names(comb2.mean))

length(olp)
14744

cor(comb1.mean[olp],comb2.mean[olp])
0.8176188






	comb1 <- data.frame()
	comb2 <- data.frame()
	gene <- "TGFB1"
	targ <- "IFNG"
	
	sts <- unique(meta[meta[,"Method"]=="Visium",1])
	for(st in sts)
	{
		signaturePath.st <- paste0(signaturePath,st,"/")
		 
		sampleNames <- meta[meta[,"Study"]==st&meta[,"Method"]=="Visium"&meta[,"Organism"]=="Human","Sample_Name"]
		for(sampleName in sampleNames)
		{	
			signaturePath.st.sample <- paste0(signaturePath.st,sampleName,"/")
			signaturePath.st.sample.singleSig <- paste0(signaturePath.st,sampleName,"/singleSig/")
			
			fileName <- paste0(signaturePath.st.sample.singleSig,gene,".tsv")
			if(file.exists(fileName))
			{
				m_sub1 <- read.table(fileName, sep="\t")
				comb1[rownames(m_sub1),paste0(st,"@",sampleName)] <- m_sub1[,1]
			}
			fileName <- paste0(signaturePath.st.sample.singleSig,targ,".tsv")
			if(file.exists(fileName))
			{
				m_sub2 <- read.table(fileName, sep="\t")
				comb2[rownames(m_sub2),paste0(st,"@",sampleName)] <- m_sub2[,1]
			}
			
		}
	}

comb1 <- comb1[rowSums(is.na(comb1))<ncol(comb1)*0.2,,drop=F]
comb2 <- comb2[rowSums(is.na(comb2))<ncol(comb2)*0.2,,drop=F]

comb1.mean <- apply(comb1,1,function(x) mean(x,na.rm=T))
comb2.mean <- apply(comb2,1,function(x) mean(x,na.rm=T))

olp <- intersect(names(comb1.mean),names(comb2.mean))

length(olp)
12954

cor(comb1.mean[olp],comb2.mean[olp],use="pairwise.complete.obs")
0.569781







	comb1 <- data.frame()
	comb2 <- data.frame()
	gene <- "TGFB1"
	targ <- "IFNG"
	
	sts <- unique(meta[meta[,"Method"]=="Visium",1])
	for(st in sts)
	{
		signaturePath.st <- paste0(signaturePath,st,"/")
		 
		sampleNames <- meta[meta[,"Study"]==st&meta[,"Method"]=="Visium"&meta[,"Organism"]=="Human","Sample_Name"]
		for(sampleName in sampleNames)
		{	
			signaturePath.st.sample <- paste0(signaturePath.st,sampleName,"/")
			signaturePath.st.sample.singleSig <- paste0(signaturePath.st,sampleName,"/singleSig/")
			
			fileName <- paste0(signaturePath.st.sample.singleSig,gene,".tsv")
			if(file.exists(fileName))
			{
				m_sub1 <- read.table(fileName, sep="\t")
				comb1[rownames(m_sub1),paste0(st,"@",sampleName)] <- m_sub1[,1]
			}
			fileName <- paste0(signaturePath.st.sample.singleSig,targ,".tsv")
			if(file.exists(fileName))
			{
				m_sub2 <- read.table(fileName, sep="\t")
				comb2[rownames(m_sub2),paste0(st,"@",sampleName)] <- m_sub2[,1]
			}
			
		}
	}

comb1.mean <- apply(comb1,1,function(x) mean(x,na.rm=T))
comb2.mean <- apply(comb2,1,function(x) mean(x,na.rm=T))

olp <- intersect(names(comb1.mean),names(comb2.mean))

length(olp)
31366

cor(comb1.mean[olp],comb2.mean[olp],use="pairwise.complete.obs")
0.4155654


comb1_r <- cor( comb1,matrix(comb1.mean,ncol=1),use="pairwise.complete.obs" )
comb2_r <- cor( comb2,matrix(comb2.mean,ncol=1),use="pairwise.complete.obs" )







smy <- read.csv(paste0("/data/rub2/project/Secretome/results/TGFB1/","smy.csv"),row.names=1)
rownames(smy) <- gsub("_TGFB1.csv","",rownames(smy))
rownames(smy) <- sapply(strsplit(rownames(smy),"$",fixed=T),function(x) return(x[2]))
	


for(ds in c("s1_r1","s0_r2","s0_r3","s0_r4"))
{
	gene <- "TGFB1"
	targ <- "IFNG"
	comb1 <- data.frame()
	comb2 <- data.frame()
	
	sts <- unique(meta[meta[,"Method"]=="Visium",1])
	for(st in sts)
	{
		signaturePath.st <- paste0(signaturePath,st,"/")
		 
		sampleNames <- meta[meta[,"Study"]==st&meta[,"Method"]=="Visium"&meta[,"Organism"]=="Human","Sample_Name"]
		for(sampleName in sampleNames)
		{	
			signaturePath.st.sample <- paste0(signaturePath.st,sampleName,"/")
			signaturePath.st.sample.singleSig <- paste0(signaturePath.st,sampleName,"/singleSig_scale_IP_",ds,"/")
			
			fileName1 <- paste0(signaturePath.st.sample.singleSig,gene,".tsv")
			fileName2 <- paste0(signaturePath.st.sample.singleSig,targ,".tsv")
			if(file.exists(fileName1)&file.exists(fileName2))
			{
				m_sub1 <- read.table(fileName1, sep="\t")
				m_sub2 <- read.table(fileName2, sep="\t")
				
				m_sub1 <- m_sub1[!is.na(m_sub1[,1]),,drop=F]
				m_sub2 <- m_sub2[!is.na(m_sub2[,1]),,drop=F]
				
				comb1[rownames(m_sub1),paste0(st,"@",sampleName)] <- m_sub1[,1]
				comb2[rownames(m_sub2),paste0(st,"@",sampleName)] <- m_sub2[,1]
				#comb1[rownames(m_sub1),paste0(st,"@",sampleName)] <- scale(m_sub1)[,1]
				#comb2[rownames(m_sub2),paste0(st,"@",sampleName)] <- scale(m_sub2)[,1]
			}
			
		}
	}

comb1.mean <- apply(comb1,1,function(x) mean(x,na.rm=T))
comb2.mean <- apply(comb2,1,function(x) mean(x,na.rm=T))



olp <- intersect(names(comb1.mean),names(comb2.mean))

length(olp)

fg.df <- data.frame(
	TGFB1_comb_sig=comb1.mean[olp],
	IFNG_comb_sig=comb2.mean[olp]
	)
cor(comb1.mean[olp],comb2.mean[olp],use="pairwise.complete.obs")

	cor_res <- cor.test(fg.df[,1],fg.df[,2],use="pairwise.complete.obs")
	rvalue <- round(cor_res$estimate,3)
	pvalue <- signif(cor_res$p.value,3)
	
	library(ggplot2)
	s2 <- ggplot(fg.df,aes(x=TGFB1_comb_sig, y=IFNG_comb_sig)) + 
			geom_point(size=0.8)+
			ggtitle(paste0("scale row IP ",ds,", r:",rvalue," p:",pvalue))+
			theme_bw()+ 
			theme(
				plot.background = element_blank(),
				panel.grid = element_blank(),
				plot.title = element_text(hjust = 0.5,size=20),
				axis.text = element_text(size=20,colour = "black"),
				axis.title = element_text(size=20,colour = "black")
			)
	ggsave(paste0("/data/rub2/project/Secretome/results/TGFB1/","TGFB1_IFNG_comb_sig_corr_",ds,"_scale_row_IP.png"), s2, width = 15, height = 15, dpi=300, units = "cm")
}	




comb1_r <- cor( comb1,matrix(comb1.mean,ncol=1),use="pairwise.complete.obs" )
comb2_r <- cor( comb2,matrix(comb2.mean,ncol=1),use="pairwise.complete.obs" )


fg.df <- data.frame(
r_TGFB1_IFNG_expr_panel3=smy[rownames(comb1_r),1],
r_TGFB1_IFNG_sig_panel4=smy[rownames(comb1_r),2],
r_TGFB1_meanSig_singleSig=comb1_r,
r_IFNG_meanSig_singleSig=comb2_r
)
	
	cor_res <- cor.test(fg.df[,2],fg.df[,3],use="pairwise.complete.obs")
	rvalue <- round(cor_res$estimate,3)
	pvalue <- signif(cor_res$p.value,3)
	
	library(ggplot2)
	s2 <- ggplot(fg.df,aes(x=r_TGFB1_IFNG_sig_panel4, y=r_TGFB1_meanSig_singleSig)) + 
			geom_point(size=0.8)+
			ggtitle(paste0("r:",rvalue," p:",pvalue))+
			theme_bw()+ 
			theme(
				plot.background = element_blank(),
				panel.grid = element_blank(),
				plot.title = element_text(hjust = 0.5,size=20),
				axis.text = element_text(size=20,colour = "black"),
				axis.title = element_text(size=20,colour = "black")
			)
	ggsave(paste0("/data/rub2/project/Secretome/results/TGFB1/","panel4_TGFB1_scale.png"), s2, width = 15, height = 15, dpi=300, units = "cm")
	
	cor_res <- cor.test(fg.df[,2],fg.df[,4],use="pairwise.complete.obs")
	rvalue <- round(cor_res$estimate,3)
	pvalue <- signif(cor_res$p.value,3)
	
	library(ggplot2)
	s2 <- ggplot(fg.df,aes(x=r_TGFB1_IFNG_sig_panel4, y=r_IFNG_meanSig_singleSig)) + 
			geom_point(size=0.8)+
			ggtitle(paste0("r:",rvalue," p:",pvalue))+
			theme_bw()+ 
			theme(
				plot.background = element_blank(),
				panel.grid = element_blank(),
				plot.title = element_text(hjust = 0.5,size=20),
				axis.text = element_text(size=20,colour = "black"),
				axis.title = element_text(size=20,colour = "black")
			)
	ggsave(paste0("/data/rub2/project/Secretome/results/TGFB1/","panel4_IFNG_scale.png"), s2, width = 15, height = 15, dpi=300, units = "cm")
	






signature.centroid <- read.table("/data/rub2/project/Secretome/code/CytoSig-master/CytoSig/signature.centroid",sep="\t",check.names=F)
signature.centroid.12 <- read.table("/data/rub2/project/Secretome/code/CytoSig-master/CytoSig/signature.centroid.12",sep="\t",check.names=F)

colnames(signature.centroid) <- paste0("c_",colnames(signature.centroid))
colnames(signature.centroid.12) <- paste0("a_",colnames(signature.centroid.12))

olp <- intersect(rownames(signature.centroid),rownames(signature.centroid.12))
x <- cbind(signature.centroid[olp,paste0("c_",c("TGFB1","IFNG"))],signature.centroid.12[olp,paste0("a_",c("TGFB1","IFNG"))])
cor(x)

signature.centroid.12 <- read.table("/data/rub2/project/Secretome/code/CytoSig-master/CytoSig/signature.centroid.12",sep="\t",check.names=F)
cor(signature.centroid.12[,c("TGFB1","IFNG")])

signature.centroid.13 <- read.table("/data/rub2/project/Secretome/code/CytoSig-master/CytoSig/signature.centroid.13",sep="\t",check.names=F)
cor(signature.centroid.13[,c("TGFB1","IFNG")])

signature.centroid.14 <- read.table("/data/rub2/project/Secretome/code/CytoSig-master/CytoSig/signature.centroid.14",sep="\t",check.names=F)
cor(signature.centroid.14[,c("TGFB1","IFNG")])

signature.centroid.15 <- read.table("/data/rub2/project/Secretome/code/CytoSig-master/CytoSig/signature.centroid.15",sep="\t",check.names=F)
cor(signature.centroid.15[,c("TGFB1","IFNG")])


}

