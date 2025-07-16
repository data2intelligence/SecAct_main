source("Secretome_s0_path.R")
source("../../SpaCE/code/ifun.R")
source("importHPA.R")

#SWARM -t 2 -g 80 --time 5:30:00
args = commandArgs(trailingOnly=TRUE)
st <- args[1]
sampleName <- args[2]

sts <- unique(meta[,1])
#for(st in sts)
#{
	signaturePath.st <- paste0(signaturePath,st,"/")
	dir.create(signaturePath.st)
	
	QCPath.st <- paste0(QCPath,st,"/")
	dir.create(QCPath.st)
	
	deconvResPath.st <- paste0(deconvResPath,st,"/")
	dir.create(deconvResPath.st)
	
	sampleNames <- meta[meta[,"Study"]==st,"Sample_Name"]
	#for(sampleName in sampleNames)
	#{	
		signaturePath.st.sample <- paste0(signaturePath.st,sampleName,"/")
		dir.create(signaturePath.st.sample)
		
		QCPath.st.sample <- paste0(QCPath.st,sampleName,"/")
		dir.create(QCPath.st.sample)
		
		deconvResPath.st.sample <- paste0(deconvResPath.st,sampleName,"/")
		dir.create(deconvResPath.st.sample)
		
		
		#for(version in c("vst","weighted","weighted2","vst_condition_logUMI_cellType","vst_condition_cellType"))
		for(version in c("vst_condition_logUMI_cellType","vst_condition_cellType"))
		{
			if(version%in%c("vst","weighted","weighted2"))
			{
				x <- scan(paste0(signaturePath.st.sample,"spatial_sig_vst_inhouse_s0_r3.tsv"))
				dims <- floor(sqrt(length(x) * 2))
				m <- matrix(NA, dims, dims)
				m[upper.tri(m, diag = TRUE)] <- x
				m <- t(m)
				m[upper.tri(m)] <- t(m)[upper.tri(m)]
				
				st.matrix.data.vst <- read.table(paste0(signaturePath.st.sample,"vst.tsv"), sep="\t", check.names=F)
				rownames(m) <- rownames(st.matrix.data.vst)
				colnames(m) <- rownames(st.matrix.data.vst)
				
				if(version=="vst")
				{
					sigmat <- m
				}else{
					weight_vec <- read.table(paste0(deconvResPath.st.sample,version,".tsv"), sep="\t")
					weight_vec <- weight_vec[,1]
					
					sigmat <- sweep(m, 1, weight_vec, "*")
				}
				
				rm(x);rm(m);gc()
				
			}else{
				st.matrix.data.vst <- as.matrix(read.table(paste0(signaturePath.st.sample,"vst_condition_cellType.tsv"), sep="\t",check.names=F))
		
				W <- calWeights(colnames(st.matrix.data.vst), r=3, diag0=TRUE)
				st.matrix.data.vst <- st.matrix.data.vst[,colnames(W)]
				
				m <- spatialCrossCorrelation(st.matrix.data.vst, W)
				
				sigmat <- m
				
				rm(m);gc()
			}




		
# TCGA
cancers <- list.files(paste0("/data/rub2/data/TCGA_firebrowse/mRNAseq_gene_exp_symbol"))

#genes <- c()
#for(cancer in cancers)
#{
#	cdata <- read.tpm(cancer)
#	genes <- c(genes, rownames(cdata))
#}
#genes <- unique(genes)
#SPs[!SPs%in%genes]

for(cancer in cancers)
{
	#for(ppl in c("Xena","CIDE"))
	for(ppl in c("Xena"))
	{
		if(ppl=="Xena")
		{
			cdata <- read.Xena(cancer)	
		}else{
			cdata <- read.CIDE(cancer)	
		}
		
		#for(geneScale in c("genomeWide","filter"))
		for(geneScale in c("filter"))
		{
			if(geneScale=="filter")
			{
				cdata <- filter.counts(cdata)
			}
			
			#for(trans in c("log","loglog"))
			for(trans in c("log"))
			{
				if(trans=="loglog")
				{
					cdata <- log2(cdata+1)
				}
				
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
				
				write.csv(smyMat,paste0(QCPath.st.sample,"TCGA_",ppl,"_",geneScale,"_",trans,"_",cancer,"_",version,".csv"),quote=F)
				
			}
		}
	}
}



# ICGC
dataPath <- "/data/rub2/data/ICGC/"
allFiles <- list.files(dataPath)
allFiles <- allFiles[grepl("seq.expression",allFiles)]
cancers <- gsub(".seq.expression.gz","",allFiles)

for(cancer in cancers)
{
	cdata_T <- as.matrix(read.csv(paste0(dataPath,cancer,".seq.expression.gz"),sep="\t",row.names=1))	
	
	#print(paste0(cancer," ",ncol(cdata_T)," ",nrow(cdata_T)," ",nrow(cdata_T_alt))) }
	#print(paste0(cancer," ",ncol(cdata_T) )) }
	
	if(ncol(cdata_T)<40) next
	
	#for(geneScale in c("genomeWide","filter"))
	for(geneScale in c("filter"))
	{
		if(geneScale=="filter")
		{
			cdata_T <- filter.counts(cdata_T)
		}
	
		rownames(cdata_T) <- transferSymbol(rownames(cdata_T))
		cdata_T <- rm_duplicates(cdata_T)
		
		cdata_T_minusBG <- cdata_T-rowMeans(cdata_T)
		cdata_T_minusBG <- as.matrix(cdata_T_minusBG)
		
		smyMat <- cor.act.exp(X=sigmat, Y=cdata_T_minusBG)
		
		write.csv(smyMat,paste0(QCPath.st.sample,"ICGC_",geneScale,"_",cancer,"_",version,".csv"),quote=F)			
		
		
	}
}




}







## CPTAC
#dataPath <- "/data/rub2/data/CPTAC/"
#cancers <- c("BRCA", "COAD", "GBM", "HNSC", "KIRC", "LUAD", "LUSC", "OV", "PDAC", "UCEC")
#
#geneAnno <- read.csv(paste0(dataPath,"Proteome_BCM_GENCODE_v34_harmonized_v1/README/Gene_annotation_and_representable_isoform_mapping_table.txt"),sep="\t",check.names=F)
#geneAnno[geneAnno[,1]=="ENSG00000284024.2",4] <- "MSANTD7"
#
##SPs[!SPs%in%geneAnno[,4]]
## [1] "AC233755.1"
#
#geneAnno <- geneAnno[,c(1,4)]
#geneAnno <- geneAnno[!duplicated(geneAnno[,1]),]
#rownames(geneAnno) <- geneAnno[,1]
#
#
## CPTAC RNA
#for(cancer in cancers)
#{
#	cdata_T <- as.matrix(read.csv(paste0(dataPath,"RNA_BCM_v1/",cancer,"_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt"),sep="\t",row.names=1,check.names=F))
#	
#	if(file.exists(paste0(dataPath,"RNA_BCM_v1/",cancer,"_RNAseq_gene_RSEM_coding_UQ_1500_log2_Normal.txt")))
#	{
#		cdata_N <- as.matrix(read.csv(paste0(dataPath,"RNA_BCM_v1/",cancer,"_RNAseq_gene_RSEM_coding_UQ_1500_log2_Normal.txt"),sep="\t",row.names=1,check.names=F))
#		
#		olp <- intersect(rownames(cdata_T),rownames(cdata_N))
#		cdata_T_minusBG <- cdata_T[olp,]-rowMeans(cdata_N[olp,])
#		
#		## print(paste0(cancer, " ",ncol(cdata_N))) }}
#		## [1] "HNSC 46"
#		## [1] "KIRC 72"
#		## [1] "LUAD 101"
#		## [1] "LUSC 94"
#		## [1] "PDAC 15"
#		
#	}else{
#		cdata_T_minusBG <- cdata_T-rowMeans(cdata_T)
#	}
#	
#	rownames(cdata_T_minusBG) <- geneAnno[rownames(cdata_T_minusBG),2]
#	rownames(cdata_T_minusBG) <- transferSymbol(rownames(cdata_T_minusBG))
#	cdata_T_minusBG <- rm_duplicates(cdata_T_minusBG)
#
#
#	X <- sigmat
#	Y <- as.matrix(cdata_T_minusBG)
#	
#	
#	olp <- intersect(rownames(X),rownames(Y))
#	X_olp <- X[olp,,drop=F]
#	Y_olp <- Y[olp,,drop=F]
#
#	
#	cc_corr_r <- cor(X_olp,Y_olp,use="pairwise.complete.obs")	
#	
#	olp <- intersect(rownames(cc_corr_r),rownames(Y_olp))
#	
#	smyMat <- data.frame()
#	for(i in olp)
#	{
#		smyMat[i,"r"] <- cor(cc_corr_r[i,],Y_olp[i,])
#	}
#	
#	
#	write.csv(smyMat,paste0(QCPath.st.sample,"/CPTAC_RNA_genomeWide_",cancer,".csv"),quote=F)	
#	
#}
#
#
## CPTAC RNA filtered by protein
#for(cancer in cancers)
#{
#	cdata_T <- as.matrix(read.csv(paste0(dataPath,"RNA_BCM_v1/",cancer,"_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt"),sep="\t",row.names=1,check.names=F))
#	cdata_T.p <- as.matrix(read.csv(paste0(dataPath,"Proteome_BCM_GENCODE_v34_harmonized_v1/",cancer,"_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt"),sep="\t",row.names=1,check.names=F))
#	cdata_T <- cdata_T[rownames(cdata_T)%in%rownames(cdata_T.p),]
#	
#	if(file.exists(paste0(dataPath,"RNA_BCM_v1/",cancer,"_RNAseq_gene_RSEM_coding_UQ_1500_log2_Normal.txt")))
#	{
#		cdata_N <- as.matrix(read.csv(paste0(dataPath,"RNA_BCM_v1/",cancer,"_RNAseq_gene_RSEM_coding_UQ_1500_log2_Normal.txt"),sep="\t",row.names=1,check.names=F))
#		
#		olp <- intersect(rownames(cdata_T),rownames(cdata_N))
#		cdata_T_minusBG <- cdata_T[olp,]-rowMeans(cdata_N[olp,])
#		
#		## print(paste0(cancer, " ",ncol(cdata_N))) }}
#		## [1] "HNSC 46"
#		## [1] "KIRC 72"
#		## [1] "LUAD 101"
#		## [1] "LUSC 94"
#		## [1] "PDAC 15"
#		
#	}else{
#		cdata_T_minusBG <- cdata_T-rowMeans(cdata_T)
#	}
#	
#	rownames(cdata_T_minusBG) <- geneAnno[rownames(cdata_T_minusBG),2]
#	rownames(cdata_T_minusBG) <- transferSymbol(rownames(cdata_T_minusBG))
#	cdata_T_minusBG <- rm_duplicates(cdata_T_minusBG)
#
#
#	X <- sigmat
#	Y <- as.matrix(cdata_T_minusBG)
#	
#	
#	olp <- intersect(rownames(X),rownames(Y))
#	X_olp <- X[olp,,drop=F]
#	Y_olp <- Y[olp,,drop=F]
#
#	
#	cc_corr_r <- cor(X_olp,Y_olp,use="pairwise.complete.obs")	
#	
#	olp <- intersect(rownames(cc_corr_r),rownames(Y_olp))
#	
#	smyMat <- data.frame()
#	for(i in olp)
#	{
#		smyMat[i,"r"] <- cor(cc_corr_r[i,],Y_olp[i,])
#	}
#	
#	
#	write.csv(smyMat,paste0(QCPath.st.sample,"/CPTAC_RNA_FilteredByProtein_",cancer,".csv"),quote=F)	
#	
#}
#
#
## CPTAC Protein
#for(cancer in cancers)
#{
#	cdata_T <- as.matrix(read.csv(paste0(dataPath,"Proteome_BCM_GENCODE_v34_harmonized_v1/",cancer,"_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt"),sep="\t",row.names=1,check.names=F))
#
#	if(file.exists(paste0(dataPath,"Proteome_BCM_GENCODE_v34_harmonized_v1/",cancer,"_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt")))
#	{
#		cdata_N <- as.matrix(read.csv(paste0(dataPath,"Proteome_BCM_GENCODE_v34_harmonized_v1/",cancer,"_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt"),sep="\t",row.names=1,check.names=F))
#		
#		olp <- intersect(rownames(cdata_T),rownames(cdata_N))
#		cdata_T_minusBG <- cdata_T[olp,]-rowMeans(cdata_N[olp,],na.rm=T)
#		
#		## print(paste0(cancer, " ",ncol(cdata_N))) }}
#		## [1] "COAD 100"
#		## [1] "HNSC 62"
#		## [1] "KIRC 80"
#		## [1] "LUAD 101"
#		## [1] "LUSC 99"
#		## [1] "OV 19"
#		## [1] "PDAC 44"		
#		
#	}else{
#		cdata_T_minusBG <- cdata_T-rowMeans(cdata_T,na.rm=T)
#	}
#	
#	rownames(cdata_T_minusBG) <- geneAnno[rownames(cdata_T_minusBG),2]
#	rownames(cdata_T_minusBG) <- transferSymbol(rownames(cdata_T_minusBG))
#	cdata_T_minusBG <- rm_duplicates(cdata_T_minusBG)
#	
#
#	X <- sigmat
#	Y <- as.matrix(cdata_T_minusBG)
#	
#	
#	olp <- intersect(rownames(X),rownames(Y))
#	X_olp <- X[olp,,drop=F]
#	Y_olp <- Y[olp,,drop=F]
#
#	
#	cc_corr_r <- cor(X_olp,Y_olp,use="pairwise.complete.obs")	
#	
#	olp <- intersect(rownames(cc_corr_r),rownames(Y_olp))
#	
#	smyMat <- data.frame()
#	for(i in olp)
#	{
#		smyMat[i,"r"] <- cor(cc_corr_r[i,],Y_olp[i,],use="pairwise.complete.obs")
#	}
#	
#	
#	write.csv(smyMat,paste0(QCPath.st.sample,"/CPTAC_Protein_",cancer,".csv"),quote=F)	
#	
#}
#
#
#
## GTEx
#dataPath <- "/data/rub2/data/GTEx/"
#cancers <- gsub(".self_subtract.gz","",list.files(dataPath))
#
#for(cancer in cancers)
#{	
#	cdata_T_minusBG <- as.matrix(read.csv(gzfile(paste0(dataPath,cancer,".self_subtract.gz")),row.names=1,sep="\t"))
#	
#	rownames(cdata_T_minusBG) <- transferSymbol(rownames(cdata_T_minusBG))
#	cdata_T_minusBG <- rm_duplicates(cdata_T_minusBG)
#	
#	X <- sigmat
#	Y <- as.matrix(cdata_T_minusBG)
#	
#	olp <- intersect(rownames(X),rownames(Y))
#	X_olp <- X[olp,olp,drop=F]
#	Y_olp <- Y[olp,,drop=F]
#
#	cc_corr_r <- cor(X_olp,Y_olp,use="pairwise.complete.obs")	
#	
#	olp <- intersect(rownames(cc_corr_r),rownames(Y_olp))
#	
#	smyMat <- data.frame()
#	for(i in olp)
#	{
#		smyMat[i,"r"] <- cor(cc_corr_r[i,],Y_olp[i,])
#	}
#	
#	write.csv(smyMat,paste0(QCPath.st.sample,"/GTEx_",cancer,".csv"),quote=F)
#}



#dataPath <- "/data/rub2/project/Secretome/data/"
#cancers <- c("GSE84976","GSE22138")
#
#for(cancer in cancers)
#{
#	cdata_T_minusBG <- as.matrix(read.csv(gzfile(paste0(dataPath,cancer,".self_subtract.trans.gz")),sep="\t"))
#	rownames(cdata_T_minusBG) <- transferSymbol(rownames(cdata_T_minusBG))
#	cdata_T_minusBG <- rm_duplicates(cdata_T_minusBG)
#		
#	
#	X <- sigmat
#	Y <- as.matrix(cdata_T_minusBG)
#	
#	
#	olp <- intersect(rownames(X),rownames(Y))
#	X_olp <- X[olp,,drop=F]
#	Y_olp <- Y[olp,,drop=F]
#
#	cc_corr <- psych::corr.test(X_olp,Y_olp,method="pearson",adjust="none",ci=FALSE)
#	cc_corr_r <- cc_corr$r
#	
#	
#	olp <- intersect(rownames(cc_corr_r),rownames(Y_olp))
#	corFun <- function(x)
#	{
#		corRes <- cor.test(cc_corr_r[x,],Y_olp[x,])
#		c(corRes$estimate,corRes$p.value)
#	}
#	
#	smy <- parallel::mclapply(olp, corFun, mc.cores=coreNo)
#	smyMat <- matrix(unlist(smy),ncol=2,byrow=T)
#	smyMat <- cbind(smyMat, padj=p.adjust(smyMat[,2],method="BH"))
#	rownames(smyMat) <- olp
#	colnames(smyMat) <- c("r","pv","padj")
#	
#	write.csv(smyMat,paste0(QCPath.st.sample,"/GEO_",cancer,".csv"),quote=F)			
#}		


	#}
#}

