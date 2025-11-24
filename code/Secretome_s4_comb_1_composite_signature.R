source("Secretome_s0_path.R")

#SWARM -t 30 -g 100 --time 02:30:00
args = commandArgs(trailingOnly=TRUE)
cancer <- args[1]
version <- args[2]


compositeSigList <- list()

# all signatures without filtering
# signatures filtered by auto correlation only
compositeSigList[["AllSig"]] <- data.frame()
compositeSigList[["AllSigFilteredBy_MoranI"]] <- data.frame()

all_stat <- data.frame()
significant_moranI_stat <- data.frame()
significant_moranI_value <- data.frame()

sts <- unique(meta[,"Study"])
for(st in sts)
{
	QCPath.st <- paste0(QCPath,st,"/")
	
	sampleNames <- meta[meta[,"Study"]==st,"Sample_Name"]
	for(sampleName in sampleNames)
	{	
		QCPath.st.sample <- paste0(QCPath.st,sampleName,"/")
		
		moranIMat <- read.csv(paste0(QCPath.st.sample, "AutoCorrelation_",version,".csv"),row.names=1)
		
		all_stat[rownames(moranIMat),paste0(st,"@",sampleName)] <- 1
		
		moranIMat <- moranIMat[moranIMat[,"p.Moran_Padj"]<0.01 & moranIMat[,"p.Moran_I"]>0,]
		
		significant_moranI_stat[rownames(moranIMat),paste0(st,"@",sampleName)] <- 1
		significant_moranI_value[rownames(moranIMat),paste0(st,"@",sampleName)] <- moranIMat[,"p.Moran_I"]
	}
}

all_stat <- as.matrix(all_stat)
all_stat[is.na(all_stat)] <- 0

significant_moranI_stat <- as.matrix(significant_moranI_stat)
significant_moranI_stat[is.na(significant_moranI_stat)] <- 0

compositeSigList[["AllSig"]] <- all_stat
compositeSigList[["AllSigFilteredBy_MoranI"]] <- significant_moranI_stat

significant_moranI_value <- as.matrix(significant_moranI_value)


# signatures filtered by gating only
# old strategy
compositeSigList[["AllSigFilteredBy_bulk_TCGA"]] <- data.frame()
compositeSigList[["AllSigFilteredBy_bulk_ICGC"]] <- data.frame()
compositeSigList[["AllSigFilteredBy_bulk_TCGA_ICGC"]] <- data.frame()

for(st in sts)
{
	QCFilterPath.st <- paste0(QCFilterPath,st,"/")
	
	sampleNames <- meta[meta[,"Study"]==st,"Sample_Name"]
	for(sampleName in sampleNames)
	{	
		QCFilterPath.st.sample <- paste0(QCFilterPath.st,sampleName,"/")
	
		if(file.info(paste0(QCFilterPath.st.sample,"Secreted_filterBar_",version,"_old.csv"))$size==1) next
		filterBar <- read.csv(paste0(QCFilterPath.st.sample,"Secreted_filterBar_",version,"_old.csv"),row.names=1,check.names=F)
		
		compositeSigList[["AllSigFilteredBy_bulk_TCGA"]] [rownames(filterBar),paste0(st,"@",sampleName)] <- rowSums(filterBar[,grepl("TCGA",colnames(filterBar)),drop=F])
		compositeSigList[["AllSigFilteredBy_bulk_ICGC"]] [rownames(filterBar),paste0(st,"@",sampleName)] <- rowSums(filterBar[,grepl("ICGC",colnames(filterBar)),drop=F])
		compositeSigList[["AllSigFilteredBy_bulk_TCGA_ICGC"]] [rownames(filterBar),paste0(st,"@",sampleName)] <- rowSums(filterBar)
	}
}

for(i in 3:5)
{
	compositeSigList[[i]] <- as.matrix(compositeSigList[[i]])
	compositeSigList[[i]][is.na(compositeSigList[[i]])] <- 0
}



# signatures filtered by both auto-correlation and gating
# old gating strategy
compositeSigList[["AllSigFilteredBy_MoranI_TCGA"]] <- data.frame()
compositeSigList[["AllSigFilteredBy_MoranI_ICGC"]] <- data.frame()
compositeSigList[["AllSigFilteredBy_MoranI_TCGA_ICGC"]] <- data.frame()

for(i in 6:8)
{
	mat1 <- compositeSigList[["AllSigFilteredBy_MoranI"]]
	mat2 <- compositeSigList[[i-3]]
	
	olp1 <- intersect(rownames(mat1),rownames(mat2))
	olp2 <- intersect(colnames(mat1),colnames(mat2))
	
	mat1_olp <- mat1[olp1,olp2]
	mat2_olp <- mat2[olp1,olp2]
	
	mat3 <- mat1_olp*mat2_olp
	
	compositeSigList[[i]] <- mat3
}



# signatures filtered by gating only
# new strategy
compositeSigList[["AllSigFilteredBy_bulk_TCGA_ICGC_0.05"]] <- data.frame()
compositeSigList[["AllSigFilteredBy_bulk_TCGA_ICGC_0.1"]] <- data.frame()
compositeSigList[["AllSigFilteredBy_bulk_TCGA_ICGC_0.15"]] <- data.frame()
compositeSigList[["AllSigFilteredBy_bulk_TCGA_ICGC_0.2"]] <- data.frame()
compositeSigList[["AllSigFilteredBy_bulk_TCGA_ICGC_0.25"]] <- data.frame()

for(st in sts)
{
	QCFilterPath.st <- paste0(QCFilterPath,st,"/")
	
	sampleNames <- meta[meta[,"Study"]==st,"Sample_Name"]
	for(sampleName in sampleNames)
	{	
		QCFilterPath.st.sample <- paste0(QCFilterPath.st,sampleName,"/")
		
		for(ratio in c(0.05 * 1:5))
		{
			if(file.info(paste0(QCFilterPath.st.sample,"Secreted_filterBar_",version,"_new_",ratio,".csv"))$size==1) next
			filterBar <- read.csv(paste0(QCFilterPath.st.sample,"Secreted_filterBar_",version,"_new_",ratio,".csv"),row.names=1,check.names=F)
			
			compositeSigList[[paste0("AllSigFilteredBy_bulk_TCGA_ICGC_",ratio)]] [rownames(filterBar),paste0(st,"@",sampleName)] <- rowSums(filterBar)
		}
	}
}

for(i in 9:13)
{
	compositeSigList[[i]] <- as.matrix(compositeSigList[[i]])
	compositeSigList[[i]][is.na(compositeSigList[[i]])] <- 0
}



# signatures filtered by both auto-correlation and gating
# old gating strategy
compositeSigList[["AllSigFilteredBy_MoranI_TCGA_ICGC_0.05"]] <- data.frame()
compositeSigList[["AllSigFilteredBy_MoranI_TCGA_ICGC_0.1"]] <- data.frame()
compositeSigList[["AllSigFilteredBy_MoranI_TCGA_ICGC_0.15"]] <- data.frame()
compositeSigList[["AllSigFilteredBy_MoranI_TCGA_ICGC_0.2"]] <- data.frame()
compositeSigList[["AllSigFilteredBy_MoranI_TCGA_ICGC_0.25"]] <- data.frame()

for(i in 14:18)
{
	mat1 <- compositeSigList[["AllSigFilteredBy_MoranI"]]
	mat2 <- compositeSigList[[i-5]]
	
	olp1 <- intersect(rownames(mat1),rownames(mat2))
	olp2 <- intersect(colnames(mat1),colnames(mat2))
	
	mat1_olp <- mat1[olp1,olp2]
	mat2_olp <- mat2[olp1,olp2]
	
	mat3 <- mat1_olp*mat2_olp
	
	compositeSigList[[i]] <- mat3
}



# 19 cancer type-specific
# 20 leave one out signature
if(cancer!="Pancancer")
{
	cancer_vec <- sapply(strsplit(colnames(compositeSigList[[18]]),"_",fixed=T),function(x) return(x[1]))

	compositeSigList[[paste0("AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_",cancer,"_filterByPan")]] <- 
		compositeSigList[[18]][,cancer_vec%in%cancer,drop=F]
			
	compositeSigList[[paste0("AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_non-",cancer,"_filterByPan")]] <- 
		compositeSigList[[18]][,!cancer_vec%in%cancer]
}



for(i in 1:length(compositeSigList))
{
	compositeSigList[[i]] <- compositeSigList[[i]] [rowSums(compositeSigList[[i]]) > 0,]
	print(dim(compositeSigList[[i]]))
	
	#dir.create(paste0(signatureCombPath,names(compositeSigList)[i]))
}






n <- ifelse(cancer=="Pancancer",18,19)

for(i in n:length(compositeSigList))
{
	dir.create(paste0(signatureCombPath,names(compositeSigList)[i]))
}

sigMean <- function(sigmat)
{
	sigmat <- sigmat[rowSums(is.na(sigmat))<ncol(sigmat)*0.2,,drop=F]
	apply(sigmat,1,function(x) mean(x,na.rm=T))
}

sigMean.weighted <- function(sigmat,weightST)
{
	sigmat <- sigmat[rowSums(is.na(sigmat))<ncol(sigmat)*0.2,,drop=F]
	sigmat <- t(t(sigmat)*weightST)
	apply(sigmat,1,function(x) sum(x,na.rm=T))
}

composite <- function(gene)
{
	if(file.info(paste0(signatureCombPath,"combSig_all_in_one_SP/",gene,"_",version,".rds"))$size>100)
	{
		sigmat <- readRDS(paste0(signatureCombPath,"combSig_all_in_one_SP/",gene,"_",version,".rds"))
		
		if(ncol(sigmat)>0) 
		{
			for(compSig in names(compositeSigList)[n:length(compositeSigList)])
			{
				compositeSig <- compositeSigList[[compSig]]
					
				if(gene%in%rownames(compositeSig))
				{
					sigSamples <- colnames(compositeSig)[compositeSig[gene,] > 0]
					sigmat_sub <- sigmat[,colnames(sigmat)%in%sigSamples,drop=F]
					
					if(ncol(sigmat_sub)>0) 
					{
						weightST <- compositeSig[gene,colnames(sigmat_sub)]
						weightST <- weightST/sum(weightST)
						
						sigmat_sub_mean <- sigMean.weighted(sigmat_sub,weightST)
						write.csv(sigmat_sub_mean,paste0(signatureCombPath,compSig,"/",gene,"_",version,".csv"),quote=F)
						write.csv(ncol(sigmat_sub),paste0(signatureCombPath,compSig,"/",gene,"_",version,"_stat.csv"),quote=F)
					}
				}
				
			}
		}
	}
}

parallel::mclapply(SPs, composite, mc.cores=20) 



for(compSig in names(compositeSigList)[n:length(compositeSigList)])
{
	genes <- list.files(paste0(signatureCombPath,compSig))
	genes <- genes[grepl(paste0(version,".csv"),genes)]
	genes <- gsub(paste0("_",version,".csv"),"",genes)
	genes <- unique(genes)
	
	genes <- genes[ !grepl("IGH",genes) ]
	genes <- genes[ !grepl("IGK",genes) ]
	genes <- genes[ !grepl("IGL",genes) ]
	
	genes <- genes[ !grepl("MT-",genes,fixed=T) ]
	
	Receptors <- readLines(paste0(HPAPath,"CIDE_Receptor_genes"))
	Receptors <- transferSymbol(Receptors)
	genes <- genes[!genes%in%Receptors]

	
	smry <- data.frame()
	stat <- data.frame()
	for(gene in genes)
	{
		sigmat_sub_mean <- read.csv(paste0(signatureCombPath,compSig,"/",gene,"_",version,".csv"),header=T,row.names=1,check.names=F)
		smry[rownames(sigmat_sub_mean),gene] <- sigmat_sub_mean[,1]
		
		sigmat_sub_stat <- read.csv(paste0(signatureCombPath,compSig,"/",gene,"_",version,"_stat.csv"),header=T,row.names=1,check.names=F)
		stat[gene,"sigCount"] <- sigmat_sub_stat[1,1]
		
		x <- compositeSigList[[compSig]][gene,]
		x_names <- names(x)[x>0]
		x_names_datasets <- sapply(strsplit(x_names,"@",fixed=T),function(x) return(x[1]))
		stat[gene,"dsCount"] <- length(unique(x_names_datasets))
	}
	
	write.csv(stat,paste0(signatureCombPath,compSig,"_",version,"_stat.csv"),quote=F)
		


	for(m in 3)
	{
		smry_m <- as.matrix(smry[, stat[,"dsCount"] >= m ])
		
		smry_m <- smry_m[apply(smry_m,1,function(x) sum(is.na(x)))==0,]
		print(dim(smry_m))
		
		write.table(smry_m,paste0(signatureCombPath,compSig,"_ds",m,"_",version,".tsv"),quote=F,sep="\t")
		
		if(cancer=="Pancancer"&version=="vst")
		{
			write.table(smry_m, gzfile(paste0(finalSignaturesPath,"SecAct.tsv.gz")), quote=F, sep="\t")
		}
		
	}
	
	
}	

