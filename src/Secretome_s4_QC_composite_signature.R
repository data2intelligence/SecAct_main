source("Secretome_s0_path.R")
source("../../SpaCE/code/ifun.R")
source("importHPA.R")


#SWARM -t 50 -g 100 --time 02:30:00
args = commandArgs(trailingOnly=TRUE)
version <- args[1]


if(version%in%c("vst","weighted","weighted2"))
{
	moranIFile <- "vst_moranI.fast.csv"
}else{
	moranIFile <- paste0(version,"_moranI.fast.csv")
}


all_stat <- data.frame()
significant_moranI_stat <- data.frame()
significant_moranI_value <- data.frame()

sts <- unique(meta[meta[,"Method"]=="Visium",1])
for(st in sts)
{
	QCPath.st <- paste0(QCPath,st,"/")
	
	sampleNames <- meta[meta[,"Study"]==st&meta[,"Method"]=="Visium"&meta[,"Organism"]=="Human","Sample_Name"]
	for(sampleName in sampleNames)
	{	
		QCPath.st.sample <- paste0(QCPath.st,sampleName,"/")
		
		moranIMat <- read.csv(paste0(QCPath.st.sample, moranIFile),row.names=1)
		
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


significant_moranI_value <- as.matrix(significant_moranI_value)



#sum(all_stat)
#
#all_stat_colSums <- colSums(all_stat==1)
#sort( all_stat_colSums)
#
#all_stat_rowSums <- rowSums(all_stat==1)
#sort( all_stat_rowSums)

#sum(significant_moranI_stat)
#
#significant_moranI_stat_colSums <- colSums(significant_moranI_stat==1)
#sort( significant_moranI_stat_colSums)
#
#significant_moranI_stat_rowSums <- rowSums(significant_moranI_stat==1)
#sort( significant_moranI_stat_rowSums)



compositeSigList <- list()

compositeSigList[["AllSig"]] <- all_stat
compositeSigList[["AllSigFilteredBy_MoranI"]] <- significant_moranI_stat





# top moran I sig
#for(x in c(1:5,10,15,20,25,30))
#{
#	significant_moranI_stat2 <- significant_moranI_stat
#	significant_moranI_stat2[!is.na(significant_moranI_stat2)] <- 0
#	
#	for(gene in rownames(significant_moranI_stat2))
#	{
#		temp <- significant_moranI_value[gene,]
#		temp <- temp[!is.na(temp)]
#		temp <- sort(temp,decreasing=T)
#		
#		if(length(temp)<x)
#		{
#			significant_moranI_stat2[gene,names(temp)] <- 1
#		}else{
#			significant_moranI_stat2[gene,names(temp)[1:x]] <- 1
#		}
#	}
#	
#	compositeSigList[[paste0("AllSigFilteredBy_MoranI_Top_",x)]] <- significant_moranI_stat2
#}






compositeSigList[["AllSigFilteredBy_bulk_TCGA"]] <- data.frame()
compositeSigList[["AllSigFilteredBy_bulk_ICGC"]] <- data.frame()
compositeSigList[["AllSigFilteredBy_bulk_TCGA_ICGC"]] <- data.frame()

for(st in sts)
{
	QCFilterPath.st <- paste0(QCFilterPath,st,"/")
	
	sampleNames <- meta[meta[,"Study"]==st&meta[,"Method"]=="Visium"&meta[,"Organism"]=="Human","Sample_Name"]
	for(sampleName in sampleNames)
	{	
		QCFilterPath.st.sample <- paste0(QCFilterPath.st,sampleName,"/")
	
		if(file.info(paste0(QCFilterPath.st.sample,"TCGA_Xena_filter_log_Secreted_filterBar_",version,"_old.csv"))$size==1) next
		filterBar <- read.csv(paste0(QCFilterPath.st.sample,"TCGA_Xena_filter_log_Secreted_filterBar_",version,"_old.csv"),row.names=1,check.names=F)
		
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




# secreted

compositeSigList[["AllSigFilteredBy_bulk_TCGA_ICGC_0.05"]] <- data.frame()
compositeSigList[["AllSigFilteredBy_bulk_TCGA_ICGC_0.1"]] <- data.frame()
compositeSigList[["AllSigFilteredBy_bulk_TCGA_ICGC_0.15"]] <- data.frame()
compositeSigList[["AllSigFilteredBy_bulk_TCGA_ICGC_0.2"]] <- data.frame()
compositeSigList[["AllSigFilteredBy_bulk_TCGA_ICGC_0.25"]] <- data.frame()

for(st in sts)
{
	QCFilterPath.st <- paste0(QCFilterPath,st,"/")
	
	sampleNames <- meta[meta[,"Study"]==st&meta[,"Method"]=="Visium"&meta[,"Organism"]=="Human","Sample_Name"]
	for(sampleName in sampleNames)
	{	
		QCFilterPath.st.sample <- paste0(QCFilterPath.st,sampleName,"/")
		
		for(ratio in c(0.05 * 1:5))
		{
			if(file.info(paste0(QCFilterPath.st.sample,"TCGA_Xena_filter_log_Secreted_filterBar_",version,"_new_",ratio,".csv"))$size==1) next
			filterBar <- read.csv(paste0(QCFilterPath.st.sample,"TCGA_Xena_filter_log_Secreted_filterBar_",version,"_new_",ratio,".csv"),row.names=1,check.names=F)
			
			compositeSigList[[paste0("AllSigFilteredBy_bulk_TCGA_ICGC_",ratio)]] [rownames(filterBar),paste0(st,"@",sampleName)] <- rowSums(filterBar)
		}
	}
}


for(i in 9:13)
{
	compositeSigList[[i]] <- as.matrix(compositeSigList[[i]])
	compositeSigList[[i]][is.na(compositeSigList[[i]])] <- 0
}


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











## membrane
#
#compositeSigList[["AllSigFilteredBy_Membrane_bulk_TCGA_ICGC_0.05"]] <- data.frame()
#compositeSigList[["AllSigFilteredBy_Membrane_bulk_TCGA_ICGC_0.1"]] <- data.frame()
#compositeSigList[["AllSigFilteredBy_Membrane_bulk_TCGA_ICGC_0.15"]] <- data.frame()
#compositeSigList[["AllSigFilteredBy_Membrane_bulk_TCGA_ICGC_0.2"]] <- data.frame()
#compositeSigList[["AllSigFilteredBy_Membrane_bulk_TCGA_ICGC_0.25"]] <- data.frame()
#
#for(st in sts)
#{
#	QCFilterPath.st <- paste0(QCFilterPath,st,"/")
#	
#	sampleNames <- meta[meta[,"Study"]==st&meta[,"Method"]=="Visium"&meta[,"Organism"]=="Human","Sample_Name"]
#	for(sampleName in sampleNames)
#	{	
#		QCFilterPath.st.sample <- paste0(QCFilterPath.st,sampleName,"/")
#		
#		for(ratio in c(0.05 * 1:5))
#		{
#			if(file.info(paste0(QCFilterPath.st.sample,"TCGA_Xena_filter_log_Membrane_filterBar_new_",ratio,".csv"))$size==1) next
#			filterBar <- read.csv(paste0(QCFilterPath.st.sample,"TCGA_Xena_filter_log_Membrane_filterBar_new_",ratio,".csv"),row.names=1,check.names=F)
#			
#			compositeSigList[[paste0("AllSigFilteredBy_Membrane_bulk_TCGA_ICGC_",ratio)]] [rownames(filterBar),paste0(st,"@",sampleName)] <- rowSums(filterBar)
#		}
#	}
#}
#
#
#for(i in 9:13)
#{
#	compositeSigList[[i]] <- as.matrix(compositeSigList[[i]])
#	compositeSigList[[i]][is.na(compositeSigList[[i]])] <- 0
#}
#
#
#compositeSigList[["AllSigFilteredBy_Membrane_MoranI_TCGA_ICGC_0.05"]] <- data.frame()
#compositeSigList[["AllSigFilteredBy_Membrane_MoranI_TCGA_ICGC_0.1"]] <- data.frame()
#compositeSigList[["AllSigFilteredBy_Membrane_MoranI_TCGA_ICGC_0.15"]] <- data.frame()
#compositeSigList[["AllSigFilteredBy_Membrane_MoranI_TCGA_ICGC_0.2"]] <- data.frame()
#compositeSigList[["AllSigFilteredBy_Membrane_MoranI_TCGA_ICGC_0.25"]] <- data.frame()
#
#
#for(i in 14:18)
#{
#	mat1 <- compositeSigList[["AllSigFilteredBy_MoranI"]]
#	mat2 <- compositeSigList[[i-5]]
#	
#	olp1 <- intersect(rownames(mat1),rownames(mat2))
#	olp2 <- intersect(colnames(mat1),colnames(mat2))
#	
#	mat1_olp <- mat1[olp1,olp2]
#	mat2_olp <- mat2[olp1,olp2]
#	
#	mat3 <- mat1_olp*mat2_olp
#	
#	compositeSigList[[i]] <- mat3
#}














#compositeSigList[["AllSigFilteredBy_bulk_TCGA_0.05"]] <- data.frame()
#compositeSigList[["AllSigFilteredBy_bulk_TCGA_0.1"]] <- data.frame()
#compositeSigList[["AllSigFilteredBy_bulk_TCGA_0.15"]] <- data.frame()
#compositeSigList[["AllSigFilteredBy_bulk_TCGA_0.2"]] <- data.frame()
#compositeSigList[["AllSigFilteredBy_bulk_TCGA_0.25"]] <- data.frame()
#
#for(st in sts)
#{
#	QCFilterPath.st <- paste0(QCFilterPath,st,"/")
#	
#	sampleNames <- meta[meta[,"Study"]==st&meta[,"Method"]=="Visium"&meta[,"Organism"]=="Human","Sample_Name"]
#	for(sampleName in sampleNames)
#	{	
#		QCFilterPath.st.sample <- paste0(QCFilterPath.st,sampleName,"/")
#		
#		for(ratio in c(0.05 * 1:5))
#		{
#			if(file.info(paste0(QCFilterPath.st.sample,"TCGA_Xena_filter_log_filterBar_new_",ratio,".csv"))$size==1) next
#			filterBar <- read.csv(paste0(QCFilterPath.st.sample,"TCGA_Xena_filter_log_filterBar_new_",ratio,".csv"),row.names=1,check.names=F)
#			filterBar <- filterBar[,grepl("TCGA",colnames(filterBar)),drop=F]
#			
#			compositeSigList[[paste0("AllSigFilteredBy_bulk_TCGA_",ratio)]] [rownames(filterBar),paste0(st,"@",sampleName)] <- rowSums(filterBar)
#		}
#	}
#}
#
#
#for(i in 19:23)
#{
#	compositeSigList[[i]] <- as.matrix(compositeSigList[[i]])
#	compositeSigList[[i]][is.na(compositeSigList[[i]])] <- 0
#}
#
#
#compositeSigList[["AllSigFilteredBy_MoranI_TCGA_0.05"]] <- data.frame()
#compositeSigList[["AllSigFilteredBy_MoranI_TCGA_0.1"]] <- data.frame()
#compositeSigList[["AllSigFilteredBy_MoranI_TCGA_0.15"]] <- data.frame()
#compositeSigList[["AllSigFilteredBy_MoranI_TCGA_0.2"]] <- data.frame()
#compositeSigList[["AllSigFilteredBy_MoranI_TCGA_0.25"]] <- data.frame()
#
#
#for(i in 24:28)
#{
#	mat1 <- compositeSigList[["AllSigFilteredBy_MoranI"]]
#	mat2 <- compositeSigList[[i-5]]
#	
#	olp1 <- intersect(rownames(mat1),rownames(mat2))
#	olp2 <- intersect(colnames(mat1),colnames(mat2))
#	
#	mat1_olp <- mat1[olp1,olp2]
#	mat2_olp <- mat2[olp1,olp2]
#	
#	mat3 <- mat1_olp*mat2_olp
#	
#	compositeSigList[[i]] <- mat3
#}







## i = 19 ~ 24
#for(cancer in c("BRCA","CRC","GBM","KIRC","LIHC","OV"))
#{
#	compositeSigList[[paste0("AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_",cancer,"_filterByPan")]] <- 
#		compositeSigList[[18]][,grepl(cancer,colnames(compositeSigList[[18]]))]
#}
#
## i = 25 ~ 30
#for(cancer in c("BRCA","CRC","GBM","KIRC","LIHC","OV"))
#{
#	if(cancer=="BRCA") TCGA.ICGC <- c("TCGA_BRCA","ICGC_BRCA-KR")
#	if(cancer=="CRC") TCGA.ICGC <- c("TCGA_COAD","TCGA_READ")
#	if(cancer=="GBM")  TCGA.ICGC <- c("TCGA_GBM")
#	if(cancer=="KIRC") TCGA.ICGC <- c("TCGA_KIRC","ICGC_RECA-EU")
#	if(cancer=="LIHC") TCGA.ICGC <- c("TCGA_LIHC","ICGC_LICA-FR","ICGC_LIRI-JP")
#	if(cancer=="OV")   TCGA.ICGC <- c("TCGA_OV","ICGC_OV-AU")
#	
#	i <- paste0("AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_",cancer,"_filterBy",cancer)
#	compositeSigList[[i]] <- data.frame()
#	
#	for(st in sts[grepl(cancer,sts)])
#	{
#		QCFilterPath.st <- paste0(QCFilterPath,st,"/")
#		
#		sampleNames <- meta[meta[,"Study"]==st&meta[,"Method"]=="Visium"&meta[,"Organism"]=="Human","Sample_Name"]
#		for(sampleName in sampleNames)
#		{	
#			QCFilterPath.st.sample <- paste0(QCFilterPath.st,sampleName,"/")
#			
#			for(ratio in c(0.25))
#			{
#				if(file.info(paste0(QCFilterPath.st.sample,"TCGA_Xena_filter_log_filterBar_new_",ratio,".csv"))$size==1) next
#				filterBar <- read.csv(paste0(QCFilterPath.st.sample,"TCGA_Xena_filter_log_filterBar_new_",ratio,".csv"),row.names=1,check.names=F)
#				filterBar <- filterBar[,colnames(filterBar)%in%TCGA.ICGC,drop=F]
#				
#				compositeSigList[[i]] [rownames(filterBar),paste0(st,"@",sampleName)] <- rowSums(filterBar)
#			}
#		}
#	}
#	
#	
#	compositeSigList[[i]] <- as.matrix(compositeSigList[[i]])
#	compositeSigList[[i]][is.na(compositeSigList[[i]])] <- 0
#	
#	
#	mat1 <- compositeSigList[["AllSigFilteredBy_MoranI"]]
#	mat2 <- compositeSigList[[i]]
#	
#	olp1 <- intersect(rownames(mat1),rownames(mat2))
#	olp2 <- intersect(colnames(mat1),colnames(mat2))
#	
#	mat1_olp <- mat1[olp1,olp2]
#	mat2_olp <- mat2[olp1,olp2]
#	
#	mat3 <- mat1_olp*mat2_olp
#	
#	compositeSigList[[i]] <- mat3
#}












for(i in 1:length(compositeSigList))
{
	compositeSigList[[i]] <- compositeSigList[[i]] [rowSums(compositeSigList[[i]]) > 0,]
	print(dim(compositeSigList[[i]]))
	
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
	if(file.info(paste0(signatureCombPath,"combSig_mean_SP_s0_r3/",gene,"_",version,".tsv"))$size>1)
	{
		sigmat <- read.table(paste0(signatureCombPath,"combSig_mean_SP_s0_r3/",gene,"_",version,".tsv"), sep="\t", check.names=F)
		#sigmat <- data.table::fread(paste0(signatureCombPath,"combSig_mean_SP_s0_r3/",gene,"_",version,".tsv"), sep = "\t",check.names=F)
		#setnames(sigmat, 1, "SampleID") 
		#rownames(sigmat) <- c(sigmat[,1])
		#sigmat <- sigmat[,-1,drop=F]
		#sigmat <- as.matrix(sigmat)
		
		if(ncol(sigmat)>0) 
		{
			for(compSig in names(compositeSigList)[18])
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



for(compSig in names(compositeSigList)[18])
{
	genes <- list.files(paste0(signatureCombPath,compSig))
	genes <- genes[grepl(paste0(version,".csv"),genes)]
	genes <- gsub(paste0("_",version,".csv"),"",genes)
	genes <- unique(genes)
	
	genes <- genes[ !grepl("IGH",genes) ]
	genes <- genes[ !grepl("IGK",genes) ]
	genes <- genes[ !grepl("IGL",genes) ]
	
	genes <- genes[ !grepl("MT-",genes,fixed=T) ]
	
	Receptors <- readLines("/data/rub2/project/Secretome/results/QCSummary/Receptor_genes")
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
		#if(m==3)
		#{
			#smry_m <- as.matrix(smry[, stat[,"dsCount"] >= 3 | (stat[,"dsCount"] >= 2 & stat[,"sigCount"] >=5)])
		#}else{
			smry_m <- as.matrix(smry[, stat[,"dsCount"] >= m ])
		#}
		
		smry_m <- smry_m[apply(smry_m,1,function(x) sum(is.na(x)))==0,]
	
		write.table(smry_m,paste0(signatureCombPath,compSig,"_ds",m,"_1k_",version,".tsv"),quote=F,sep="\t")
		print(dim(smry_m))
		
		
		
		cor_cutoff <- 0.9
		mat <- smry_m
		dis <- as.dist(1-cor(mat,method="pearson"))
		hc <- hclust(dis,method="complete")
		
		group_labels <- cutree(hc, h = 1-cor_cutoff)
			
		newsig <- data.frame()
		
		for(j in unique(group_labels))
		{
			geneGroups <- names(group_labels)[group_labels==j]
			
			newsig[rownames(mat),paste0(geneGroups,collapse="|")] <- rowMeans(mat[,geneGroups,drop=F])
		}
		
		print(dim(newsig))
		write.table(newsig,paste0(signatureCombPath,compSig,"_ds",m,"_1k_",version,"_",cor_cutoff,".tsv"),quote=F,sep="\t")
		
	}
	
	
}	

#########################################HLA.E


# 0 CytoSig
# 1 1 sample
# 2 2 sample
# 3 3 sample
# 4 4 sample
# 5 5 sample
# 6 0.05 FDR
# 7 0.10 FDR
# 8 0.15 FDR
# 9 0.20 FDR
# 10 0.25 FDR
# 11 NicheNet
# 112 NicheNet v2
# 12 0.25 FDR, < 3 datasets
# 13 all mean-centrolized
# 14 MoranI mean-centrolized
# 15 MoranI+TCGA mean-centrolized
# 16 MoranI+TCGA z scale
# 17 MoranI mean-centrolized for filtered SPs
# 18 MoranI+TCGA mean-centrolized for filtered SPs
# 19 exponential decay + SP in sig
# 20 immuneDic
# 37 secact new add two more cancer types from TCGA and ICGC in filtering step

# 618 618 samples filter by TCGA and ICGC
# 6182 618 samples filter by TCGA only


