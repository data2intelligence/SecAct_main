source("Secretome_s0_path.R")
source("../../SpaCE/code/ifun.R")
source("importHPA.R")

#SWARM -t 50 -g 100 --time 01:30:00
args = commandArgs(trailingOnly=TRUE)
st <- args[1]
sampleName <- args[2]

sts <- unique(meta[,1])
#for(st in sts)
#{
	preprocessDataPath.st <- paste0(preprocessDataPath,st,"/")
	
	signaturePath.st <- paste0(signaturePath,st,"/")
	dir.create(signaturePath.st)
	
	deconvResPath.st <- paste0(deconvResPath,st,"/")
	dir.create(deconvResPath.st)
	
	sampleNames <- meta[meta[,"Study"]==st,"Sample_Name"]
	#for(sampleName in sampleNames)
	#{	
		signaturePath.st.sample <- paste0(signaturePath.st,sampleName,"/")
		dir.create(signaturePath.st.sample)
		
		deconvResPath.st.sample <- paste0(deconvResPath.st,sampleName,"/")
		dir.create(deconvResPath.st.sample)
		
		# Version 1 condition logUMI
		if(FALSE)
		{
			st.matrix.data <- as.matrix(read.table(paste0(preprocessDataPath.st,sampleName,"_counts.tsv"),check.names=FALSE))
			rownames(st.matrix.data) <- transferSymbol(rownames(st.matrix.data))
			st.matrix.data <- rm_duplicates(st.matrix.data)
			
			set.seed(123456)
			st.matrix.data.vst <- sctransform::vst(st.matrix.data, min_cells=5)$y
			st.matrix.data.vst <- round(st.matrix.data.vst,3)
			
			spots <- colnames(st.matrix.data.vst)
			coordi <- round(t(matrix(as.numeric(unlist(strsplit(spots,"x"))),nrow=2)),0)
			rownames(coordi) <- spots
			coordi <- coordi[order(coordi[,1],coordi[,2]),]
			
			st.matrix.data.vst.sorted <- st.matrix.data.vst[,rownames(coordi)]
			write.table(st.matrix.data.vst.sorted, paste0(signaturePath.st.sample,"vst.tsv"), quote=F, sep="\t")
			
			runMoranI <- 
				paste0(
				"./pairwise_moran_I/Release/pairwise_moran_I",
				" -i ",
				paste0(signaturePath.st.sample,"vst.tsv"),
				" -o ",
				paste0(signaturePath.st.sample,"spatial_sig_vst_inhouse_s0_r3.tsv"),
				" -s 0",
				" -r 3"
			)
			system(runMoranI)
		}
		
		
		# Version 2 vst condition logUMI+cellType
		if(FALSE)
		{
			st.matrix.data <- as.matrix(read.table(paste0(preprocessDataPath.st,sampleName,"_counts.tsv"),check.names=FALSE))
			rownames(st.matrix.data) <- transferSymbol(rownames(st.matrix.data))
			st.matrix.data <- rm_duplicates(st.matrix.data)
			st.matrix.data <- st.matrix.data[,colSums(st.matrix.data)>100]
			
			propMat <- as.matrix(read.csv(paste0(deconvResPath.st.sample,"SpaCET_res.csv"),row.names=1,check.names=FALSE))
			
			## option 1
			#propMat <- propMat[1:(which(rownames(propMat)=="Unidentifiable")-1),]
			#propMat <- propMat[rowMeans(propMat)>0.01,]
			#meta_data <- as.data.frame(t(propMat))
			#colnames(meta_data) <- gsub(" ","",colnames(meta_data))		
			#colnames(meta_data) <- gsub("y","",colnames(meta_data)) # there is a bug in vst, when a word with y, vst will remove it.
			
			# option 2
			meta_data <- data.frame()
			meta_data[colnames(st.matrix.data),"logUMI"] <- log1p(colSums(st.matrix.data))  # log(UMI + 1)
			meta_data$Malignant <- round(propMat["Malignant",],3)
			meta_data$Stromal <- round(colSums(propMat[c("CAF","Endothelial"),]),3)
			meta_data$ImmuneL <- round(colSums(propMat[c("B cell","T CD4","T CD8","NK","Plasma"),]),3)
			meta_data$ImmuneM <- round(colSums(propMat[c("cDC","pDC","Macrophage","Mast","Neutrophil"),]),3)
			meta_data <- meta_data[,c(TRUE,colMeans(meta_data[,2:5])>0.01),drop=F]
			
			set.seed(123456)
			st.matrix.data.vst <- sctransform::vst(st.matrix.data, cell_attr = meta_data, latent_var = colnames(meta_data), method="glmGamPoi", min_cells=5)$y
			st.matrix.data.vst <- round(st.matrix.data.vst,3)
			
			write.table(st.matrix.data.vst, paste0(signaturePath.st.sample,"vst_condition_logUMI_cellType.tsv"), quote=F, sep="\t")
		}
		
		
		# Version 3 vst condition cellType
		if(TRUE)
		{
			st.matrix.data <- as.matrix(read.table(paste0(preprocessDataPath.st,sampleName,"_counts.tsv"),check.names=FALSE))
			rownames(st.matrix.data) <- transferSymbol(rownames(st.matrix.data))
			st.matrix.data <- rm_duplicates(st.matrix.data)
			st.matrix.data <- st.matrix.data[,colSums(st.matrix.data)>100]
			
			propMat <- as.matrix(read.csv(paste0(deconvResPath.st.sample,"SpaCET_res.csv"),row.names=1,check.names=FALSE))
			
			# option 2
			meta_data <- data.frame()
			meta_data[colnames(st.matrix.data),"Malignant"] <- 0
			meta_data$Malignant <- round(propMat["Malignant",],3)
			meta_data$Stromal <- round(colSums(propMat[c("CAF","Endothelial"),]),3)
			meta_data$ImmuneL <- round(colSums(propMat[c("B cell","T CD4","T CD8","NK","Plasma"),]),3)
			meta_data$ImmuneM <- round(colSums(propMat[c("cDC","pDC","Macrophage","Mast","Neutrophil"),]),3)
			meta_data <- meta_data[,colMeans(meta_data)>0.01,drop=F]
			
			set.seed(123456)
			st.matrix.data.vst <- sctransform::vst(st.matrix.data, cell_attr = meta_data, latent_var = colnames(meta_data), method="glmGamPoi", min_cells=5)$y
			st.matrix.data.vst <- round(st.matrix.data.vst,3)
			
			write.table(st.matrix.data.vst, paste0(signaturePath.st.sample,"vst_condition_cellType.tsv"), quote=F, sep="\t")
		}
		
		
		
		
		
		#set.seed(123456)
		#st.matrix.data.vst <- sctransform::vst(st.matrix.data, min_cells=5)$y
		#st.matrix.data.vst <- round(st.matrix.data.vst,3)
		#
		#W <- calWeights(colnames(st.matrix.data.vst), r=3, diag0=TRUE)
		#st.matrix.data.vst <- st.matrix.data.vst[,colnames(W)]
		#
		#Moran_I_paired_genes <- spatialCrossCorrelation(st.matrix.data.vst, W)
		#
		#
		#
		#propMat <- as.matrix(read.csv(paste0(deconvResPath.st.sample,"SpaCET_res.csv"),row.names=1,check.names=FALSE))
		#propMat <- propMat[1:(which(rownames(propMat)=="Unidentifiable")-1),]
		#propMat <- propMat[rowMeans(propMat)>0.01,]
		#
		#olp <- intersect(colnames(st.matrix.data.vst), colnames(propMat))
		#
		#weight_mat <- WGCNA::cor(t(st.matrix.data.vst[,olp,drop=F]),t(propMat[,olp,drop=F]) )
		#weight_vec <- apply(weight_mat, 1, function(x) max(x,na.rm = TRUE))
		#weight_vec <- 1 - weight_vec
		#
		#Moran_I_paired_genes_weighted <- sweep(Moran_I_paired_genes, 1, weight_vec, "*")
		#Moran_I_paired_genes_weighted <- round(Moran_I_paired_genes_weighted, 5)
		#
		#
		#write.csv(Moran_I_paired_genes, paste0(signaturePath.st.sample,"Moran_I_paired_genes.tsv"), quote=F)
		#saveRDS(Moran_I_paired_genes, file = paste0(signaturePath.st.sample,"Moran_I_paired_genes.rds"))
		#
		#
		#
		#
		#
		#
		#
		#
		#
		#f2 <- read.table("/data/rub2/project/Secretome/code/CytoSig-master/CytoSig/signature.centroid",row.names=1,header=T)
		#colnames(f2) <- transferSymbol(colnames(f2))
		#
		#
		#
		#st.matrix.data <- as.matrix(read.table(paste0(preprocessDataPath.st,sampleName,"_counts.tsv"),check.names=FALSE))
		#rownames(st.matrix.data) <- transferSymbol(rownames(st.matrix.data))
		#st.matrix.data <- rm_duplicates(st.matrix.data)
		#st.matrix.data <- st.matrix.data[rowSums(st.matrix.data>0)>=5,]
		#
		#
		#meta_data <- data.frame()
		#meta_data[colnames(st.matrix.data),"logUMI"] <- log1p(colSums(st.matrix.data))  # log(UMI + 1)
		#
		#propMat <- as.matrix(read.csv(paste0(deconvResPath.st.sample,"SpaCET_res.csv"),row.names=1,check.names=FALSE))
		#propMat <- propMat[1:13,]
		#propMat <- propMat[rowMeans(propMat)>0.01,]
		#propMat <- t(propMat)
		#
		#propMat_s <- scale(propMat)
		#colnames(propMat_s) <- paste0(colnames(propMat_s),"_s")
		#
		#meta_data <- cbind(meta_data, propMat)
		#meta_data <- cbind(meta_data, propMat_s)
		#colnames(meta_data) <- gsub(" ","",colnames(meta_data))
		
		
		#meta_data$Malignant <- round(propMat["Malignant",],3)
		#meta_data$Stromal <- round(colSums(propMat[c("CAF","Endothelial"),]),3)
		#meta_data$ImmuneL <- round(colSums(propMat[c("B cell","T CD4","T CD8","NK","Plasma"),]),3)
		#meta_data$ImmuneM <- round(colSums(propMat[c("cDC","pDC","Macrophage","Mast","Neutrophil"),]),3)
		#meta_data$CAF <- round(propMat["CAF",],3)
		#meta_data$Macrophage <- round(propMat["Macrophage",],3)
		
		
		
		#set.seed(123456)
		#st.matrix.data.log <- t(t(st.matrix.data)*1e5/colSums(st.matrix.data))
		#st.matrix.data.vst <- log2(st.matrix.data.log + 1)
		#st.matrix.data.vst <- round(st.matrix.data.vst,3)
		#
		#W <- calWeights(colnames(st.matrix.data.vst), r=3, diag0=TRUE)
		#st.matrix.data.vst <- st.matrix.data.vst[,colnames(W)]
		#Moran_I_paired_genes <- spatialCrossCorrelation(st.matrix.data.vst, W)
		#write.csv(Moran_I_paired_genes["TGFB1",], paste0(signaturePath.st.sample,"TGFB1_log.tsv"), quote=F)
		#
		#olp_targets <- intersect(rownames(f2),rownames(Moran_I_paired_genes))
		#olp_SP <- intersect(colnames(f2),colnames(Moran_I_paired_genes))
		#cors <- diag(cor(f2[olp_targets,olp_SP], Moran_I_paired_genes[olp_targets,olp_SP]))
		#write.csv(cors, paste0(signaturePath.st.sample,"CytoSig_log.tsv"), quote=F)
		
		
		
		#set.seed(123456)
		#st.matrix.data.vst <- sctransform::vst(st.matrix.data, min_cells=5)$y
		#st.matrix.data.vst <- round(st.matrix.data.vst,3)
		
		#W <- calWeights(colnames(st.matrix.data.vst), r=3, diag0=TRUE)
		#st.matrix.data.vst <- st.matrix.data.vst[,colnames(W)]
		#Moran_I_paired_genes <- spatialCrossCorrelation(st.matrix.data.vst, W)
		#write.csv(Moran_I_paired_genes["TGFB1",], paste0(signaturePath.st.sample,"TGFB1_vst_umi_default.tsv"), quote=F)
		#
		#olp_targets <- intersect(rownames(f2),rownames(Moran_I_paired_genes))
		#olp_SP <- intersect(colnames(f2),colnames(Moran_I_paired_genes))
		#cors <- diag(cor(f2[olp_targets,olp_SP], Moran_I_paired_genes[olp_targets,olp_SP]))
		#write.csv(cors, paste0(signaturePath.st.sample,"CytoSig_vst_umi_default.tsv"), quote=F)
		
		
		#weight_mat <- WGCNA::cor(t(st.matrix.data.vst),propMat)
		#weight_vec <- apply(weight_mat, 1, max)
		#weight_vec <- 1 - weight_vec^2
		
		
		#Moran_I_paired_genes <- Moran_I_paired_genes * weight_vec
		#
		#olp_targets <- intersect(rownames(f2),rownames(Moran_I_paired_genes))
		#olp_SP <- intersect(colnames(f2),colnames(Moran_I_paired_genes))
		#cors <- diag(cor(f2[olp_targets,olp_SP], Moran_I_paired_genes[olp_targets,olp_SP]))
		#write.csv(cors, paste0(signaturePath.st.sample,"CytoSig_vst_umi_default_weighted.tsv"), quote=F)
		#
		#
		#
		#
		#
		#
		#
		#Moran_I_paired_genes["TGFB1",]
		#
		#old <- Moran_I_paired_genes["TGFB1",]
		#new <- Moran_I_paired_genes["TGFB1",]*weight_vec
		#
		#
		#head(sort(old, decreasing=T),100)
		#head(sort(new, decreasing=T),100)
		#
		#cor(f2[olp_targets,"TGFB1"],old[olp_targets])
		#cor(f2[olp_targets,"TGFB1"],new[olp_targets])
		
		#spot_coordinates <- t(matrix(as.numeric(unlist(strsplit(colnames(st.matrix.data),"x"))),nrow=2))
		#colnames(spot_coordinates) <- c("x","y")
		#rownames(spot_coordinates) <- colnames(st.matrix.data)
		#spot_coordinates <- data.frame(spot_coordinates)
		#
		#library(SPARK)
		#spark_obj <- CreateSPARKObject(
		#  counts = st.matrix.data[c(rownames(st.matrix.data)[1:100],"C1QA","C1QB","C1QC","TGFB1","OLR1","CTSL","RND1","GADD45B","TGFBR1","POSTN","IGFL1","SMAD7","COMP","SERPINE1"),],
		#  location = spot_coordinates,  # matrix of x, y for each spot
		#  percentage = 0.1,             # keep top genes
		#  min_total_counts = 0
		#)
		#
		## Fit model (no covariates here, since VST already removed cell type effects)
		#spark_obj <- spark.vc(spark_obj, covariates=propMat, num_core = 18)
		#
		## Test spatial variability
		#spark_obj <- spark.test(spark_obj)
		#
		#
		#write.csv(spark_obj@res_mtest, paste0(signaturePath.st.sample,"spark.tsv"), quote=F)
		
		
		
#		W <- calWeights(colnames(st.matrix.data.vst), r=3, diag0=TRUE)
#		st.matrix.data.vst <- st.matrix.data.vst[,colnames(W)]
#		Moran_I_paired_genes <- spatialCrossCorrelation(st.matrix.data.vst, W)
#		write.csv(Moran_I_paired_genes["TGFB1",], paste0(signaturePath.st.sample,"TGFB1_vst_umi_default_spark.tsv"), quote=F)
#		
#		olp_targets <- intersect(rownames(f2),rownames(Moran_I_paired_genes))
#		olp_SP <- intersect(colnames(f2),colnames(Moran_I_paired_genes))
#		cors <- diag(cor(f2[olp_targets,olp_SP], Moran_I_paired_genes[olp_targets,olp_SP]))
#		write.csv(cors, paste0(signaturePath.st.sample,"CytoSig_vst_umi_default_spark.tsv"), quote=F)
#		
#		
#		
#		
#		
#		
#		library(SPARK)
#		spark_obj <- CreateSPARKObject(
#		  counts = vst_output$umi,
#		  location = spot_coordinates,  # matrix of x, y for each spot
#		  percentage = 0.1,             # keep top genes
#		  min_total_counts = 10
#		)
#		
#		# Fit model (no covariates here, since VST already removed cell type effects)
#		spark_obj <- spark.vc(spark_obj, num_core = 4)
#		
#		# Test spatial variability
#		spark_obj <- spark.test(spark_obj)
#		
#		
#		
#		
#		
#		latent_var_list <- list(
#			vst_highCellFraction=colnames(meta_data)[2:9],
#			vst_highCellFractionScale=colnames(meta_data)[10:17]
#			#vst_Malignant=c("Malignant"),
#			#vst_Malignant_s=c("Malignant_s"),
#			#vst_Stromal=c("Stromal"),
#			#vst_ImmuneL=c("ImmuneL"),
#			#vst_ImmuneM=c("ImmuneM"),
#			#vst_cellFraction=c("Malignant","Stromal","ImmuneL","ImmuneM"),
#			#vst_cellFraction_umi_my=c("logUMI","Malignant","Stromal","ImmuneL","ImmuneM"),
#			#vst_Macrophage_umi_my=c("logUMI","Macrophage"),
#			#vst_Macrophage=c("Macrophage"),
#			#vst_CAF_umi_my=c("logUMI","CAF")
#		)
#		
#		for(i in 1:length(latent_var_list))
#		{
#			set.seed(123456)
#			st.matrix.data.vst <- sctransform::vst(st.matrix.data, cell_attr = meta_data, latent_var = latent_var_list[[i]], min_cells=5, method="glmGamPoi")$y
#			st.matrix.data.vst <- round(st.matrix.data.vst,3)
#			
#			W <- calWeights(colnames(st.matrix.data.vst), r=3, diag0=TRUE)
#			st.matrix.data.vst <- st.matrix.data.vst[,colnames(W)]
#			Moran_I_paired_genes <- spatialCrossCorrelation(st.matrix.data.vst, W)
#			write.csv(Moran_I_paired_genes["TGFB1",], paste0(signaturePath.st.sample,"TGFB1_",names(latent_var_list)[i],".tsv"), quote=F)
#			
#			olp_targets <- intersect(rownames(f2),rownames(Moran_I_paired_genes))
#			olp_SP <- intersect(colnames(f2),colnames(Moran_I_paired_genes))
#			cors <- diag(cor(f2[olp_targets,olp_SP], Moran_I_paired_genes[olp_targets,olp_SP]))
#			write.csv(cors, paste0(signaturePath.st.sample,"CytoSig_",names(latent_var_list)[i],".tsv"), quote=F)
#		}
		
		
		
		
		
		
#		gene <- "TGFB1"
#		
#		
#		if(gene%in%rownames(Moran_I_paired_genes))
#		{
#			write.csv(Moran_I_paired_genes[gene,], paste0(signaturePath.st.sample,gene,"_raw.tsv"), quote=F)
#		}
#		
#		
#		propMat <- as.matrix(read.csv(paste0(deconvResPath.st.sample,"SpaCET_res.csv"),row.names=1,check.names=FALSE))
#		#non_CAF_spots <- colnames(propMat)[!propMat[c("Macrophage","Endothelial","cDC"),] > 0.3]
#		
#		for(target in c("SERPINE1","LTBP2"))
#		{
#			non_CAF_spots <- names(st.matrix.data.vst[target,])[st.matrix.data.vst[target,]<0]
#			
#			
#			W_CAF <- W
#			W_CAF[,colnames(W_CAF)%in%non_CAF_spots] <- 0
#			
#			Moran_I_paired_genes <- spatialCrossCorrelation(st.matrix.data.vst, W_CAF)
#			
#			if(gene%in%rownames(Moran_I_paired_genes))
#			{
#				write.csv(Moran_I_paired_genes[gene,], paste0(signaturePath.st.sample,gene,"_",target,".tsv"), quote=F)
#			}
#		
#		}
#		
		#corMat <- cor(propMat)
		#diag(corMat) <- 0
		#
		#logical_W <- W > 0
		#logical_W[logical_W == FALSE] <- NA
		#
		#corMat_W <- corMat*logical_W
		#corMat_W_vec <- c(corMat_W)
		#corMat_W_vec <- corMat_W_vec[!is.na(corMat_W_vec)]
		
#> quantile(corMat_W_vec,probs = seq(0, 1, 0.05))
#        0%         5%        10%        15%        20%        25%        30% 
#-0.1030399  0.6544652  0.7518672  0.8085847  0.8481711  0.8780432  0.9031741 
#       35%        40%        45%        50%        55%        60%        65% 
# 0.9240329  0.9408678  0.9545844  0.9660453  0.9747229  0.9813410  0.9863157 
#       70%        75%        80%        85%        90%        95%       100% 
# 0.9901633  0.9930020  0.9951858  0.9968227  0.9982096  0.9993469  1.0000000 
		
		
		
		
		
		
		
		## SP
		#sigmat <- Moran_I_paired_genes[,colnames(Moran_I_paired_genes)%in%SPs]		
		#sigmat <- sigmat - rowMeans(sigmat,na.rm=TRUE)
		#fname <- "singleSig_mean_centralized_SP_s0_r3"
		#
		#signaturePath.st.sample.singleSig <- paste0(signaturePath.st,sampleName,"/",fname,"/")
		#dir.create(signaturePath.st.sample.singleSig)
		#
		#for(gene in colnames(sigmat))
		#{
		#	write.table(sigmat[,gene,drop=F], paste0(signaturePath.st.sample.singleSig,gene,".tsv"), quote=F, sep="\t")
		#}
		#
		#
		## MP
		#sigmat <- Moran_I_paired_genes[,colnames(Moran_I_paired_genes)%in%MPs]		
		#sigmat <- sigmat - rowMeans(sigmat,na.rm=TRUE)
		#fname <- "singleSig_mean_centralized_MP_s0_r3"
		#
		#signaturePath.st.sample.singleSig <- paste0(signaturePath.st,sampleName,"/",fname,"/")
		#dir.create(signaturePath.st.sample.singleSig)
		#
		#for(gene in colnames(sigmat))
		#{
		#	write.table(sigmat[,gene,drop=F], paste0(signaturePath.st.sample.singleSig,gene,".tsv"), quote=F, sep="\t")
		#}
		#
		#
		## compute self Moran I
		#Moran_I_self_genes <- spatialAutoCorrelation(st.matrix.data.vst[rownames(st.matrix.data.vst)%in%SPs,], W)
		#
		#write.csv(Moran_I_self_genes, paste0(signaturePath.st.sample,"MoranI_SPs.csv"), quote=F)

		
#	}
#}




if(FALSE)
{

source("Secretome_s0_path.R")
source("../../SpaCE/code/ifun.R")
source("importHPA.R")

sts <- unique(meta[,1])
for(st in sts)
{
	preprocessDataPath.st <- paste0(preprocessDataPath,st,"/")
	
	signaturePath.st <- paste0(signaturePath,st,"/")
	dir.create(signaturePath.st)
	
	deconvResPath.st <- paste0(deconvResPath,st,"/")
	dir.create(deconvResPath.st)
	
	QCPath.st <- paste0(QCPath,st,"/")
	dir.create(QCPath.st)
	
	QCFilterPath.st <- paste0(QCFilterPath,st,"/")
	dir.create(QCFilterPath.st)
		
	sampleNames <- meta[meta[,"Study"]==st,"Sample_Name"]
	for(sampleName in sampleNames)
	{	
		signaturePath.st.sample <- paste0(signaturePath.st,sampleName,"/")
		dir.create(signaturePath.st.sample)
		
		deconvResPath.st.sample <- paste0(deconvResPath.st,sampleName,"/")
		dir.create(deconvResPath.st.sample)
		
		QCPath.st.sample <- paste0(QCPath.st,sampleName,"/")
		dir.create(QCPath.st.sample)
		
		QCFilterPath.st.sample <- paste0(QCFilterPath.st,sampleName,"/")
		dir.create(QCFilterPath.st.sample)
		
		
		#if(!file.exists(paste0(signaturePath.st.sample,"vst_condition_cellType.tsv"))) print(paste0(st,"@",sampleName))
		#if(!file.exists(paste0(signaturePath.st.sample,"spatial_sig_vst_inhouse_s0_r3.tsv"))) print(paste0(st,"@",sampleName))
		#if(!file.exists(paste0(deconvResPath.st.sample,"cor.txt"))) print(paste0(st,"@",sampleName))
		#if(!file.exists(paste0(deconvResPath.st.sample,"SpaCET_res.csv"))) print(paste0(st,"@",sampleName))
		#if(!file.exists(paste0(deconvResPath.st.sample,"weighted.tsv"))) print(paste0(st,"@",sampleName))
		#if(!file.exists(paste0(QCPath.st.sample,"moranI.fast.csv"))) print(paste0(st,"@",sampleName))
		#if(!file.exists(paste0(QCPath.st.sample,"vst_condition_logUMI_cellType_moranI.fast.csv"))) print(paste0(st,"@",sampleName))
		#if(!file.exists(paste0(QCPath.st.sample,"ICGC_filter_RECA-EU_vst_condition_cellType.csv"))) print(paste0(st,"@",sampleName))
		if(!file.exists(paste0(QCFilterPath.st.sample,"TCGA_Xena_filter_log_Secreted_filterBar_vst_condition_cellType_FDR.png"))) print(paste0(st,"@",sampleName))
		
		
		#runMoranI <- 
		#	paste0(
		#	"mv ",
		#	paste0(signaturePath.st.sample ,"/singleSig_mean_SP_s0_r3"),
		#	" ",
		#	paste0(signaturePath.st.sample ,"/singleSig_mean_SP_s0_r3_vst")
		#)
		#system(runMoranI)
	}		
}		




sts <- unique(meta[meta[,"Method"]=="Visium",1])
for(st in sts)
{
	signaturePath.st <- paste0(signaturePath,st,"/")
	system(paste0("mkdir /data/Jiang_Lab/Data/Beibei/ST/MoranI/", st))
	
	sampleNames <- meta[meta[,"Study"]==st&meta[,"Method"]=="Visium"&meta[,"Organism"]=="Human","Sample_Name"]
	for(sampleName in sampleNames)
	{	
		signaturePath.st.sample <- paste0(signaturePath.st,sampleName,"/")
		system(paste0("mkdir /data/Jiang_Lab/Data/Beibei/ST/MoranI/", st, "/",sampleName))
		
		system(paste0("cp ", signaturePath.st.sample, "*tsv /data/Jiang_Lab/Data/Beibei/ST/MoranI/", st, "/",sampleName))
		
	}
}



}