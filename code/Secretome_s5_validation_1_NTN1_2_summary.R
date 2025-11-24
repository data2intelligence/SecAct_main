source("Secretome_s0_path.R")

inputPath <- paste0(dataValPath,"NTN1/")
outputPath <- paste0(validationPath,"NTN1/")

# precess UCEC ST data

library(SpaCET)
library(SecAct)

sampleNames <- c("01-034_C1D1","01-034_C3D1","01-039_C1D1","01-039_C3D1")

###########
# fig s4b #
###########

for(sampleName in sampleNames)
{
	if(sampleName=="01-039_C1D1"){
		widthValue=8
		heightValue=8.2
	}else if(sampleName=="01-039_C3D1"){
		widthValue=10
		heightValue=10
	}else if(sampleName=="01-034_C1D1"){
		widthValue=7.8
		heightValue=7.8
	}else if(sampleName=="01-034_C3D1"){
		widthValue=9.8
		heightValue=9.6
	}else{
		widthValue=14
		heightValue=12
	}
	
	SpaCET_obj <- readRDS(paste0(outputPath,sampleName,".rds"))
	
	p1 <- SpaCET.visualize.spatialFeature(
  			  SpaCET_obj,
  			  spatialType = "CellFraction",
  			  spatialFeatures=c("Malignant"),
  			  imageBg=FALSE
  			)
  	p1 <- p1 + theme(
  		panel.background = element_blank(),
  		plot.background = element_rect(fill = "transparent", color = NA)
  	)
  	
	ggsave(paste0(outputPath,sampleName,"_malignant.png"), p1, width = widthValue, height = heightValue, dpi=200, units = "cm",limitsize = FALSE)

	write.csv(SpaCET_obj@results$deconvolution$propMat["Malignant",],paste0(outputPath,sampleName,"_malignant.csv"),quote=F)
}


######################
# fig 3b and fig s4a #
######################

for(sampleName in sampleNames)
{
	if(sampleName=="01-039_C1D1"){
		widthValue=8.5
		heightValue=8.2
	}else if(sampleName%in%c("01-039_C3D1")){
		widthValue=10
		heightValue=10
	}else if(sampleName=="01-034_C1D1"){
		widthValue=7.8
		heightValue=7.8
	}else if(sampleName%in%c("01-034_C3D1")){
		widthValue=9.8
		heightValue=9.6
	}else{
		widthValue=14
		heightValue=12
	}
	
	SpaCET_obj <- readRDS(paste0(outputPath,sampleName,".rds"))
	write.csv(SpaCET_obj@results$SecAct_output$SecretedProteinActivity$zscore["NTN1",],paste0(outputPath,sampleName,"_NTN1_act.csv"),quote=F)

	act <- SpaCET_obj@results$SecAct_output$SecretedProteinActivity$zscore
	act["NTN1",act["NTN1",]< -2] <- -2
	act["NTN1",act["NTN1",]>6] <- 6
	act -> SpaCET_obj@results$SecAct_output$SecretedProteinActivity$zscore
	
	p1 <- SpaCET.visualize.spatialFeature(
  			  SpaCET_obj,
  			  spatialType = "SecretedProteinActivity",
  			  spatialFeatures=c("NTN1"),
  			  imageBg=FALSE,
  			  colors=c("#b8e186","#c51b7d")
  			)

	ggsave(paste0(outputPath,sampleName,"_NTN1_act.png"), p1, width = widthValue, height = heightValue, dpi=200, units = "cm",limitsize = FALSE)

	p1 <- SpaCET.visualize.spatialFeature(
  			  SpaCET_obj,
  			  spatialType = "GeneExpression",
  			  spatialFeatures=c("NTN1"),
  			  imageBg=FALSE)

	ggsave(paste0(outputPath,sampleName,"_NTN1_expr.png"), p1, width = widthValue, height = heightValue, dpi=200, units = "cm",limitsize = FALSE)
}




###########
# fig 3c #
###########

for(patient in c("01-034","01-039"))
{
	SpaCET_obj1 <- readRDS(paste0(outputPath,patient,"_C1D1.rds"))
	SpaCET_obj2 <- readRDS(paste0(outputPath,patient,"_C3D1.rds"))
	
	if(patient=="01-039")
	{
		xpos = 1.8
		ypos = 5.5
		yposn = -8.8
	}
	if(patient=="01-034")
	{
		xpos = 1.8
		ypos = 5
		yposn = -6.8
	}
	
	
	pre <- data.frame(
		Treatment="Pre",
		Activity=SpaCET_obj1@results$SecAct_output$SecretedProteinActivity$zscore["NTN1",],
		Fraction=SpaCET_obj1@results$deconvolution$propMat["Malignant",]
	)
	on <- data.frame(
		Treatment="On",
		Activity=SpaCET_obj2@results$SecAct_output$SecretedProteinActivity$zscore["NTN1",],
		Fraction=SpaCET_obj2@results$deconvolution$propMat["Malignant",]
	)
	
	fg.df <- rbind(pre,on)
	fg.df <- fg.df[fg.df[,3]>0.5,]
	fg.df[,1] <- factor(fg.df[,1], levels=c("Pre","On"))

	wilcox_res <- wilcox.test(fg.df[fg.df[,1]=="Pre",2], fg.df[fg.df[,1]=="On",2])
	pv <- signif(wilcox_res$p.value,2)
	
	library(ggplot2)
	p2 <- ggplot(fg.df,aes(x=Treatment,y=Activity)) + 
		geom_violin(aes(group=Treatment, fill=Treatment),trim=FALSE)+
		geom_boxplot(width=0.15,outlier.shape=NA)+
		scale_fill_manual( values=c("#f2cecf","#c8bedf") )+
		annotate("text", x = xpos, y=ypos, label = paste0("p = ",pv))+
		#annotate("text", x = 1, y=yposn, label = paste0("n = ",table(fg.df[,1])[1]))+
		#annotate("text", x = 2, y=yposn, label = paste0(table(fg.df[,1])[2]))+
		scale_x_discrete(labels= c(paste0("Pre\n",table(fg.df[,1])[1]),paste0("On\n",table(fg.df[,1])[2])))+
		ylab("NTN1 Activity")+
		xlab("Treatment")+
		theme_classic()+ 
		theme(
			panel.grid = element_blank(),
	  		panel.background = element_blank(),
	  		plot.title = element_text(hjust = 0.5),
    		axis.title = element_text(colour = "black"),
			axis.text = element_text(colour = "black", size=12),
			axis.title.x = element_blank(),
			legend.position="none"
		)
	ggsave(paste0(outputPath,patient,"_compare.png"), p2, width = 5.8, height = 6.2, dpi=300, units = "cm",limitsize = FALSE)
	write.csv(fg.df,paste0(outputPath,patient,"_compare.csv"),quote=F)
}







if(FALSE)
{
	# precess UCEC bulk RNAseq data

	NTN_bulk <- as.matrix(read.csv(paste0(inputPath,"GSE225687_Patient_Raw_Count.csv"),sep=";",row.names=1,header=T))
	rownames(NTN_bulk) <- substr(rownames(NTN_bulk),1,15)
	NTN_bulk <- rm_duplicates(NTN_bulk)
	
	gene_anno <- as.matrix(read.csv(gzfile(paste0(inputPath,"GSE225689_RAW/C1D1/features.tsv.gz")),as.is=T,header=F,sep="\t"))
	
	olp <- intersect(rownames(NTN_bulk),gene_anno[,1])
	
	NTN_bulk <- NTN_bulk[olp,]
	rownames(NTN_bulk) <- gene_anno[match(rownames(NTN_bulk),gene_anno[,1]),2]
	rownames(NTN_bulk) <- transferSymbol(rownames(NTN_bulk))
	NTN_bulk <- rm_duplicates(NTN_bulk)
	
	expr <- NTN_bulk
	
	expr <- expr[rowSums(expr>0)>=5,]
	expr.scaled <- t(t(expr)*1e6/colSums(expr))
	expr.scaled.log <- log2(expr.scaled + 1 )
			
	f3 <- expr.scaled.log[,(1:12)*2] - expr.scaled.log[,(1:12)*2 - 1]
	
	write.table(f3,paste0(dataValPath,"blocking/NTN1_GSE225691.diff"),quote=FALSE,sep="\t")
	
	
	# Mouse_Pten data
	
	gene_anno <- as.matrix(read.csv(gzfile(paste0("/data/rub2/project/Melanoma/data/SpaceRanger_Result_Files/Sample_ST_10_EPG/outs/filtered_feature_bc_matrix/features.tsv.gz")),as.is=T,header=F,sep="\t"))
	
	NTN_bulk <- as.matrix(read.csv("/data/rub2/project/Secretome/data/GSE225688_PTEN_Raw_count.csv",sep=";",row.names=1,header=T))
	rownames(NTN_bulk) <- substr(rownames(NTN_bulk),1,18)
	NTN_bulk <- rm_duplicates(NTN_bulk)
	
	olp <- intersect(rownames(NTN_bulk),gene_anno[,1])
	
	NTN_bulk <- NTN_bulk[olp,]
	rownames(NTN_bulk) <- gene_anno[match(rownames(NTN_bulk),gene_anno[,1]),2]
	
	expr <- NTN_bulk
	
	expr.scaled <- t(t(expr)*1e6/colSums(expr))
	expr.scaled.log <- log2(expr.scaled + 1 )
	
	cdata <- expr.scaled.log
	
	
	rownames(cdata) <- mouseGenes
	cdata <- cdata[!rownames(cdata)%in%c("humanGeneNotExist","humanGeneMultiple"),]
	
	
	expr.scaled.log <- cdata
	
	
	sig_single <- expr.scaled.log[,4:6] - expr.scaled.log[,1:3]
	sig_mean <- rowMeans(sig_single)
	
	cdata_T_minusBG <- cbind(sig_single,sig_mean)
	
	library(SecAct)
		
	res <- SecAct.inference(Y=cdata_T_minusBG, SigMat=ref, lambda=1000000, nrand=1000)
		
		
	NTN_SecAct_mouse <- res$zscore
	NTN_SecAct_mouse[c("NTN1"),,drop=F]


	# scRNA-seq
	
	library(Seurat)
	library(DoubletFinder)
	
	dataPath <- "/data/rub2/project/Secretome/data/GSE225691_NTN1/GSE225689_RAW/"
	
	names <- c("C1D1","C3D1")
	list.samples.0 <- list()
	
	for(i in 1:length(names))
	{
	  list.samples.0[[i]] <- Read10X(data.dir = paste0(dataPath,names[[i]],"/"))
	  list.samples.0[[i]] <- CreateSeuratObject(list.samples.0[[i]])
	  list.samples.0[[i]]$sample <- names[[i]]
	}
	
	names(list.samples.0) <- names
	lapply(list.samples.0, dim)
	
	## $C1D1
	## [1] 36601 12385
	## 
	## $C3D1
	## [1] 36601 10168
	
	
	object.doublet <- list.samples.0
	for(i in 1:length(object.doublet)){
	  
	  temp1 <- SCTransform(object.doublet[[i]])
	  temp1 <- RunPCA(temp1)
	  temp1 <- RunUMAP(temp1, dims = 1:10)
	  
	  nExp_poi <- round(0.075*nrow(temp1@meta.data))  ## Assuming 7.5% doublet formation rate
	  
	  object.doublet[[i]] <- doubletFinder(temp1, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
	  
	  rm(temp1)
	  
	}
	
	
	
	for(i in 1:length(object.doublet)){
	  temp1 <- object.doublet[[i]]
	  colnames(temp1@meta.data)[8] = "doublets"
	  object.doublet[[i]] <- temp1
	}
	
	sc.merged <- merge(object.doublet[[1]], y = object.doublet[[2]], 
	                   add.cell.ids = c("NP3","NP4"), project = "NP137")
	
	table(sc.merged$doublets, sc.merged$sample)
	
	##           C1D1  C3D1
	##  Doublet   929   763
	##  Singlet 11456  9405
	
	library(dplyr)
	sc.merged <- PercentageFeatureSet(sc.merged, pattern = "^MT-", col.name = "percent.mt") %>%
	  PercentageFeatureSet(pattern = "^RP[SL]", col.name = "percent.ribo")
	
	sc.merged <- subset(sc.merged, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & doublets == "Singlet" & percent.mt < 25) %>%
	  SCTransform(vars.to.regress = "percent.mt") %>% RunPCA() %>% FindNeighbors(dims = 1:30) %>% 
	  FindClusters(resolution = 0.5) %>% RunUMAP(dims = 1:30) %>% RunTSNE()
	
	
	library(ggplot2)
	g <- FeaturePlot(sc.merged, reduction = "umap", features = c("EPCAM","PTPRC","ACTA2","PECAM1"), ncol =4)
	ggsave(paste0("sc_marker.jpg"), g, width = 18*4, height = 16, dpi=300, units = "cm", limitsize = FALSE)
	
	g <- DimPlot(sc.merged, reduction = "umap", group.by = "seurat_clusters",label = TRUE)
	ggsave(paste0("sc_group.by_clustering.jpg"), g, width = 18, height = 16, dpi=300, units = "cm")
	
	saveRDS(sc.merged, file = paste0("sc.merged.rds"))
	
	sc.merged <- readRDS(paste0("sc.merged.rds"))
	
	Idents(sc.merged) <-  "seurat_clusters"
	sc.merged.c <- subset(sc.merged, idents = c(0,4,8,11:13,21))
	
	cdata <- as.matrix(sc.merged.c@ assays$ SCT@ counts)
	cdata <- sctransform::vst(cdata, min_cells=5)$y
	cdata <- round(cdata,3)
	
	cdata <- cdata - rowMeans(cdata)
	
	#cdata <- t(t(cdata)*1e4/colSums(cdata))
	#cdata <- log2(cdata + 1)
	
	#cdata.c <- cdata[,1:4470]
	#cdata.t <- cdata[,4471:5581]
	#cdata.t <- cdata.t - rowMeans(cdata.c)
	
	
	library(SecAct)
	SpaCET_obj <- SecAct.inference(cdata, SigMat=ref, lambda=1e6, nrand=1000)
			
	zscore <- SpaCET_obj$zscore
	
	mean(zscore["NTN1",1:4470])
	mean(zscore["NTN1",4471:5581])
}