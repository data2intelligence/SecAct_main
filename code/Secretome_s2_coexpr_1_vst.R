source("Secretome_s0_path.R")

#SWARM -t 2 -g 30 --time 02:00:00
args = commandArgs(trailingOnly=TRUE)
st <- args[1]
sampleName <- args[2]

sts <- unique(meta[,"Study"])
#for(st in sts)
#{
	visiumPath.st <- paste0(visiumPath,st,"/")
	dir.create(visiumPath.st)
	
	signaturePath.st <- paste0(signaturePath,st,"/")
	dir.create(signaturePath.st)
	
	sampleNames <- meta[meta[,"Study"]==st,"Sample_Name"]
	#for(sampleName in sampleNames)
	#{	
		visiumPath.st.sample <- paste0(visiumPath.st,sampleName,"/")
		dir.create(visiumPath.st.sample)
		
		signaturePath.st.sample <- paste0(signaturePath.st,sampleName,"/")
		dir.create(signaturePath.st.sample)
		
		st.matrix.data <- as.matrix(read.table(gzfile(paste0(visiumPath.st.sample,"counts.tsv.gz")),check.names=FALSE))
		rownames(st.matrix.data) <- transferSymbol(rownames(st.matrix.data))
		st.matrix.data <- rm_duplicates(st.matrix.data)
		st.matrix.data <- st.matrix.data[,colSums(st.matrix.data)>100]
		
		
		
		# Version 1 vst condition logUMI
		set.seed(123456)
		st.matrix.data.vst <- sctransform::vst(st.matrix.data, min_cells=5)$y
		st.matrix.data.vst.1 <- round(st.matrix.data.vst,3)
		
		write.table(st.matrix.data.vst.1, gzfile(paste0(signaturePath.st.sample,"vst.tsv.gz")), quote=F, sep="\t")
		
		
		
		# Version 2 vst condition logUMI+cellType
		library(SpaCET)
		counts <- as(st.matrix.data, "CsparseMatrix")
		spotCoordinates <- t(matrix(as.numeric(unlist(strsplit(colnames(st.matrix.data),"x"))),nrow=2))
		rownames(spotCoordinates) <- colnames(st.matrix.data)
		colnames(spotCoordinates) <- c("array_row","array_col")
		
		SpaCET_obj <- create.SpaCET.object(counts=counts, spotCoordinates=spotCoordinates, imagePath=NA, platform="Visium")
			
		cancerType <- sapply(strsplit(st,"_",fixed=T),function(x) return(x[1]))
		if(cancerType%in%c("EPN","GIST","HB","HNAS","IPMN","MB","PCNSL","PanIN","OSCC","OPC","PNST","WT","OS","sRCC","GBC","MM"))
		{
			cancerType <- "PANCAN"
		}
		
		SpaCET_obj <- SpaCET.deconvolution(SpaCET_obj, cancerType=cancerType, coreNo=8)  
			
		propMat <- SpaCET_obj@results$deconvolution$propMat
		
		
		meta_data <- data.frame()
		meta_data[colnames(st.matrix.data),"logUMI"] <- log1p(colSums(st.matrix.data))  # log(UMI + 1)
		meta_data$Malignant <- round(propMat["Malignant",],3)
		meta_data$Stromal <- round(colSums(propMat[c("CAF","Endothelial"),]),3)
		meta_data$ImmuneL <- round(colSums(propMat[c("B cell","T CD4","T CD8","NK","Plasma"),]),3)
		meta_data$ImmuneM <- round(colSums(propMat[c("cDC","pDC","Macrophage","Mast","Neutrophil"),]),3)
		meta_data <- meta_data[,c(TRUE,colMeans(meta_data[,2:5])>0.01),drop=F]
		
		set.seed(123456)
		st.matrix.data.vst <- sctransform::vst(st.matrix.data, cell_attr = meta_data, latent_var = colnames(meta_data), method="glmGamPoi", min_cells=5)$y
		st.matrix.data.vst.2 <- round(st.matrix.data.vst,3)
		
		write.table(st.matrix.data.vst.2, gzfile(paste0(signaturePath.st.sample,"vst_condition_logUMI_cellType.tsv.gz")), quote=F, sep="\t")
		

		
		# Version 1 vs Version 2
		smyMat <- cor.two.matrix.same.rows(st.matrix.data.vst.1, st.matrix.data.vst.2)
		smyMat.median <- median(smyMat)
		
		writeLines(as.character(smyMat.median),paste0(signaturePath.st.sample,"cor_of_two_vst_version.txt"))
		
		
		
		# Version 0 log-transformed
		set.seed(123456)
		st.matrix.data <- st.matrix.data[rownames(st.matrix.data.vst.1),]
		st.matrix.data.vst <- t(t(st.matrix.data)*1e5/colSums(st.matrix.data))
		st.matrix.data.vst <- log2(st.matrix.data.vst + 1 )
		st.matrix.data.vst.0 <- round(st.matrix.data.vst,3)
		
		write.table(st.matrix.data.vst.0, gzfile(paste0(signaturePath.st.sample,"vst_free.tsv.gz")), quote=F, sep="\t")
#	}
#}
