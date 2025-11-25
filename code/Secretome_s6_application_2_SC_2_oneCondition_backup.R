source("Secretome_s0_path.R")

library(anndata)
library(SecAct)
library(Seurat)
library(ggplot2)
library(patchwork)

#SWARM -t 20 -g 200 --time 36:00:00
#SWARM -t 2 -g 200 --time 5:00:00
args = commandArgs(trailingOnly=TRUE)
cancer <- args[1]


if(TRUE)
{
	adata <- read_h5ad(paste0(applicationPath,"PanCancer/",cancer,".h5ad"))
		
	meta <- adata$obs
	meta$newCluster <- "Others"
	
	meta[meta$majorCluster=="CD4T","newCluster"] <- "CD4T"
	meta[meta$majorCluster=="CD8T","newCluster"] <- "CD8T"
	meta[meta$majorCluster=="Endothelial","newCluster"] <- "Endothelial"
	
	meta[meta$subCluster=="CD4T08_Treg_FOXP3","newCluster"] <- "Treg"
	meta[meta$subCluster=="CD8T07_Tex_CXCL13","newCluster"] <- "Tex"
	
	meta[sapply(strsplit(as.character(meta$subCluster),"_",fixed=T),function(x) return(x[1]))%in%c("B01","B02","B03","B04","B05","B06","B07","B013","B014"),"newCluster"] <- "B"
	meta[sapply(strsplit(as.character(meta$subCluster),"_",fixed=T),function(x) return(x[1]))%in%c("Epi"),"newCluster"] <- "Tumor"
	
	meta[sapply(strsplit(as.character(meta$subCluster),"_",fixed=T),function(x) return(x[2]))%in%c("CD16hiNK","CD16loNK"),"newCluster"] <- "NK"
	meta[sapply(strsplit(as.character(meta$subCluster),"_",fixed=T),function(x) return(x[2]))%in%c("cDC","cDC1","cDC2"),"newCluster"] <- "cDC"
	meta[sapply(strsplit(as.character(meta$subCluster),"_",fixed=T),function(x) return(x[2]))%in%c("pDC"),"newCluster"] <- "pDC"
	meta[sapply(strsplit(as.character(meta$subCluster),"_",fixed=T),function(x) return(x[2]))%in%c("Mph"),"newCluster"] <- "Macrophage"
	meta[sapply(strsplit(as.character(meta$subCluster),"_",fixed=T),function(x) return(x[2]))%in%c("Fb","iCAF","myCAF","apCAF"),"newCluster"] <- "Fibroblast"
	
	meta$newCluster <- factor(meta$newCluster)
	
	Seurat_obj <- CreateSeuratObject(counts = t(as.matrix(adata$X)), meta.data = meta)
	Seurat_obj <- subset(Seurat_obj, subset = newCluster != "Others")
	Seurat_obj <- subset(Seurat_obj, subset = tissue == "Tumor")
	
	Seurat_obj <- SecAct.CCC.scRNAseq(
	  Seurat_obj, 
	  cellType_meta="newCluster",
	  condition_meta=NULL
	)  
	
	saveRDS(Seurat_obj, file = paste0(applicationPath,"PanCancer/",cancer,"_single_condition_CCC_Seurat_obj.rds"))
}


if(FALSE)
{
	cancers <- list.files(paste0(applicationPath,"PanCancer/"))
	cancers <- cancers[grepl("CCC",cancers)]
	cancers <- gsub("_single_condition_CCC_Seurat_obj.rds","",cancers)
	
	smy <- data.frame()
	for(cancer in cancers)
	{
		Seurat_obj <- readRDS(paste0(applicationPath,"PanCancer/",cancer,"_single_condition_CCC_Seurat_obj.rds"))
		
		ccc <- Seurat_obj @misc $SecAct_output $SecretedProteinCCC
		
		SP_stat <- table(ccc[,2]) 
		unique_SP <- names(SP_stat)[SP_stat==1]
		
		ccc_unique <- ccc[ccc[,2]%in%unique_SP,]
		
		if(length(rownames(ccc_unique))==0) next
		
		smy[paste0(cancer,rownames(ccc_unique)),"cancer"] <- cancer
		smy[paste0(cancer,rownames(ccc_unique)),"ccc"] <- rownames(ccc_unique)
	}
	
	
	smy_stat <- table(smy[,"ccc"]) 
	smy_stat[grepl("SAA1",names(smy_stat))]
	
	
	smy_stat_consensus <- smy_stat[smy_stat>=2]
	
	fg.df <- as.data.frame(smy_stat_consensus)
	fg.df <- fg.df[order(fg.df[,2],decreasing=T),]
	fg.df[,1] <- factor(fg.df[,1], levels=fg.df[,1])
	
	p1 <- ggplot(fg.df,aes(x = Var1, y = Freq, fill=Freq))+
			geom_bar(stat = "identity", position=position_dodge())+
			ggtitle(" ")+
			xlab(" ")+
			ylab("Number of scRNA-Seq cohorts")+
			theme_bw()+ 
			theme(
			  panel.grid.major.x = element_blank(),
			  panel.background = element_blank(),
			  axis.text = element_text(size=11,colour = "black"),
			  axis.title = element_text(size=12,colour = "black"),
			  axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1 ),
			  axis.text.y = element_text(angle = 90, hjust=0.5),
			  legend.position="none"
			) 
	ggsave(paste0(applicationPath,"PanCancer/twoCellTypeOnly_stat.png"), p1, width = 26, height =15, dpi=200, units = "cm",limitsize = FALSE)
	write.csv(fg.df, paste0(applicationPath,"PanCancer/twoCellTypeOnly_stat.csv"), quote=FALSE)
	write.csv(smy, paste0(applicationPath,"PanCancer/twoCellTypeOnly.csv"), quote=FALSE)


	smy_unique <- sapply(strsplit(names(smy_stat_consensus),"_",fixed=T),function(x) return(x[2]))
	
	for(gene in smy_unique)
	{
		cat(gene)
		cat(names(smy_stat)[grepl(gene, names(smy_stat))])
		cat("\n")
	}
	
	
	for(SP in c("MDK","LGALS3"))
	{
		smy_SP <- smy[grepl(SP,smy[,"ccc"]),]
		cancers <- smy_SP[,1]
	
	
		for(cancer in cancers)
		{
			Seurat_obj <- readRDS(paste0(applicationPath,"PanCancer/",cancer,"_single_condition_CCC_Seurat_obj.rds"))
			
			act <- Seurat_obj @misc $SecAct_output $SecretedProteinActivity $zscore
			print(cancer)
			print(act[SP,])
			
			fg.df1 <- data.frame(cellType=colnames(act), value=act[SP,], group=act[SP,]>0)
			
			p1 <- ggplot(fg.df1,aes(x = cellType, y = value, fill=group))+
				geom_bar(stat = "identity", position=position_dodge())+
				scale_fill_manual(values=c("grey","brown"))+
				ggtitle(" ")+
				xlab(" ")+
				ylab("activity")+
				theme_bw()+ 
				theme(
				  panel.grid = element_blank(),
				  panel.background = element_blank(),
				  axis.text = element_text(size=10,colour = "black"),
				  axis.title = element_text(size=10,colour = "black"),
				  legend.position="none"
				) 
			
			fg.df2 <- data.frame()
			for(cellType in colnames(act))
			{
				fg.df2[cellType,"cellType"] <- cellType
				fg.df2[cellType,"value"] <- Seurat_obj @misc $SecAct_output $ SecretedProteinExpression[[cellType]] [SP,"exp_logFC"]
				fg.df2[cellType,"group"] <- Seurat_obj @misc $SecAct_output $ SecretedProteinExpression[[cellType]] [SP,"exp_logFC"]>0
			}
			
			p2 <- ggplot(fg.df2,aes(x = cellType, y = value, fill=group))+
				geom_bar(stat = "identity", position=position_dodge())+
				scale_fill_manual(values=c("grey","brown"))+
				ggtitle(paste0(cancer, " ------ ", smy_SP[smy_SP[,1]==cancer,2]))+
				xlab(" ")+
				ylab("expression (logFC)")+
				theme_bw()+ 
				theme(
				  panel.grid = element_blank(),
				  panel.background = element_blank(),
				  axis.text = element_text(size=10,colour = "black"),
				  axis.title = element_text(size=10,colour = "black"),
				  legend.position="none"
				) 
			
			ggsave(paste0(applicationPath,"PanCancer/",cancer,"_",SP,".png"), p2/p1, width = 22, height =15, dpi=200, units = "cm",limitsize = FALSE)
			
			
			source("ifun.R")
	
			v1 <- fg.df1[,2]
			v2 <- fg.df2[,2]
			names(v1) <- fg.df1[,1]
			names(v2) <- fg.df2[,1]
			
			v1 <- v1[sort(names(v1))]
			v2 <- v2[sort(names(v2))]
			
			v1[v1<0] <- 0
			v2[v2<0] <- 0
			
			
			p3 <- two_col_heatmap_with_arrow_legends(
			  v2, v1, 
			  left_title = paste0(SP,"\nExpression"), 
			  right_title = paste0(SP,"\nActivity"), 
			  mid_title = " ",
			  title = cancer, gap_x = 3
			)
			
			ggsave(paste0(applicationPath,"PanCancer/",cancer,"_",SP,"_net.png"), p3, width = 15, height =12, dpi=200, units = "cm",limitsize = FALSE)
		}
	}
	
}






### two conditional ccc start

if(FALSE)
{
	adata <- read_h5ad(paste0(applicationPath,"PanCancer/",cancer,".h5ad"))
		
	meta <- adata$obs
	meta$newCluster <- "Others"
	
	meta[meta$majorCluster=="CD4T","newCluster"] <- "CD4T"
	meta[meta$majorCluster=="CD8T","newCluster"] <- "CD8T"
	meta[meta$majorCluster=="Endothelial","newCluster"] <- "Endothelial"
	
	meta[meta$subCluster=="CD4T08_Treg_FOXP3","newCluster"] <- "Treg"
	meta[meta$subCluster=="CD8T07_Tex_CXCL13","newCluster"] <- "Tex"
	
	meta[sapply(strsplit(as.character(meta$subCluster),"_",fixed=T),function(x) return(x[1]))%in%c("B01","B02","B03","B04","B05","B06","B07","B013","B014"),"newCluster"] <- "B"
	meta[sapply(strsplit(as.character(meta$subCluster),"_",fixed=T),function(x) return(x[1]))%in%c("Epi"),"newCluster"] <- "Tumor"
	
	meta[sapply(strsplit(as.character(meta$subCluster),"_",fixed=T),function(x) return(x[2]))%in%c("CD16hiNK","CD16loNK"),"newCluster"] <- "NK"
	meta[sapply(strsplit(as.character(meta$subCluster),"_",fixed=T),function(x) return(x[2]))%in%c("cDC","cDC1","cDC2"),"newCluster"] <- "cDC"
	meta[sapply(strsplit(as.character(meta$subCluster),"_",fixed=T),function(x) return(x[2]))%in%c("pDC"),"newCluster"] <- "pDC"
	meta[sapply(strsplit(as.character(meta$subCluster),"_",fixed=T),function(x) return(x[2]))%in%c("Mph"),"newCluster"] <- "Macrophage"
	meta[sapply(strsplit(as.character(meta$subCluster),"_",fixed=T),function(x) return(x[2]))%in%c("Fb","iCAF","myCAF","apCAF"),"newCluster"] <- "Fibroblast"
	
	meta$newCluster <- factor(meta$newCluster)
	
	Seurat_obj <- CreateSeuratObject(counts = t(as.matrix(adata$X)), meta.data = meta)
	Seurat_obj <- subset(Seurat_obj, subset = newCluster != "Others")
	
	Seurat_obj <- SecAct.CCC.scRNAseq(
	  Seurat_obj, 
	  cellType_meta="newCluster",
	  condition_meta="tissue", 
	  conditionCase="Metastasis", 
	  conditionControl="Tumor"
	)  
	
	saveRDS(Seurat_obj, file = paste0(applicationPath,"PanCancer/",cancer,"_Seurat_obj.rds"))
}


if(FALSE)
{
	smy <- data.frame()
	smy2 <- data.frame()
	
	cancers <- c("OV_Zheng2023","LUAD_GSE131907","TNBC_GSE169246","PDAC_GSE205013","PTC_GSE184362")
	for(cancer in cancers)
	{
		if(cancer=="OV_Zheng2023")
		{
			Seurat_obj <- readRDS(paste0(applicationPath,"OV/OV_Metastatic_Seurat_obj.rds"))
		}else{
			Seurat_obj <- readRDS(paste0(applicationPath,"PanCancer/",cancer,"_Seurat_obj.rds"))
		}
		
		ccc <- Seurat_obj@misc$SecAct_output$SecretedProteinCCC
		ccc <- ccc[!is.na(ccc[,1]),]
		
		if(cancer=="OV_Zheng2023")
		{	
			ccc[ccc[,1]=="DC",1] <- "cDC"
			ccc[ccc[,3]=="DC",3] <- "cDC"
		}
		
		SP_stat <- table(ccc[,"secretedProtein"])
		SPs <- names(SP_stat)[SP_stat==1]
		
		smy[rownames(ccc)[ccc[,"secretedProtein"]%in%SPs],cancer] <- 1
		
		smy2[rownames(ccc),cancer] <- 1
	}
	smy2[is.na(smy2)] <- 0
	
	smy2_rowSum <- rowSums(smy2)
	
	smy[is.na(smy)] <- ""
	write.csv(smy, paste0(applicationPath,"PanCancer/twoCellTypesOnly.csv"), quote=F)
}

### two conditional ccc end



if(FALSE)
{
	adata <- read_h5ad(paste0(applicationPath,"PanCancer/",cancer,".h5ad"))
	cdata_T <- t(as.matrix(adata$X))
	
	
	rownames(cdata_T) <- transferSymbol(rownames(cdata_T))
	cdata_T <- rm_duplicates_sparse(cdata_T)
	
	cdata_T <- sweep(cdata_T, 2, Matrix::colSums(cdata_T), "/") *1e5
	cdata_T <- log2(cdata_T + 1)
	
	
	# option 1: normalize across all cell types
	cdata_T <- cdata_T-rowMeans(cdata_T, na.rm=T)
	
	
	ref <- paste0(signatureCombPath,"AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3_1k_vst_condition_logUMI_cellType_0.9.tsv")
	a <- 5e+05
		
	res <- SecAct.inference(Y=cdata_T, SigMat=ref, lambda=a, nrand=1000)
	save(res, file = paste0(applicationPath,"PanCancer/",cancer,".RData"))
}	


if(FALSE)
{
	adata <- read_h5ad(paste0(applicationPath,"PanCancer/",cancer,".h5ad"))
		
	meta <- adata$obs[,c("majorCluster","subCluster")]
	meta$newCluster <- "Others"
	
	meta[meta$majorCluster=="CD4T","newCluster"] <- "CD4T"
	meta[meta$majorCluster=="CD8T","newCluster"] <- "CD8T"
	meta[meta$majorCluster=="Endothelial","newCluster"] <- "Endothelial"
	
	meta[sapply(strsplit(as.character(meta$subCluster),"_",fixed=T),function(x) return(x[1]))%in%c("B01","B02","B03","B04","B05","B06","B07","B013","B014"),"newCluster"] <- "B"
	meta[sapply(strsplit(as.character(meta$subCluster),"_",fixed=T),function(x) return(x[1]))%in%c("Epi"),"newCluster"] <- "Malignant"
	
	meta[sapply(strsplit(as.character(meta$subCluster),"_",fixed=T),function(x) return(x[2]))%in%c("CD16hiNK","CD16loNK"),"newCluster"] <- "NK"
	meta[sapply(strsplit(as.character(meta$subCluster),"_",fixed=T),function(x) return(x[2]))%in%c("cDC","cDC1","cDC2"),"newCluster"] <- "cDC"
	meta[sapply(strsplit(as.character(meta$subCluster),"_",fixed=T),function(x) return(x[2]))%in%c("pDC"),"newCluster"] <- "pDC"
	meta[sapply(strsplit(as.character(meta$subCluster),"_",fixed=T),function(x) return(x[2]))%in%c("Mo"),"newCluster"] <- "Monocyte"
	meta[sapply(strsplit(as.character(meta$subCluster),"_",fixed=T),function(x) return(x[2]))%in%c("Mph"),"newCluster"] <- "Macrophage"
	meta[sapply(strsplit(as.character(meta$subCluster),"_",fixed=T),function(x) return(x[2]))%in%c("immNeu","mNeu"),"newCluster"] <- "Neutrophil"
	meta[sapply(strsplit(as.character(meta$subCluster),"_",fixed=T),function(x) return(x[2]))%in%c("Fb","iCAF","myCAF","apCAF"),"newCluster"] <- "Fibroblast"
	
	meta$newCluster <- factor(meta$newCluster)
	
	cellType_vec <- meta$newCluster
	cellType_stat <- table(cellType_vec)
	
			
	load(paste0(applicationPath,"PanCancer/",cancer,".RData"))
	sp <- as.matrix(res$zscore)
	sp <- expand_rows(sp)
	
	cdata_T <- t(as.matrix(adata$X))
	rownames(cdata_T) <- transferSymbol(rownames(cdata_T))
	cdata_T <- rm_duplicates_sparse(cdata_T)
	cdata_T <- sweep(cdata_T, 2, Matrix::colSums(cdata_T), "/") *1e5
	cdata_T <- log2(cdata_T + 1)
	cdata_T <- cdata_T[rownames(cdata_T)%in%rownames(sp),]
	
	
	avg_exp <- sapply(levels(cellType_vec), function(ct) {
	  rowMeans(cdata_T[, cellType_vec == ct, drop = FALSE])
	  })
	  
	avg_act <- sapply(levels(cellType_vec), function(ct) {
	  rowMeans(sp[, cellType_vec == ct, drop = FALSE])
	  })
	
	avg_exp <- avg_exp[,cellType_stat>=100,drop=F]
	avg_act <- avg_act[,cellType_stat>=100,drop=F]
	
	write.csv(avg_exp, paste0(applicationPath,"PanCancer/",cancer,"_cellType_exp.csv"), quote=F)
	write.csv(avg_act, paste0(applicationPath,"PanCancer/",cancer,"_cellType_act.csv"), quote=F)
}	




if(FALSE)
{
	genes <- c("TGFB1","IFNG","FLT3LG","AOAH","LY86")
	
	avg_act <- read.csv(paste0(applicationPath,"PanCancer/",cancer,"_cellType_act.csv"), header=TRUE, row.names=1)
	genes <- rownames(avg_act)
	
	smy <- data.frame()
	for(gene in genes)
	{	
		smy_exp <- data.frame()
		smy_act <- data.frame()
		
		fileNames <- list.files(paste0(applicationPath,"PanCancer/"))
		fileNames <- fileNames[grepl("_cellType_act.csv",fileNames)]
		cancers <- gsub("_cellType_act.csv","",fileNames)
	
		for(cancer in cancers)
		{
			avg_exp <- read.csv(paste0(applicationPath,"PanCancer/",cancer,"_cellType_exp.csv"), header=TRUE, row.names=1)
			avg_act <- read.csv(paste0(applicationPath,"PanCancer/",cancer,"_cellType_act.csv"), header=TRUE, row.names=1)
			
			avg_exp <- avg_exp[,colnames(avg_exp)!="Others",drop=F]
			avg_act <- avg_act[,colnames(avg_act)!="Others",drop=F]
			
			avg_exp <- t(scale(t(avg_exp)))
			
			if(!gene%in%rownames(avg_exp)) next
			smy_exp[colnames(avg_exp),cancer] <- unlist(avg_exp[gene,])
			smy_act[colnames(avg_act),cancer] <- unlist(avg_act[gene,])
		}
		if(ncol(smy_exp)==0) next
		
		smy_exp <- smy_exp[c("Malignant","Fibroblast","Endothelial","B","CD4T","CD8T","NK","cDC","pDC","Monocyte","Macrophage","Neutrophil"),]
		smy_act <- smy_act[c("Malignant","Fibroblast","Endothelial","B","CD4T","CD8T","NK","cDC","pDC","Monocyte","Macrophage","Neutrophil"),]
		
		smy_exp_median <- apply(smy_exp,1,function(x) median(x, na.rm=T))
		smy_act_median <- apply(smy_act,1,function(x) median(x, na.rm=T))
		
		names(smy_exp_median) <- paste0(names(smy_exp_median),"_exp")
		names(smy_act_median) <- paste0(names(smy_act_median),"_act")
		
		smy_exp_act_median <- c(smy_exp_median, smy_act_median)
		smy[names(smy_exp_act_median),gene] <- smy_exp_act_median
		
		#library(ComplexHeatmap)
		#
		#mat <- t(as.matrix(smy_exp))
		#png(paste0(applicationPath,"PanCancer/",gene,"_exp.png"), width = 18, height = 20, res=200, units = "cm")
		#
		#column_ha <- columnAnnotation(
		#	Expr = anno_boxplot(mat, height = unit(3, "cm"), annotation_name_side = "left")
		#)
		#		
		#ht <- Heatmap(mat,
		#	name = "expr",
		#	column_title = paste0(gene,"_expr"),
		#	rect_gp = gpar(col = "white", lwd = 2),
		#	col = circlize::colorRamp2(c(-1, 0,1), c("#91bfdb", "white", "#fc8d59")),
		#	row_names_max_width = max_text_width(
		#        rownames(mat), 
		#        gp = gpar(fontsize = 12)
		#        ),
		#    column_names_max_height = max_text_width(
		#        colnames(mat), 
		#        gp = gpar(fontsize = 12)
		#        ),
		#    column_names_rot = 38,
		#    top_annotation = column_ha,
		# 	cluster_rows = FALSE,
		#	cluster_columns = FALSE
		#)
		#draw(ht)
		#
		#dev.off()
		#
		#
		#mat <- t(as.matrix(smy_act))
		#png(paste0(applicationPath,"PanCancer/",gene,"_act.png"), width = 18, height = 20, res=200, units = "cm")
		#
		#column_ha <- columnAnnotation(
		#	Act = anno_boxplot(mat, height = unit(3, "cm"), annotation_name_side = "left")
		#)
		#		
		#ht <- Heatmap(mat,
		#	name = "act",
		#	column_title = paste0(gene,"_act"),
		#	rect_gp = gpar(col = "white", lwd = 2),
		#	col = circlize::colorRamp2(c(-10, 0,10), c("#91bfdb", "white", "#fc8d59")),
		#	row_names_max_width = max_text_width(
		#        rownames(mat), 
		#        gp = gpar(fontsize = 12)
		#        ),
		#    column_names_max_height = max_text_width(
		#        colnames(mat), 
		#        gp = gpar(fontsize = 12)
		#        ),
		#    column_names_rot = 38,
		#    top_annotation = column_ha,
		# 	cluster_rows = FALSE,
		#	cluster_columns = FALSE
		#)
		#draw(ht)
		#
		#dev.off()
			
	}
	
	smy <- t(smy)
	smy <- round(smy,2)
	write.csv(smy, paste0(applicationPath,"PanCancer/smy.csv"), quote=F)
	
	
	smy <- as.matrix(read.csv(paste0(applicationPath,"PanCancer/smy.csv"), as.is=T, row.names=1, header=T))
	smy <- smy[,c(-10,-12,-22,-24)]
	
	ccc <- data.frame()
	for(gene in rownames(smy))
	{
		expr <- smy[gene,1:10]
		act  <- smy[gene,11:20]
		
		expr_sorted <- sort(expr, decreasing=T)
		act_sorted <- sort(act, decreasing=T)
		
		expr_sd <- sd(expr_sorted[2:10])
		act_sd  <- sd(act_sorted[2:10])
		
		
		expr_z <- (expr_sorted[1]-expr_sorted[2]) / expr_sd
		act_z <- (act_sorted[1]-act_sorted[2]) / act_sd

		if(expr_z>2&act_z>2)
		{
			ccc[gene,"expr"] <- names(expr_sorted)[1]
			ccc[gene,"gene"] <- gene
			ccc[gene,"act"] <- names(act_sorted)[1]
		}
	}
	ccc[,1] <- sapply(strsplit(ccc[,1],"_",fixed=T),function(x) return(x[1]))
	ccc[,3] <- sapply(strsplit(ccc[,3],"_",fixed=T),function(x) return(x[1]))
	ccc <- ccc[ccc[,1]!=ccc[,3],]
	ccc[order(ccc[,3]),]
	
}
