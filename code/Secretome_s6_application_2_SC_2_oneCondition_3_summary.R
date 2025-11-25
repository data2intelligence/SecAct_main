source("Secretome_s0_path.R")

library(anndata)
library(SecAct)
library(Seurat)
library(ggplot2)
library(patchwork)

outputPath <- paste0(applicationPath,"scRNAseq_PanCancer/")

cancers <- list.files(outputPath)
cancers <- cancers[grepl("CCC",cancers)]
cancers <- gsub("_single_condition_CCC_Seurat_obj.rds","",cancers)

smy <- data.frame()
for(cancer in cancers)
{
	Seurat_obj <- readRDS(paste0(outputPath,cancer,"_single_condition_CCC_Seurat_obj.rds"))
	
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
ggsave(paste0(outputPath,"twoCellTypeOnly_stat.png"), p1, width = 26, height =15, dpi=200, units = "cm",limitsize = FALSE)
write.csv(fg.df, paste0(outputPath,"twoCellTypeOnly_stat.csv"), quote=FALSE)
write.csv(smy, paste0(outputPath,"twoCellTypeOnly.csv"), quote=FALSE)


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
		Seurat_obj <- readRDS(paste0(outputPath,cancer,"_single_condition_CCC_Seurat_obj.rds"))
		
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
		
		ggsave(paste0(outputPath,cancer,"_",SP,".png"), p2/p1, width = 22, height =15, dpi=200, units = "cm",limitsize = FALSE)
		
		
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
		
		ggsave(paste0(outputPath,cancer,"_",SP,"_net.png"), p3, width = 15, height =12, dpi=200, units = "cm",limitsize = FALSE)
	}
}


