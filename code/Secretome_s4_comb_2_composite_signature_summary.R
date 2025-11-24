source("Secretome_s0_path.R")

# five parts
# 1 count the number of signatures for the same secreted proteins
# 2 compare pan-cancer and cancer type-specific signatures
# 3 compare cancer type-specific signatures
# 4 compare cancer type-specific and non-cancer type-specific signatures
# 5 compare signature count from other methods, e.g., CytoSig, NicheNet, ImmuneDic


###########
# fig 1i #
###########

sigs <- read.table(paste0(signatureCombPath,"AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3_vst.tsv"),sep="\t",check.names=F)
stat <- read.csv(paste0(signatureCombPath,"AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_vst_stat.csv"),as.is=T,row.names=1,header=T)

stat <- stat[colnames(sigs),,drop=F]

library(ggplot2)
p0 <- ggplot(stat,aes(x=sigCount)) + 
	geom_histogram(fill="#74c7b9", colour="white", alpha=0.8, position = "stack")+ #fill="#fc8d59"
	ylab("Number of\nSecreted Proteins")+
	xlab("Number of Signatures")+
	theme_classic()+ 
	theme(
		panel.background = element_blank(),
		panel.grid = element_blank(),
		axis.title = element_text(colour = "black"),
		axis.text = element_text(colour = "black"),
		axis.text.x = element_text(colour = "black",hjust=0.8),
		legend.position="none"
	)
ggsave(paste0(signatureCombPath,"AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3_vst_stat.png"), p0, width = 6.6, height = 5.5, dpi=400, units = "cm", limitsize =FALSE)



###########
# fig s3a #
###########

topcancers <- c("BRCA","PDAC","GBM","CRC","LUAD","KIRC","OV","LIHC","PRAD")
names(topcancers) <- c("Breast","Pancreatic","Glioblastoma","Colorectal","Lung-Adeno","Kidney","Ovarian","Liver","Prostate")

cancersFromPRECOG <- c("PRAD","PDAC","LUAD")

for(topcancer in topcancers)
{
	if(!topcancer%in%cancersFromPRECOG)
	{
		dataPath <- paste0(GEOPath,topcancer,"/")
		cancers <- list.files(dataPath)
		cancers <- cancers[grepl("subtract",cancers)]
		cancers <- cancers[!grepl("MFS",cancers)]
		names(cancers) <- sapply(strsplit(cancers,".",fixed=T),function(x) return(x[1]))
	}else{
		dataPath <- paste0(GEOPath,"PRECOG/")
		cancers <- list.files(dataPath)
		cancers <- cancers[grepl(names(topcancers)[topcancers==topcancer],cancers)]
		cancers <- cancers[grepl("expression",cancers)]
		cancers <- cancers[grepl("GSE",cancers)]
		names(cancers) <- sapply(strsplit(cancers,"_",fixed=T),function(x) return(x[2]))		
		cancers <- cancers[nchar(names(cancers))==8]
	}
		
	sigs <- c(
		"AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3_vst",
		paste0("AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_",topcancer,"_filterByPan_ds3_vst"),
		paste0("AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_non-",topcancer,"_filterByPan_ds3_vst")
	)
	names(sigs) <- c("Pan-cancer","Cancer-type-specific","Leave-one-cancer-type-out")
	

	comb <- data.frame()
	for(cancer in cancers)
	{
		if(!topcancer%in%cancersFromPRECOG)
		{
			cdata_T_minusBG <- as.matrix(read.csv(paste0(dataPath,cancer),sep="\t",row.names=1))
		}else{
			cdata_T <- as.matrix(read.csv(paste0(dataPath,cancer),sep="\t",row.names=1))
			#print(cancer)
			
			if(grepl("EntrezCDF",cancer)|grepl("GSE11117_GPL6650",cancer)|grepl("GSE21034_GPL10264",cancer)) cdata_T <- log2(cdata_T+1)
			#print(cdata_T[1:6,1:6])

			cdata_T_minusBG <- cdata_T-rowMeans(cdata_T)
		}
		
		if(topcancer=="BRCA")
		{
			id2symbol <- read.csv(gzfile(paste0(NCBIPath,"NCBI_20251008_gene_result_id2symbol.csv.gz")))
			cdata_T_minusBG <- cdata_T_minusBG[rownames(cdata_T_minusBG)%in%id2symbol[,1],]
			rownames(cdata_T_minusBG) <- id2symbol[match(rownames(cdata_T_minusBG),id2symbol[,1]),2]
		}
		
		rownames(cdata_T_minusBG) <- transferSymbol(rownames(cdata_T_minusBG))
		cdata_T_minusBG <- rm_duplicates(cdata_T_minusBG)
		
		comb.sameSPs <- data.frame()
		for(i in 1:length(sigs))
		{
			sigmat <- as.matrix(read.table(paste0(signatureCombPath,sigs[i],".tsv"),sep="\t",check.names=F))
			smyMat <- cor.act.exp(X=sigmat, Y=cdata_T_minusBG)
			smyMat <- as.matrix(smyMat,ncol=1)
			
			smyMat <- smyMat[!is.na(smyMat[,1]),,drop=F]
			
			comb.sameSPs[paste0(cancer,sigs[i],rownames(smyMat)),"cancer"] <- names(cancers)[cancers%in%cancer]
			comb.sameSPs[paste0(cancer,sigs[i],rownames(smyMat)),"sig"] <- names(sigs)[i]
			comb.sameSPs[paste0(cancer,sigs[i],rownames(smyMat)),"gene"] <- rownames(smyMat)
			comb.sameSPs[paste0(cancer,sigs[i],rownames(smyMat)),"value"] <- smyMat[,1]
		}
		
		genes.df <- as.data.frame(table(comb.sameSPs[,"gene"]))
		genes <- genes.df[genes.df[,2]==length(unique(comb.sameSPs[,"sig"])),1]
		comb.sameSPs <- comb.sameSPs[comb.sameSPs[,"gene"]%in%genes,]
		comb <- rbind(comb,comb.sameSPs)
	}
	
	comb[,"sig"] <- factor(comb[,"sig"], levels=names(sigs))
	
	library(ggplot2)
	library(ggsignif)
	p2 <- ggplot(comb,aes(x = sig, y = value))+
		geom_hline(yintercept=0, color = "grey", linewidth=0.6, linetype="dashed")+
		geom_violin(aes(group=sig, fill=sig), width=0.7, trim=FALSE, alpha=0.7)+
		geom_boxplot(width=0.2, outlier.shape=NA)+
		geom_signif( comparisons = list(c("Cancer-type-specific","Pan-cancer"),c("Cancer-type-specific","Leave-one-cancer-type-out")) ,y_position = c(0.75, 0.82), test = "wilcox.test")+
		scale_fill_manual(values=c("#9163b6","#e0598b","#e2975d"))+
		ggtitle(names(topcancers)[topcancers==topcancer])+
		ylim(-0.5,1)+
		ylab("Pearson r value")+
		labs(fill="Composite signature")+
		theme_bw()+ 
		theme(
		  panel.grid = element_blank(),
		  panel.background = element_blank(),
		  plot.title = element_text(hjust = 0.5, size=18),
		  axis.text = element_text(size=15,colour = "black"),
		  axis.title = element_text(size=15,colour = "black"),
		  axis.text.x = element_blank(),
		  axis.ticks.x = element_blank(),
		  axis.title.x = element_blank(),
		  legend.position="bottom",
		  legend.text=element_text(size=15),
		  legend.title=element_text(size=15),
		  strip.text=element_text(size=13)
		) + facet_wrap(~ cancer, nrow =1) #scales = "free"
		
	ggsave(paste0(signatureCombPath,"cancer_specific_",topcancer,".png"), p2, width = 5*length(cancers) + 1, height = 8.2, dpi=200, units = "cm", limitsize = FALSE)
}



###########
# fig s3c #
###########

median_vec <- c()

smy <- data.frame()
for(topcancer in topcancers)
{
	sig1 <- read.table(paste0(signatureCombPath, "AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_",topcancer,"_filterByPan_ds3_vst.tsv"))
	sig2 <- read.table(paste0(signatureCombPath, "AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_non-",topcancer,"_filterByPan_ds3_vst.tsv"))

	olp_col <- intersect(colnames(sig1), colnames(sig2))
	olp_row <- intersect(rownames(sig1), rownames(sig2))
	
	sig1 <- sig1[olp_row,olp_col]
	sig2 <- sig2[olp_row,olp_col]
	
	cor_vec <- diag(cor(sig1, sig2))
	median_vec <- c(median_vec, median(cor_vec))
	
	smy[paste0(topcancer,"_",names(cor_vec)),"cancer1"] <- names(topcancers)[topcancers==topcancer]
	smy[paste0(topcancer,"_",names(cor_vec)),"gene"] <- names(cor_vec)
	smy[paste0(topcancer,"_",names(cor_vec)),"cor"] <- cor_vec
}


txt <- c(
	median=round(median(median_vec),3), 
	sd=round(sd(median_vec),3)
	)
txt
write.csv(txt, paste0(signatureCombPath,"compare_cancer_non-cancer_sig.txt"), quote=F)


smy <- cbind(smy, x="x")


library(ggplot2)
p2 <- ggplot(smy,aes(x = x, y = cor))+
	geom_hline(yintercept=0, color = "grey", linewidth=0.6, linetype="dashed")+
	geom_violin(aes(group=x), fill = "#b0d992", trim=FALSE)+
	geom_boxplot(width=0.1, outlier.shape=NA)+
	ggtitle(" ")+
	ylab("Pearson r value")+
	theme_bw()+
	theme(
	  panel.grid = element_blank(),
	  panel.background = element_blank(),
	  plot.title = element_text(hjust = 0.5, size=18),
	  axis.text = element_text(size=15,colour = "black"),
	  axis.title = element_text(size=15,colour = "black"),
	  axis.text.x = element_blank(),
	  axis.ticks.x = element_blank(),
	  axis.title.x = element_blank(),
	  legend.position="right",
	  legend.text=element_text(size=15),
	  legend.title=element_text(size=15),
	  strip.background = element_rect(fill = "white"),
	  strip.text = element_text(size = 15)
	) + facet_wrap(~cancer1, ncol =3) 
ggsave(paste0(signatureCombPath,"compare_cancer_non-cancer_sig.png"), p2, width = 18, height = 18, dpi=200, units = "cm", limitsize = FALSE)
write.csv(smy,paste0(signatureCombPath, "compare_cancer_non-cancer_sig.csv"),quote=F)



###########
# fig s3b #
###########

topcancers <- c('Pan-cancer'="Pan-cancer",topcancers)

sigList <- list()
for(topcancer in topcancers)
{	
	if(topcancer=="Pan-cancer")
	{
		sigList[[topcancer]] <- read.table(paste0(signatureCombPath, "AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3_vst.tsv"))
	}else{
		sigList[[topcancer]] <- read.table(paste0(signatureCombPath, "AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_",topcancer,"_filterByPan_ds3_vst.tsv"))
	}
}

median_vec <- c()

smy <- data.frame()
for(i in 1:length(topcancers))
{
	for(j in 1:length(topcancers))
	{
		if(i>=j) next
		
		sig1 <- sigList[[i]]
		sig2 <- sigList[[j]]
		
		olp_col <- intersect(colnames(sig1), colnames(sig2))
		olp_row <- intersect(rownames(sig1), rownames(sig2))
		
		sig1 <- sig1[olp_row,olp_col]
		sig2 <- sig2[olp_row,olp_col]
		
		cor_vec <- diag(cor(sig1, sig2))
		median_vec <- c(median_vec, median(cor_vec))

		smy[paste0(i,"_",j,"_",names(cor_vec)),"cancer1"] <- names(topcancers)[topcancers==topcancers[i]]
		smy[paste0(i,"_",j,"_",names(cor_vec)),"cancer2"] <- names(topcancers)[topcancers==topcancers[j]]
		smy[paste0(i,"_",j,"_",names(cor_vec)),"gene"] <- names(cor_vec)
		smy[paste0(i,"_",j,"_",names(cor_vec)),"cor"] <- cor_vec
	}
}

txt <- c(
	median=round(median(median_vec),3), 
	sd=round(sd(median_vec),3)
	)
txt
write.csv(txt, paste0(signatureCombPath,"compare_cancer-specific_sig.txt"), quote=F)

smy <- cbind(smy, x="x")
smy[[1]] <- factor(smy[[1]], levels=names(topcancers))
smy[[2]] <- factor(smy[[2]], levels=names(topcancers))

library(ggplot2)
gg <- ggplot(smy,aes(x = x, y = cor))+
	geom_hline(yintercept=0, color = "grey", linewidth=0.6, linetype="dashed")+
	geom_violin(aes(group=x), fill = "#7ac7e2", trim=FALSE)+
	geom_boxplot(width=0.1, outlier.shape=NA)+
	ylab("Pearson r value")+
	facet_grid(cancer2 ~ cancer1, drop = TRUE)+
	theme_bw()+
	theme(
	  panel.grid = element_blank(),
	  panel.background = element_blank(),
	  axis.text = element_text(size=15,colour = "black"),
	  axis.title = element_text(size=15,colour = "black"),
	  axis.text.x = element_blank(),
	  axis.ticks.x = element_blank(),
	  axis.title.x = element_blank(),
	  axis.title.y = element_text(hjust = 0.99),
	  legend.position="right",
	  legend.text=element_text(size=15),
	  legend.title=element_text(size=15),
	  strip.background = element_rect(fill = "white"),
	  strip.text = element_text(size = 18),
	  strip.text.y.right = element_text(angle = 0)
	)

#ggsave(paste0(signatureCombPath,"compare_cancer-specific_sig.png"), gg, width = 55, height = 55, dpi=200, units = "cm", limitsize = FALSE)
write.csv(smy,paste0(signatureCombPath, "compare_cancer-specific_sig.csv"),quote=F)


library(grid);
library(gtable);

grob <- ggplotGrob(gg);

# Remove facets
idx <- which(grob$layout$name %in% c(
	"panel-2-1", 
	"panel-3-1", "panel-3-2", 
	"panel-4-1", "panel-4-2", "panel-4-3", 
	"panel-5-1", "panel-5-2", "panel-5-3", "panel-5-4", 
	"panel-6-1", "panel-6-2", "panel-6-3", "panel-6-4", "panel-6-5", 
	"panel-7-1", "panel-7-2", "panel-7-3", "panel-7-4", "panel-7-5", "panel-7-6",
	"panel-8-1", "panel-8-2", "panel-8-3", "panel-8-4", "panel-8-5", "panel-8-6", "panel-8-7",
	"panel-9-1", "panel-9-2", "panel-9-3", "panel-9-4", "panel-9-5", "panel-9-6", "panel-9-7", "panel-9-8"));
for (i in idx) grob$grobs[[i]] <- nullGrob();

# Move x axes up
# axis-b-1 needs to move up 4 rows
# axis-b-2 needs to move up 2 rows
# idx <- which(grob$layout$name %in% c("axis-b-1", "axis-b-2"));
# grob$layout[idx, c("t", "b")] <- grob$layout[idx, c("t", "b")] - c(4, 2);
#
# Move y axes right
# axis-l-2 needs to move 2 columns to the right
# axis-l-3 needs ot move 4 columns to the right
idx <- which(grob$layout$name %in% c("axis-l-2", "axis-l-3", "axis-l-4", "axis-l-5", "axis-l-6", "axis-l-7", "axis-l-8", "axis-l-9"));
grob$layout[idx, c("l", "r")] <- grob$layout[idx, c("l", "r")] + c(2, 4, 6, 8, 10, 12, 14, 16);

png(paste0(signatureCombPath,"compare_cancer-specific_sig.png"), width = 48, height = 41, res=200, units = "cm")
grid.newpage();
grid.draw(grob);
dev.off()



###########
# fig 3a #
###########
f1 <- read.table(paste0(finalSignaturesPath,"SecAct.tsv.gz"))
f2 <- read.table(paste0(finalSignaturesPath,"CytoSig.tsv.gz"))
f3 <- read.table(paste0(finalSignaturesPath,"NicheNet.v1.tsv.gz"))
f4 <- read.table(paste0(finalSignaturesPath,"NicheNet.v2.tsv.gz"))
f5 <- read.table(paste0(finalSignaturesPath,"ImmuneDic.tsv.gz"))

f1_unique <- length(setdiff(colnames(f1), c(colnames(f2), colnames(f3), colnames(f4), colnames(f5))))
f2_unique <- length(setdiff(colnames(f2), c(colnames(f1), colnames(f3), colnames(f4), colnames(f5))))
f3_unique <- length(setdiff(colnames(f3), c(colnames(f1), colnames(f2), colnames(f5))))
f4_unique <- length(setdiff(colnames(f4), c(colnames(f1), colnames(f2), colnames(f5))))
f5_unique <- length(setdiff(colnames(f5), c(colnames(f1), colnames(f2), colnames(f3), colnames(f4))))

f1_shared <- ncol(f1)-f1_unique
f2_shared <- ncol(f2)-f2_unique
f3_shared <- ncol(f3)-f3_unique
f4_shared <- ncol(f4)-f4_unique
f5_shared <- ncol(f5)-f5_unique


fg.df <- data.frame(
	sig = rep(c("SecAct","CytoSig","NicheNet.v1","NicheNet.v2","ImmuneDic"),each=2),
	group = rep(c("Unique","Shared"),5),
	count = c(f1_unique,f1_shared, f2_unique,f2_shared, f3_unique,f3_shared, f4_unique,f4_shared, f5_unique,f5_shared)
)
fg.df[,1] <- factor(fg.df[,1],levels=c("SecAct","NicheNet.v2","NicheNet.v1","ImmuneDic","CytoSig"))

library(ggplot2)
p1 <- ggplot(fg.df, aes(x=sig, y=count, fill=group)) +
  geom_bar(stat="identity", color="white", width=1, alpha=0.55) +
  scale_fill_manual( values=c("#379E38","#EF7F13") )+
  scale_y_continuous(breaks = seq(0, 1250, by = 250))+
  ylab("# Secreted Proteins")+
	xlab(" ")+
	theme_classic()+ 
	theme(
		panel.background = element_blank(),
		panel.grid = element_blank(),
		axis.title = element_text(colour = "black"),
		axis.text = element_text(colour = "black"),
		axis.text.x = element_text(angle = 38,hjust = 1),
		legend.position = c(0.75, 0.85),
		legend.title=element_blank()
	)
ggsave(paste0(signatureCombPath,"compare_signature_unqiue_shared.png"), p1, width = 5.5, height = 5.5*1.33, dpi=500, units = "cm", limitsize =FALSE)
write.csv(fg.df,paste0(signatureCombPath, "compare_signature_unqiue_shared.csv"),quote=F)

