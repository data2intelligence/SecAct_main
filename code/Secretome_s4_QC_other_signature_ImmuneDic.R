source("Secretome_s0_path.R")
source("../../SpaCE/code/ifun.R")


#SWARM -t 2 -g 200 --time 30:00:00
args = commandArgs(trailingOnly=TRUE)
cancer <- args[1]


# immune-dictionary
# cancers <- c("B_cell","cDC1","cDC2","eTAC","ILC","Macrophage","MigDC","Monocyte","Neutrophil","NK_cell","pDC","T_cell_CD4","T_cell_CD8","T_cell_gd","Treg")
# http://www.immune-dictionary.org/app/home
# It is so wired to download the data
# 1. directly clicking the link works well
# 2. right click and open a new tab does not work
# the latter one, file can not readRDS

rds <- readRDS(paste0("/data/rub2/data/immune_dictionary/ref_data_",cancer,".RDS"))

meta <- data.frame(
	a=rds@meta.data$celltype,
	b=rds@meta.data$sample,
	c=rds@meta.data$channel_hashtag
)

new_symbol <- as.matrix(read.csv("../../Secretome/data/new_symbol_immunedictionary.txt",as.is=T,sep="\t"))
meta[meta[,"b"]%in%new_symbol[,1],"b"] <- new_symbol[match(meta[meta[,"b"]%in%new_symbol[,1],"b"],new_symbol[,1]),3]


cdata <- as.matrix(rds@assays$RNA@counts)
colnames(cdata) <- meta[,"b"]

mgi_ann = read.table("/data/rub2/project/Secretome/results/z_LIHC/HOM_MouseHumanSequence.rpt.gz", sep = "\t", header = TRUE)
mgi_ann = mgi_ann[, c(1,2,4)]
colnames(mgi_ann) = c("homoloid", "org", "symbol")
mgi_ann$org = gsub(", laboratory", "", mgi_ann$org)

mouseGenes <- rownames(cdata)
for(i in 1:length(mouseGenes))
{
	id <- mgi_ann[mgi_ann[,3]==mouseGenes[i],1]
	if(length(id)==0)
	{
		mouseGenes[i] <- "humanGeneNotExist"
	}else{
		hGenes <- mgi_ann[mgi_ann[,1]%in%id&mgi_ann[,2]=="human",3]
		if(length(hGenes)==0)
		{
			mouseGenes[i] <- "humanGeneNotExist"
		}else{
			if(length(hGenes)==1)
			{
				mouseGenes[i] <- hGenes
			}else{
				if(toupper(mouseGenes[i])%in%hGenes)
				{
					mouseGenes[i] <- toupper(mouseGenes[i])
				}else{
					mouseGenes[i] <- "humanGeneMultiple"
				}
			}
		}
	}
}
mouseGenes -> rownames(cdata)
cdata <- cdata[!rownames(cdata)%in%c("humanGeneNotExist","humanGeneMultiple"),]

rownames(cdata) <- transferSymbol(rownames(cdata))
cdata <- rm_duplicates(cdata)

cdata <- t(t(cdata)*1e5/colSums(cdata))
cdata <- log2(cdata+1)


smy <- data.frame()
stat <- data.frame()
for(gene in setdiff(unique(colnames(cdata)),"PBS"))
{
	smy[rownames(cdata),gene]<- 
		rowMeans(cdata[,colnames(cdata)==gene,drop=F]) - 
		rowMeans(cdata[,colnames(cdata)=="PBS",drop=F])
	stat["count",gene] <- sum(colnames(cdata)==gene)
}
write.table(smy, paste0(validationPath,"singlecell3/",cancer,".sig"),sep="\t",quote=F)
write.table(stat, paste0(validationPath,"singlecell3/",cancer,".stat"),sep="\t",quote=F)


celltype <- cancer

# TCGA
cancers <- list.files(paste0("/data/rub2/data/TCGA_firebrowse/mRNAseq_gene_exp_symbol"))
cancers <- setdiff(cancers, c("DLBC","LAML","THYM"))


qc <- data.frame()
for(cancer in cancers)
{
	for(ppl in c("Xena"))
	{
		if(ppl=="Xena")
		{
			cdata <- read.Xena(cancer)	
		}else{
			cdata <- read.CIDE(cancer)	
		}
		
		for(geneScale in c("filter"))
		{
			if(geneScale=="filter")
			{
				cdata <- filter.counts(cdata)
			}
			
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
				
				smyMat <- cor.act.exp(X=smy, Y=cdata_T_minusBG)
				
				qc[rownames(smyMat),cancer] <- smyMat[,1]

				write.table(qc, paste0(validationPath,"singlecell3/",celltype,".qc"),sep="\t",quote=F)

			}
		}
	}
}




if(FALSE)
{


cancers <- c("B_cell","cDC1","cDC2","eTAC","ILC",
"Macrophage","MigDC","Monocyte","Neutrophil","NK_cell",
"pDC","T_cell_CD4","T_cell_CD8","T_cell_gd","Treg")

sigMean <- function(sigmat)
{
	sigmat <- sigmat[rowSums(sigmat==0)<ncol(sigmat)*0.9,,drop=F]
	apply(sigmat,1,function(x) median(x,na.rm=T))
}

comb <- data.frame()
stat <- data.frame()
qcs <- data.frame()
for(cancer in cancers)
{
	smy <- read.table(paste0(validationPath,"singlecell3/",cancer,".sig"),sep="\t")
	stt <- read.table(paste0(validationPath,"singlecell3/",cancer,".stat"),sep="\t")
	qc <- read.table(paste0(validationPath,"singlecell3/",cancer,".qc"),sep="\t")
	
	#n<- which(cancers==cancer) + 20
	#write.table(smy,paste0("/data/rub2/project/Secretome/code/CytoSig-master/CytoSig/signature.centroid.",n),quote=F,sep="\t")
	
	for(gene in colnames(smy))
	{
		comb[rownames(smy),paste0(cancer,"@",gene)] <- smy[,gene]
	}
	
	stat[names(stt[1,]),cancer] <- stt[1,]
	
	qcmean <- apply(qc, 1, function(x) mean(x,na.rm=T))
	qcs[names(qcmean),cancer] <- qcmean
}
comb <- as.matrix(comb)

genes <- sapply(strsplit(colnames(comb),"@",fixed=T),function(x) return(paste0(x[2])))

combcomb <- data.frame()
for(gene in unique(genes)) 
{
	combsub <- comb[,genes%in%gene]
	colnames(combsub) <- sapply(strsplit(colnames(combsub),"@",fixed=T),function(x) return(paste0(x[1])))
	
	if(gene%in%rownames(qcs))
	{
		celltypes <- colnames(qcs)[qcs[gene,]>0.1]
		if(length(celltypes) < 3)
		{
			temp <- unlist(qcs[gene,])
			temps <- sort(temp,decreasing=T)
			celltypes <- names(temps)[1:3]
		}
	}else{
		celltypes <- colnames(qcs)
	}
	
	temp <- sigMean(combsub[,celltypes])
	combcomb[names(temp),gene] <- temp
}

combcomb <- combcomb[apply(combcomb,1,function(x) sum(is.na(x)))==0,]
write.table(combcomb,paste0(signatureCombPath,"ImmuneDicQC.tsv"),quote=F,sep="\t")





















cancers <- c("B_cell","cDC1","cDC2","eTAC","ILC","Macrophage","MigDC","Monocyte","Neutrophil","NK_cell","pDC","T_cell_CD4","T_cell_CD8","T_cell_gd","Treg")

smy <- data.frame()
for(cancer in cancers)
{
	act <- as.matrix(read.table(paste0(validationPath,"singlecell3/",cancer,".Zscore"),sep="\t",check.names=F)	)
	colnames(act) <- sapply(strsplit(colnames(act),".",fixed=T),function(x) return(x[1]))
	
	act.m <- reshape2::melt(act)
	act.m <- act.m[as.character(act.m[,1])==as.character(act.m[,2]),]
	
	smy[paste0(cancer,1:nrow(act.m)),"cellType"] <- cancer
	smy[paste0(cancer,1:nrow(act.m)),"gene"] <- act.m[,1]
	smy[paste0(cancer,1:nrow(act.m)),"value"] <- act.m[,3]
}

smy <- smy[smy[,2]!="ADIPOQ",]

library(ggplot2)
p1 <- ggplot(smy,aes(x = cellType, y = value))+
	geom_hline(yintercept=0, color = "grey1", linewidth=0.8, linetype = 'dotted')+
	geom_violin(aes(group=cellType, fill=cellType),trim=FALSE)+
	geom_boxplot(width=0.15,outlier.shape=NA)+
	ggtitle("Immune Dictionary (single cell cytokine treatment)")+
	xlab(" ")+
	ylab("Act Diff")+
	theme_bw()+ 
	theme(
	  panel.grid = element_blank(),
	  panel.background = element_blank(),
	  plot.title = element_text(hjust = 0.5),
	  axis.text = element_text(size=10,colour = "black"),
	  axis.title = element_text(size=10,colour = "black"),
	  legend.position="none"
	)+facet_wrap(~gene, ncol = 3)

ggsave(paste0(validationPath,"validation_singlecell3_violin.png"), p1, width = 100, height = 50, dpi=200, units = "cm",limitsize = FALSE)


smy2 <- data.frame()
for(i in unique(smy[,1]))
{
	for(j in unique(smy[,2]))
	{
		smy2[paste0(i,j),"cellType"] <- i
		smy2[paste0(i,j),"gene"] <- j
		smy2[paste0(i,j),"value"] <- median(smy[as.character(smy[,1])==i&as.character(smy[,2])==j,3])
	}
}

smy2[,"group"] <- smy2[,"value"] > 0

library(ggplot2)
p1 <- ggplot(smy2,aes(x = cellType, y = value, fill=group))+
	geom_bar(stat='identity')+
	geom_hline(yintercept=0, color = "grey1", linewidth=0.8, linetype = 'dotted')+
	ggtitle("Immune Dictionary (single cell cytokine treatment)")+
	xlab(" ")+
	ylab("Act Diff\n(median)")+
	theme_bw()+ 
	theme(
	  panel.grid = element_blank(),
	  panel.background = element_blank(),
	  plot.title = element_text(hjust = 0.5),
	  axis.text = element_text(size=10,colour = "black"),
	  axis.title = element_text(size=10,colour = "black"),
	  legend.position="none"
	)+facet_wrap(~gene, ncol = 3)

ggsave(paste0(validationPath,"validation_singlecell3_bar.png"), p1, width = 100, height = 50, dpi=200, units = "cm",limitsize = FALSE)


new_symbol <- as.matrix(read.csv("../../Secretome/data/new_symbol_immunedictionary.txt",as.is=T,sep="\t"))
fig1d <- read.csv("../../Secretome/data/immunedictionary_Fig1d.tsv",as.is=T,sep="\t")
fig1d[fig1d[,2]=="ADSF",2] <- "ADSF (Resistin)"
fig1d[fig1d[,2]=="IGF-1",2] <- "IGF-I"

smy3 <- data.frame()
for(i in unique(smy[,2]))
{
	smy3[i,"gene"] <- i
	smy3[i,"act_diff"] <- median(smy2[as.character(smy2[,2])==i,3])
	
	i_alt <- new_symbol[new_symbol[,"new_symbol"]==i,"suppl_symbol"]
	smy3[i,"DEG_number"] <- median(fig1d[as.character(fig1d[,2])==i_alt,3])
}
cor_res <- cor.test(smy3[,2],smy3[,3])
rv <- round(cor_res$estimate,3)
pv <- signif(cor_res$p.value,3)

smy3[smy3[,2]<0.03,1] <- ""
library(ggplot2)
library(ggrepel)
p1 <- ggplot(smy3,aes(x = DEG_number, y = act_diff, label=gene))+
	geom_hline(yintercept=0, color = "grey1", linewidth=0.8, linetype = 'dotted')+
	geom_point(size=0.8)+
	geom_text_repel()+
	annotate("text", x=255, y=2.8, label=paste0("r=",rv,"\np=",pv) ) +
	#ggtitle(paste0("Immune Dictionary (single cell cytokine treatment), r=",rv,", p=",pv))+
	xlab("Differentially expressed gene number")+
	ylab("Activity Difference (median z score)")+
	theme_bw()+ 
	theme(
	  panel.grid = element_blank(),
	  panel.background = element_blank(),
	  plot.title = element_text(hjust = 0.5),
	  axis.text = element_text(size=10,colour = "black"),
	  axis.title = element_text(size=10,colour = "black"),
	  legend.position="none"
	)+geom_smooth(method="lm",se=F)

ggsave(paste0(validationPath,"validation_singlecell3_scatter.png"), p1, width = 10, height = 10, dpi=200, units = "cm",limitsize = FALSE)




}

