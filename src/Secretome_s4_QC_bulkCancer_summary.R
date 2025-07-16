source("Secretome_s0_path.R")
source("../../SpaCE/code/ifun.R")
source("importHPA.R")

dataPath <- "/data/rub2/data/TCGA_firebrowse/mRNAseq_gene_exp_symbol"
cancers <- list.files(dataPath)
cancers <- setdiff(cancers, "LAML")
cancers1 <- cancers

dataPath <- "/data/rub2/data/ICGC/"
allFiles <- list.files(dataPath)
allFiles <- allFiles[grepl("seq.expression",allFiles)]
cancers <- gsub(".seq.expression.gz","",allFiles)
cancers <- setdiff(cancers,"CLLE-ES") 
cancers <- setdiff(cancers,c("BPLL-FR","PAEN-AU","PRAD-FR")) # sample <40
cancers2 <- cancers
#CLLE-ES	Chronic Lymphocytic Leukemia - ES	Spain
#MALY-DE	Malignant Lymphoma - DE	Germany

cancers <- c(cancers1,cancers2)

version <- "vst_condition_logUMI_cellType"


comb <- data.frame()
for(cancer in cancers)
{
	combCancer <- data.frame()
	
	sts <- unique(meta[meta[,"Method"]=="Visium",1])

	for(st in sts)
	{
		QCPath.st <- paste0(QCPath,"/",st,"/")
		dir.create(QCPath.st)
		
		sampleNames <- meta[meta[,"Study"]==st&meta[,"Method"]=="Visium"&meta[,"Organism"]=="Human","Sample_Name"]
		for(sampleName in sampleNames)
		{	
			QCPath.st.sample <- paste0(QCPath.st,sampleName,"/")
			dir.create(QCPath.st.sample)
			
			#if(paste0(st,"@",sampleName)%in%smry[,1]) next
			
			if(grepl("-",cancer))
			{
				#if(!file.exists(paste0(QCPath.st.sample,"ICGC_filter_",cancer,".csv"))) print(QCPath.st.sample)
				smyMat <- read.csv(paste0(QCPath.st.sample,"ICGC_filter_",cancer,"_",version,".csv"),row.names=1)
				cancer_alt <- paste0("ICGC_",cancer)
			}else{
				#if(!file.exists(paste0(QCPath.st.sample,"TCGA_Xena_filter_log_",cancer,"_weighted.csv"))) print(QCPath.st.sample)}}}}
				smyMat <- read.csv(paste0(QCPath.st.sample,"TCGA_Xena_filter_log_",cancer,"_",version,".csv"),row.names=1)
				cancer_alt <- paste0("TCGA_",cancer)
			}
			combCancer[rownames(smyMat),paste0(st,"@",sampleName)] <- smyMat[,1]
		}
	}
	combCancer <- apply(combCancer,1,function(x) mean(x,na.rm=T) )
	
	comb[names(combCancer),cancer_alt] <- combCancer
}


comb.m <- reshape2::melt(as.matrix(comb))
comb.m <- comb.m[!is.na(comb.m[,3]),]
comb.m <- cbind(comb.m,Group="NA")

comb.m[comb.m[,1]%in%IPs,"Group"] <- "Intracellular"
comb.m[comb.m[,1]%in%MPs,"Group"] <- "Membrane"
comb.m[comb.m[,1]%in%SPs,"Group"] <- "Secreted"

comb.m <- comb.m[!comb.m[,"Group"]%in%c("NA"),]
comb.m <- comb.m[!comb.m[,"Group"]%in%c("Membrane","NA"),]


library(ggplot2)
library(ggsignif)
p2 <- ggplot(comb.m,aes(x = Group, y = value))+
	geom_hline(yintercept=0, color = "grey", linewidth=0.6, linetype="dashed")+
	geom_violin(aes(group=Group, fill=Group), trim=FALSE)+
	geom_boxplot(width=0.1, outlier.shape=NA)+
	#geom_signif( comparisons = list(c("Secreted", "Intracellular"),c("Membrane", "Intracellular")) ,y_position = c(0.72,-0.5), test = "wilcox.test", test.args ="greater")+
	#scale_fill_manual(values = c("skyblue","green","#ff8080"))+
	geom_signif( comparisons = list(c("Secreted", "Intracellular")) ,y_position = c(0.72), test = "wilcox.test", test.args ="greater")+
	scale_fill_manual(values = c("skyblue","#ff8080"))+
	ggtitle("Correlation between secreted proteins' expression and signature scores")+
	ylab("Pearson r value")+
	theme_bw()+ 
	theme(
	  panel.grid = element_blank(),
	  panel.background = element_blank(),
	  plot.title = element_text(hjust = 0.5, size=18),
	  axis.text = element_text(size=15,colour = "black"),
	  axis.text.x = element_blank(),
	  axis.ticks.x = element_blank(),
	  axis.title.x = element_blank(),
	  axis.title = element_text(size=15,colour = "black"),
	  legend.position="bottom",
	  legend.text=element_text(size=15),
	  legend.title=element_text(size=15),
	  strip.text = element_text(size = 13)
	) + facet_wrap(~ Var2, ncol =11) #scales = "free"
	
#ggsave(paste0(QCSummaryPath,"Summary_TCGA.png"), p2, width = 55, height = 19, dpi=300, units = "cm", limitsize = FALSE)
ggsave(paste0(QCSummaryPath,"Summary_TCGA_ICGC_1k_",version,".png"), p2, width = 55, height = 23, dpi=200, units = "cm", limitsize = FALSE)












comb <- data.frame()

sts <- unique(meta[meta[,"Method"]=="Visium",1])

for(st in sts)
{
	QCPath.st <- paste0(QCPath,"/",st,"/")
	dir.create(QCPath.st)
	
	QCFilterPath.st <- paste0(QCFilterPath,"/",st,"/")
	dir.create(QCFilterPath.st)
	
	sampleNames <- meta[meta[,"Study"]==st&meta[,"Method"]=="Visium"&meta[,"Organism"]=="Human","Sample_Name"]
	for(sampleName in sampleNames)
	{	
		QCPath.st.sample <- paste0(QCPath.st,sampleName,"/")
		dir.create(QCPath.st.sample)
		QCFilterPath.st.sample <- paste0(QCFilterPath.st,sampleName,"/")
		dir.create(QCFilterPath.st.sample)
					
		smyMat <- read.csv(paste0(QCFilterPath.st.sample,"TCGA_Xena_filter_log_Secreted_filterBar_comp.csv"),row.names=1)
		
		comb[rownames(smyMat),paste0(st,"@",sampleName)] <- smyMat[,"new_0.25"]
	}
}
	
comb[!is.na(comb)] <- 1
comb <- as.matrix(comb)

apply(comb,2,function(x) sum(x,na.rm=T))

sort(apply(comb,2,function(x) sum(x,na.rm=T)),decreasing=T)


x <- comb["TCGA_GBM",]
names(x) <- NULL
x

sum(x[1:618],na.rm=T)/618
sum(x[619:1032],na.rm=T)/410








vsos <- c(
"TCGA_Xena_genomeWide_log","TCGA_Xena_genomeWide_loglog",
"TCGA_Xena_filter_log","TCGA_Xena_filter_loglog",
"TCGA_CIDE_genomeWide_log","TCGA_CIDE_genomeWide_loglog",
"TCGA_CIDE_filter_log","TCGA_CIDE_filter_loglog"
)

for(vso in vsos[3])
{

comb <- data.frame()
for(cancer in cancers1)
{
	combCancer <- data.frame()
	
	sts <- unique(meta[meta[,"Method"]=="Visium",1])
	for(st in sts)
	{
		QCPath.st <- paste0(QCPath,"/",st,"/")
		dir.create(QCPath.st)
		
		sampleNames <- meta[meta[,"Study"]==st&meta[,"Method"]=="Visium"&meta[,"Organism"]=="Human","Sample_Name"]
		for(sampleName in sampleNames)
		{	
			QCPath.st.sample <- paste0(QCPath.st,sampleName,"/")
			dir.create(QCPath.st.sample)
			
			if(!file.exists(paste0(QCPath.st.sample,vso,"TCGA_Xena_filter_log_",cancer,".csv"))) next
			smyMat <- read.csv(paste0(QCPath.st.sample,vso,"TCGA_Xena_filter_log_",cancer,".csv"),row.names=1)
			combCancer[rownames(smyMat),paste0(st,"@",sampleName)] <- smyMat[,1]
		}
	}
	combCancer <- apply(combCancer,1,function(x) mean(x,na.rm=T) )
	
	comb[names(combCancer),cancer] <- combCancer
}


comb.m <- reshape2::melt(as.matrix(comb))
comb.m <- comb.m[!is.na(comb.m[,3]),]
comb.m <- cbind(comb.m,Group="NA")

comb.m[comb.m[,1]%in%anno[["IPs"]],"Group"] <- "Intracellular"
comb.m[comb.m[,1]%in%anno[["MPs"]],"Group"] <- "Membrane"
comb.m[comb.m[,1]%in%anno[["SPs"]],"Group"] <- "Secreted"
comb.m <- comb.m[comb.m[,"Group"]!="NA",]


library(ggplot2)
library(ggsignif)
p2 <- ggplot(comb.m,aes(x = Group, y = value))+
	geom_hline(yintercept=0, color = "grey", linewidth=0.6, linetype="dashed")+
	geom_violin(aes(group=Group, fill=Group),trim=FALSE)+
	geom_boxplot(width=0.15,outlier.shape=NA)+
	geom_signif( comparisons = list(c("Secreted", "Intracellular")) ,y_position = c(0.75) )+
	xlab("Group")+
	ylab("Correlation between activity and expression")+
	theme_bw()+ 
	theme(
	  panel.grid = element_blank(),
	  panel.background = element_blank(),
	  axis.text = element_text(size=10,colour = "black"),
	  axis.text.x = element_blank(),
	  axis.title = element_text(size=10,colour = "black"),
	  legend.position="right",
	  strip.text = element_text(size = 10)
	) + facet_wrap(~ Var2, ncol =10) #scales = "free"
	
ggsave(paste0(QCSummaryPath,"Summary_",vso,".png"), p2, width = 50, height = 22, dpi=200, units = "cm", limitsize = FALSE)

}






cancers <- list.files("/data/rub2/project/Secretome/results/QC/BRCA_10x_Datasets/Version1.0.0_Breast.Cancer_rep1/")
cancers <- gsub(".csv","",cancers)
cancers <- cancers[grepl("_filtered_",cancers)]
cancers <- cancers[c(31:60,1:10)]
cancers <- sapply(strsplit(cancers,"_filtered_",fixed=T),function(x) return(x[2]))
cancers30 <- cancers[1:30]
cancers10 <- cancers[31:40]

comb1 <- data.frame()
comb2 <- data.frame()
comb3 <- data.frame()
for(cancer in cancers)
{
	combCancer1 <- data.frame()
	combCancer2 <- data.frame()
	combCancer3 <- data.frame()
	
	sts <- unique(meta[meta[,"Method"]=="Visium",1])
	for(st in sts)
	{
		QCPath.st <- paste0(QCPath,"/",st,"/")
		dir.create(QCPath.st)
		
		sampleNames <- meta[meta[,"Study"]==st&meta[,"Method"]=="Visium"&meta[,"Organism"]=="Human","Sample_Name"]
		for(sampleName in sampleNames)
		{	
			QCPath.st.sample <- paste0(QCPath.st,sampleName,"/")
			dir.create(QCPath.st.sample)
			
			if(cancer%in%cancers30)
			{
				smyMat1 <- read.csv(paste0(QCPath.st.sample,"TCGA_filtered_",cancer,".csv"),row.names=1)
				smyMat2 <- read.csv(paste0(QCPath.st.sample,"TCGA_SigNormSPs_filtered_",cancer,".csv"),row.names=1)
				smyMat3 <- read.csv(paste0(QCPath.st.sample,"TCGA_SigNormIPs_filtered_",cancer,".csv"),row.names=1)
			}else{
				smyMat1 <- read.csv(paste0(QCPath.st.sample,"ICGC_filtered_",cancer,".csv"),row.names=1)
				smyMat2 <- read.csv(paste0(QCPath.st.sample,"ICGC_SigNormSPs_filtered_",cancer,".csv"),row.names=1)
				smyMat3 <- read.csv(paste0(QCPath.st.sample,"ICGC_SigNormIPs_filtered_",cancer,".csv"),row.names=1)
			}
			
			combCancer1[rownames(smyMat1),paste0(st,"@",sampleName)] <- smyMat1[,1]
			combCancer2[rownames(smyMat2),paste0(st,"@",sampleName)] <- smyMat2[,1]
			combCancer3[rownames(smyMat3),paste0(st,"@",sampleName)] <- smyMat3[,1]
		}
	}
	
	combCancer1 <- apply(combCancer1,1,function(x) mean(x,na.rm=T) )
	combCancer2 <- apply(combCancer2,1,function(x) mean(x,na.rm=T) )
	combCancer3 <- apply(combCancer3,1,function(x) mean(x,na.rm=T) )
	
	comb1[names(combCancer1),cancer] <- combCancer1
	comb2[names(combCancer2),cancer] <- combCancer2
	comb3[names(combCancer3),cancer] <- combCancer3
}

comb <- comb1

comb.m <- reshape2::melt(as.matrix(comb))
comb.m <- comb.m[!is.na(comb.m[,3]),]
comb.m <- cbind(comb.m,Group="NA")

comb.m[comb.m[,1]%in%anno[["IPs"]],"Group"] <- "Intracellular"
comb.m[comb.m[,1]%in%anno[["MPs"]],"Group"] <- "Membrane"
comb.m[comb.m[,1]%in%anno[["SPs"]],"Group"] <- "Secreted"
comb.m <- comb.m[comb.m[,"Group"]!="NA",]

comb.m1 <- comb.m

comb <- comb2

comb.m <- reshape2::melt(as.matrix(comb))
comb.m <- comb.m[!is.na(comb.m[,3]),]
comb.m <- cbind(comb.m,Group="NA")

comb.m[comb.m[,1]%in%anno[["IPs"]],"Group"] <- "Intracellular"
comb.m[comb.m[,1]%in%anno[["MPs"]],"Group"] <- "Membrane"
comb.m[comb.m[,1]%in%anno[["SPs"]],"Group"] <- "Secreted"
comb.m <- comb.m[comb.m[,"Group"]!="NA",]

comb.m2 <- comb.m

comb <- comb3

comb.m <- reshape2::melt(as.matrix(comb))
comb.m <- comb.m[!is.na(comb.m[,3]),]
comb.m <- cbind(comb.m,Group="NA")

comb.m[comb.m[,1]%in%anno[["IPs"]],"Group"] <- "Intracellular"
comb.m[comb.m[,1]%in%anno[["MPs"]],"Group"] <- "Membrane"
comb.m[comb.m[,1]%in%anno[["SPs"]],"Group"] <- "Secreted"
comb.m <- comb.m[comb.m[,"Group"]!="NA",]

comb.m3 <- comb.m


comb.m1[,"Group"] <- paste0("NormalNo_",comb.m1[,"Group"])
comb.m2[,"Group"] <- paste0("NormalSP_",comb.m2[,"Group"])
comb.m3[,"Group"] <- paste0("NormalIP_",comb.m3[,"Group"])


comb.m <- rbind(comb.m1,comb.m2)
comb.m <- rbind(comb.m,comb.m3)


library(ggplot2)
library(ggsignif)
p2 <- ggplot(comb.m,aes(x = Group, y = value))+
	geom_hline(yintercept=0, color = "grey", linewidth=0.6, linetype="dashed")+
	geom_violin(aes(group=Group, fill=Group),trim=FALSE)+
	geom_boxplot(width=0.15,outlier.shape=NA)+
	geom_signif( comparisons = list(c("Secreted", "Intracellular")) ,y_position = c(0.75) )+
	xlab("Group")+
	ylab("Correlation between activity and expression")+
	theme_bw()+ 
	theme(
	  panel.grid = element_blank(),
	  panel.background = element_blank(),
	  axis.text = element_text(size=10,colour = "black"),
	  axis.text.x = element_blank(),
	  axis.title = element_text(size=10,colour = "black"),
	  legend.position="right",
	  strip.text = element_text(size = 10)
	) + facet_wrap(~ Var2, ncol =10) #scales = "free"
	
ggsave(paste0(QCSummaryPath,"Summary_TCGA_ICGC_normal.png"), p2, width = 50, height = 22, dpi=200, units = "cm", limitsize = FALSE)











LRdb <- data.frame()

data_list <- rio::import_list("/data/rub2/project/SpaCE/data/LR/Receptor_Cytokine.xlsx")

for(i in names(data_list))
{
	temp <- data_list[[i]]
	
	for(j in 1:nrow(temp))
	{	
		items_1 <- unlist( strsplit(temp[j,2],"+",fixed=T) )
		items_1 <- unlist( strsplit(items_1,",",fixed=T) )
		items_1 <- gsub(" ","",items_1,fixed=T)
		items_1 <- unique(items_1)
		
		items_2 <- unlist( strsplit(temp[j,3],"+",fixed=T) )
		items_2 <- unlist( strsplit(items_2,",",fixed=T) )
		items_2 <- gsub(" ","",items_2,fixed=T)
		items_2 <- unique(items_2)
		
		items <- as.matrix(expand.grid(items_1,items_2))
		
		LRdb[paste0(items[,1],"_",items[,2]),"Ligand"] <- items[,1]
		LRdb[paste0(items[,1],"_",items[,2]),"Receptor"] <- items[,2]
	}
}

LRdb1 <- LRdb

LRdb <- data.frame()

data_list <- rio::import_list("/data/rub2/project/SpaCE/data/LR/Receptor_Chemokine.xlsx")

for(i in names(data_list))
{
	temp <- data_list[[i]]
	
	for(j in 1:nrow(temp))
	{	
		items_1 <- unlist( strsplit(temp[j,2],"+",fixed=T) )
		items_1 <- unlist( strsplit(items_1,",",fixed=T) )
		items_1 <- gsub(" ","",items_1,fixed=T)
		items_1 <- gsub("-","",items_1,fixed=T)
		items_1 <- unique(items_1)
		
		items_2 <- unlist( strsplit(temp[j,3],"+",fixed=T) )
		items_2 <- unlist( strsplit(items_2,",",fixed=T) )
		items_2 <- gsub(" ","",items_2,fixed=T)
		items_2 <- gsub("-","",items_2,fixed=T)
		items_2 <- unique(items_2)
				
		items <- as.matrix(expand.grid(items_1,items_2))
		
		LRdb[paste0(items[,1],"_",items[,2]),"Ligand"] <- items[,1]
		LRdb[paste0(items[,1],"_",items[,2]),"Receptor"] <- items[,2]
	}
}

LRdb2 <- LRdb

LRdb <- data.frame()

data_list <- rio::import_list("/data/rub2/project/SpaCE/data/LR/Receptor_Growth_Factor.xlsx")

for(i in names(data_list))
{
	temp <- data_list[[i]]
	
	for(j in 1:nrow(temp))
	{	
		items_1 <- unlist( strsplit(temp[j,2],"+",fixed=T) )
		items_1 <- unlist( strsplit(items_1,",",fixed=T) )
		items_1 <- gsub(" ","",items_1,fixed=T)
		items_1 <- gsub("-","",items_1,fixed=T)
		items_1 <- unique(items_1)
		
		items_2 <- unlist( strsplit(temp[j,3],"+",fixed=T) )
		items_2 <- unlist( strsplit(items_2,",",fixed=T) )
		items_2 <- gsub(" ","",items_2,fixed=T)
		items_2 <- gsub("-","",items_2,fixed=T)
		items_2 <- unique(items_2)
				
		items <- as.matrix(expand.grid(items_1,items_2))
		
		LRdb[paste0(items[,1],"_",items[,2]),"Ligand"] <- items[,1]
		LRdb[paste0(items[,1],"_",items[,2]),"Receptor"] <- items[,2]
	}
}

LRdb3 <- LRdb



SIGs <- c("SIGLEC1","CD22","CD33","MAG", paste0("SIGLEC",5:17))


geneGroups <- list(
	ligand_cytokine = unique(transferSymbol(LRdb1[,1])),
	ligand_chemokine = unique(transferSymbol(LRdb2[,1])),
	ligand_growthfactor = unique(transferSymbol(LRdb3[,1])),
	receptor_cytokine = unique(transferSymbol(LRdb1[,2])),
	receptor_chemokine = unique(transferSymbol(LRdb2[,2])),
	receptor_growthfactor = unique(transferSymbol(LRdb3[,2])),
	SIGLECs = unique(transferSymbol(SIGs))
)




cancers <- list.files("/data/rub2/project/Secretome/results/QC/BRCA_10x_Datasets/Version1.0.0_Breast.Cancer_rep1/")
cancers <- gsub(".csv","",cancers)
cancers <- cancers[grepl("_filtered_",cancers)]
cancers <- cancers[c(11:40,1:10)]


comb <- data.frame()
for(cancer in cancers)
{
	combCancer <- data.frame()
	
	sts <- unique(meta[meta[,"Method"]=="Visium",1])
	for(st in sts)
	{
		QCPath.st <- paste0(QCPath,"/",st,"/")
		dir.create(QCPath.st)
		
		sampleNames <- meta[meta[,"Study"]==st&meta[,"Method"]=="Visium"&meta[,"Organism"]=="Human","Sample_Name"]
		for(sampleName in sampleNames)
		{	
			QCPath.st.sample <- paste0(QCPath.st,sampleName,"/")
			dir.create(QCPath.st.sample)
			
			if(!file.exists(paste0(QCPath.st.sample,cancer,".csv"))) next
			smyMat <- read.csv(paste0(QCPath.st.sample,cancer,".csv"),row.names=1)
			combCancer[rownames(smyMat),paste0(st,"@",sampleName)] <- smyMat[,1]
		}
	}
	combCancer <- apply(combCancer,1,function(x) mean(x,na.rm=T) )
	
	comb[names(combCancer),cancer] <- combCancer
}


comb.m <- reshape2::melt(as.matrix(comb))
comb.m <- comb.m[!is.na(comb.m[,3]),]
comb.m <- cbind(comb.m,Group="NA")

comb.m[comb.m[,1]%in%anno[["IPs"]],"Group"] <- "Intracellular"
comb.m[comb.m[,1]%in%anno[["MPs"]],"Group"] <- "Membrane"
comb.m[comb.m[,1]%in%anno[["SPs"]],"Group"] <- "Secreted"

comb.m_exp <- comb.m

for(i in names(geneGroups))
{
	genes <- geneGroups[[i]]
	
	comb.m_sub <- comb.m[comb.m[,1]%in%genes,]
	comb.m_sub[,4] <- i
	
	comb.m_exp <- rbind(comb.m_exp, comb.m_sub) 
}

comb.m_exp -> comb.m
comb.m <- comb.m[comb.m[,"Group"]!="NA",]


comb.m[,"Group"] <- factor(comb.m[,"Group"], 
	levels=c("Intracellular","Membrane","receptor_cytokine","receptor_chemokine","receptor_growthfactor",
	"Secreted","ligand_cytokine","ligand_chemokine","ligand_growthfactor","SIGLECs")
)


library(ggplot2)
library(ggsignif)
p2 <- ggplot(comb.m,aes(x = Group, y = value))+
	geom_hline(yintercept=0, color = "grey", linewidth=0.6, linetype="dashed")+
	geom_violin(aes(group=Group, fill=Group),trim=FALSE)+
	scale_fill_manual(name="Group", values = c(RColorBrewer::brewer.pal(9,"Set1"), "black") )+
	geom_boxplot(width=0.15,outlier.shape=NA)+
	#geom_signif( comparisons = list(c("Secreted", "Intracellular")) ,y_position = c(0.75) )+
	xlab("Group")+
	ylab("Correlation between activity and expression")+
	theme_bw()+ 
	theme(
	  panel.grid = element_blank(),
	  panel.background = element_blank(),
	  axis.text = element_text(size=10,colour = "black"),
	  axis.text.x = element_blank(),
	  axis.title = element_text(size=10,colour = "black"),
	  legend.position="right",
	  strip.text = element_text(size = 10)
	) + facet_wrap(~ Var2, ncol =8) #scales = "free"
	
ggsave(paste0(QCSummaryPath,"Summary_TCGA_ICGC_expand.png"), p2, width = 80, height = 40, dpi=200, units = "cm", limitsize = FALSE)





























cancers <- list.files("/data/rub2/project/Secretome/results/QC/BRCA_10x_Datasets/Version1.0.0_Breast.Cancer_rep1/")
cancers <- setdiff(cancers,c("filterBar.csv","filterBar.png"))
cancers <- gsub(".csv","",cancers)
cancers <- cancers[61:126]

sts <- unique(meta[meta[,"Method"]=="Visium",1])

comb <- data.frame()
for(cancer in cancers)
{
	combCancer <- data.frame()
	
	for(st in sts)
	{
		QCPath.st <- paste0(QCPath,"/",st,"/")
		dir.create(QCPath.st)
		
		sampleNames <- meta[meta[,"Study"]==st&meta[,"Method"]=="Visium"&meta[,"Organism"]=="Human","Sample_Name"]
		for(sampleName in sampleNames)
		{	
			QCPath.st.sample <- paste0(QCPath.st,sampleName,"/")
			dir.create(QCPath.st.sample)
			
			if(!file.exists(paste0(QCPath.st.sample,cancer,".csv"))) next
			smyMat <- read.csv(paste0(QCPath.st.sample,cancer,".csv"),row.names=1)
			combCancer[rownames(smyMat),paste0(st,"@",sampleName)] <- smyMat[,1]
		}
	}
	combCancer <- apply(combCancer,1,function(x) mean(x,na.rm=T) )
	
	comb[names(combCancer),cancer] <- combCancer
}


cancers <- sapply(strsplit(cancers[1:33],"_",),function(x) return(x[3]))

for(cancer in cancers)
{
	combCancer <- data.frame()
	
	for(st in sts)
	{
		QCPath.st <- paste0(QCPath,"/",st,"/")
		dir.create(QCPath.st)
		
		sampleNames <- meta[meta[,"Study"]==st&meta[,"Method"]=="Visium"&meta[,"Organism"]=="Human","Sample_Name"]
		for(sampleName in sampleNames)
		{	
			QCPath.st.sample <- paste0(QCPath.st,sampleName,"/")
			dir.create(QCPath.st.sample)
			
			if(!file.exists(paste0(QCPath.st.sample,"TCGA_genomeWide_",cancer,".csv"))) next
			if(!file.exists(paste0(QCPath.st.sample,"TCGA_filtered_",cancer,".csv"))) next
			smyMat_g <- read.csv(paste0(QCPath.st.sample,"TCGA_genomeWide_",cancer,".csv"),row.names=1)
			smyMat_f <- read.csv(paste0(QCPath.st.sample,"TCGA_filtered_",cancer,".csv"),row.names=1)
			
			smyMat <- smyMat_g[rownames(smyMat_g)%in%rownames(smyMat_f),,drop=F]
			combCancer[rownames(smyMat),paste0(st,"@",sampleName)] <- smyMat[,1]
		}
	}
	combCancer <- apply(combCancer,1,function(x) mean(x,na.rm=T) )
	
	comb[names(combCancer),paste0("TCGA_filteredAfterSig_",cancer)] <- combCancer
}


comb.m <- reshape2::melt(as.matrix(comb))
comb.m <- comb.m[!is.na(comb.m[,3]),]
comb.m <- cbind(comb.m,Group="NA")

source("importHPA.R")
comb.m[comb.m[,1]%in%anno[["IPs"]],"Group"] <- "Intracellular"
comb.m[comb.m[,1]%in%anno[["MPs"]],"Group"] <- "Membrane"
comb.m[comb.m[,1]%in%anno[["SPs"]],"Group"] <- "Secreted"

comb.m <- comb.m[comb.m[,"Group"]!="NA",]


library(ggplot2)
library(ggsignif)
p2 <- ggplot(comb.m,aes(x = Group, y = value))+
	geom_hline(yintercept=0, color = "grey", linewidth=0.6, linetype="dashed")+
	#geom_jitter(size=0.6)+
	geom_violin(aes(group=Group, fill=Group),trim=FALSE)+
	geom_boxplot(width=0.15,outlier.shape=NA)+
	geom_signif( comparisons = list(c("Secreted", "Intracellular")) ,y_position = c(0.75) )+
	xlab("Group")+
	ylab("Correlation between activity and expression")+
	theme_bw()+ 
	theme(
	  panel.grid = element_blank(),
	  panel.background = element_blank(),
	  axis.text = element_text(size=10,colour = "black"),
	  axis.text.x = element_blank(),
	  axis.title = element_text(size=10,colour = "black"),
	  legend.position="right",
	  strip.text = element_text(size = 10)
	) + facet_wrap(~ Var2, ncol =33) #scales = "free"
	
ggsave(paste0(QCSummaryPath,"Summary_TCGA.png"), p2, width = 200, height =30, dpi=200, units = "cm", limitsize = FALSE)
#ggsave(paste0(QCSummaryPath,"Summary_ICGC.png"), p2, width = 88, height =30, dpi=200, units = "cm", limitsize = FALSE)




#combStat <- data.frame()
#for(cancer in cancers)
#{
#	comb.m_sub <- comb.m[comb.m[,2]==paste0("CPTAC_Protein_",cancer),]
#	
#	x1 <- comb.m_sub[comb.m_sub[,"Group"]=="Secreted",3]
#	x2 <- comb.m_sub[comb.m_sub[,"Group"]=="Intracellular",3]
#	pv <- signif(wilcox.test(x1, x2)$p.value,2)
#	medianDif <- median(x1)-median(x2)
#	
#	combStat[paste0(cancer,"ST"),"cancer"]<- cancer
#	combStat[paste0(cancer,"ST"),"sig"]<- "ST"
#	combStat[paste0(cancer,"ST"),"ratio"]<- sum(x1>quantile(x2,0.95))/length(x1)
#}




source("Secretome_s0_path.R")
source("../../SpaCE/code/ifun.R")

cancers <- list.files("/data/rub2/project/Secretome/results/QC/BRCA_10x_Datasets/Version1.0.0_Breast.Cancer_rep1/")
cancers <- setdiff(cancers,c("filterBar.csv","filterBar.png"))
cancers <- gsub(".csv","",cancers)
cancers <- cancers[c(11:12,23:30)]
cancers <- substring(cancers,11)

sts <- unique(meta[meta[,"Method"]=="Visium",1])
for(cancer in cancers)
{
	for(st in sts)
	{
		QCPath.st <- paste0(QCPath,"/",st,"/")
		dir.create(QCPath.st)
		
		sampleNames <- meta[meta[,"Study"]==st&meta[,"Method"]=="Visium"&meta[,"Organism"]=="Human","Sample_Name"]
		for(sampleName in sampleNames)
		{	
			QCPath.st.sample <- paste0(QCPath.st,sampleName,"/")
			dir.create(QCPath.st.sample)
			
			ifile <- paste0(QCPath.st.sample,"CPTAC_RNA_",cancer,".csv")
			ofile <- paste0(QCPath.st.sample,"CPTAC_RNA_GenomeWide_",cancer,".csv")
			
			system(paste0("mv ",ifile," ",ofile))
		}
	}
}




smry <- data.frame()

sts <- unique(meta[,1])
for(st in sts)
{
	QCFilterPath.st <- paste0(QCFilterPath,st,"/")
	dir.create(QCFilterPath.st)
	
	sampleNames <- meta[meta[,"Study"]==st,"Sample_Name"]
	for(sampleName in sampleNames)
	{	
		QCFilterPath.st.sample <- paste0(QCFilterPath.st,sampleName,"/")
		dir.create(QCFilterPath.st.sample)
		
		combTest <- read.csv(paste0(QCFilterPath.st.sample,"TCGA_Xena_filter_log_","filterBar_Test.csv"), row.names=1, as.is=T)
		
		smry[rownames(combTest),paste0(st,"@",sampleName)] <- combTest[,3]
	}
}

smry <- as.matrix(smry)
smry <- smry=="Yes"

rowSums(smry)





comb <- data.frame()

sts <- unique(meta[meta[,"Method"]=="Visium",1])
for(st in sts)
{
	QCPath.st <- paste0(QCPath,st,"/")
	dir.create(QCPath.st)
	
	sampleNames <- meta[meta[,"Study"]==st&meta[,"Method"]=="Visium"&meta[,"Organism"]=="Human","Sample_Name"]
	for(sampleName in sampleNames)
	{	
		QCPath.st.sample <- paste0(QCPath.st,sampleName,"/")
		dir.create(QCPath.st.sample)
		
		combFrac <- read.csv(paste0(QCPath.st.sample,"filterBar_ratio.csv"), row.names=1, header=T)
		
		comb[paste0(st,"@",sampleName),rownames(combFrac)] <- combFrac[,1]
	}
}

