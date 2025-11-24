source("Secretome_s0_path.R")

#SWARM -t 2 -g 200 --time 50:00:00
args = commandArgs(trailingOnly=TRUE)
sigName <- args[1]
cancer <- args[2]

inputPath <- paste0(dataValPath,"scRNAseq/")
outputPath <- paste0(validationPath,"scRNAseq/")

cdata_T <- as.matrix(read.csv(gzfile(paste0(outputPath,cancer,"/expr.csv.gz")),row.names=1,check.names=F))	
rownames(cdata_T) <- transferSymbol(rownames(cdata_T))
cdata_T <- rm_duplicates(cdata_T)


#cdata_T <- cdata_T - rowMeans(cdata_T, na.rm=TRUE)


celltypes <- sapply(strsplit(colnames(cdata_T),".",fixed=T),function(x) return(x[1]))
celltypes <- gsub(",","",celltypes)

for(celltype in unique(celltypes))
{
	cdata_T_cellType <- cdata_T[,celltypes==celltype, drop=F]
	
	if(ncol(cdata_T_cellType)<10) next
	
	cdata_T_cellType_minusBG <- cdata_T_cellType - rowMeans(cdata_T_cellType, na.rm=T)
	
	cdata_T[,celltypes==celltype] <- cdata_T_cellType_minusBG
}

cdata_T <- round(cdata_T,3)


library(SecAct)

ref <- paste0(finalSignaturesPath,sigName,".tsv.gz")
lambda <- as.numeric(readLines(paste0(finalSignaturesPath,sigName,"_lambda.txt")))

res <- SecAct.activity.inference(
	inputProfile = cdata_T,
	is.differential = TRUE,
	sigMatrix = ref, 
	is.group.sig = ifelse(sigName=="SecAct",TRUE,FALSE),
	lambda = lambda,
	nrand = 1000
)

save(res, file = paste0(outputPath,cancer,"/",sigName,".RData"))







###################
# ROC calculation #
###################
inputPath <- paste0(dataValPath,"scRNAseq/")
outputPath <- paste0(validationPath,"scRNAseq/")

sp_tf <- list(
	IFN.STAT1=list(c("IFN1","IFNG","IFNL","IFNL1","IL21","IL27") , "STAT1"),
	#IFN1.STAT2=list(c("IFNAR2") , "STAT2"),
	IL6and10fam.STAT3=list(c("IL6","LIF","OSM","IL10","IL22","IL9","IL21","IL27","IL17A","IL23A") , "STAT3"),
	IL12fam.STAT4=list(c("IL12","IL23A","IL27") , "STAT4"), 
	IL2fam.STAT5=list(c("IL2","IL3","IL7","IL9","IL15","IL27") , "STAT5"), 
	#IL4.STAT6=list(c("IL4R") , "STAT6"), #IL4, IL13
	IL1andTNF.RELA=list(c("IL1A","IL1B","TNF") , "RELA"),
	TNFSFfam.RELB=list(c("LTA","CD40LG","TNFSF11","TNFSF12","TNFSF13B","TNFSF14") , "RELB"),
	BMPfam.SMAD1=list(c("BMP2","BMP3","BMP4","BMP5","BMP6","BMP7","BMP8A") , "SMAD1"),
	TGFBfam.SMAD2=list(c("INHBA","GDF11","TGFB1","TGFB2","TGFB3") , "SMAD2"),
	TGFBfam.SMAD3=list(c("INHBA","GDF11","TGFB1","TGFB2","TGFB3") , "SMAD3"),
	BMPandTGFBfam.SMAD4=list(c("BMP2","BMP3","BMP4","BMP5","BMP6","BMP7","BMP8A","INHBA","GDF11","TGFB1","TGFB2","TGFB3") , "SMAD4")
	#TGFB.SMAD5=list(c("BMP2","BMP3","BMP4","BMP5","BMP6","BMP7","BMP8A") , "SMAD5") # few cells have SMAD5 activity
)

load(paste0(outputPath,cancer,"/",sigName,".RData"))
sp <- res$zscore

tf <- read.table(paste0(inputPath,cancer,".Rabit.Cistrome.t.feature.gz"),sep="\t",check.names=F)
tf <- t(tf)
rownames(tf) <- sapply(strsplit(rownames(tf),".",fixed=T),function(x) return(x[2]))

olp <- intersect(colnames(sp),colnames(tf))
sp_olp <- sp[,olp]
tf_olp <- tf[,olp]


celltypes <- sapply(strsplit(olp,".",fixed=T),function(x) return(x[1]))
celltypes <- gsub(",","",celltypes)

summaryTable <- data.frame()
for(celltype in unique(celltypes))
{
	sp_olp_sub <- sp_olp[,celltypes==celltype]
	tf_olp_sub <- tf_olp[,celltypes==celltype,drop=F]
	if(ncol(tf_olp_sub)<2) next
	
	tf_olp_sub <- tf_olp_sub[rowSums(tf_olp_sub)!=0,]
	
	for(i in names(sp_tf))
	{
		if(sum(rownames(sp_olp_sub)%in%sp_tf[[i]][[1]])==0) next
		
		sp_act <- colMeans(sp_olp_sub[rownames(sp_olp_sub)%in%sp_tf[[i]][[1]],,drop=F])
		
		if(!sp_tf[[i]][[2]]%in%rownames(tf_olp_sub)) next
		tf_act <- colMeans(tf_olp_sub[sp_tf[[i]][[2]],,drop=F])
		
		
		if(sum(tf_act>0)==0 | sum(tf_act<=0)==0) next
		fg.df <- data.frame(sp_act=sp_act,tf_act=tf_act)
		fg.df <- fg.df[fg.df[,2]!=0,]
		fg.df[,2] <- fg.df[,2]>0
		
		
		if(length(unique(fg.df[,2]))==1) next
		if(sum(fg.df[,2]==TRUE)<3) next
		if(sum(fg.df[,2]==FALSE)<3) next
		if(nrow(fg.df)<10) next
		
		
		library(ROCR)
		predM <- prediction(fg.df$sp_act, fg.df$tf_act)
		roc = performance(predM, measure = "auc")
		roc2 = performance(predM, measure = "aucpr")
		AUPRC_baseline <- sum(fg.df[,2]==TRUE)/nrow(fg.df)
		
		summaryTable[paste0(cancer,celltype,i),"cancer"] <- cancer
		summaryTable[paste0(cancer,celltype,i),"celltype"] <- celltype
		summaryTable[paste0(cancer,celltype,i),"SP_TF"] <- i
		summaryTable[paste0(cancer,celltype,i),"AUC"] <- roc@ y.values
		summaryTable[paste0(cancer,celltype,i),"AUPRC"] <- roc2@ y.values
		summaryTable[paste0(cancer,celltype,i),"AUPRC_baseline"] <- AUPRC_baseline
	}
}

write.csv(summaryTable,paste0(outputPath,cancer,"/",sigName,"_SP_TF.csv"),quote=F)
