source("Secretome_s0_path.R")

inputPath <- paste0(dataValPath,"blocking/")
outputPath <- paste0(validationPath,"blocking/")

#SWARM -t 2 -g 50 --time 3:00:00
args = commandArgs(trailingOnly=TRUE)
sigName <- args[1]
item <- args[2]

dir.create(paste0(outputPath,item))

gene <- strsplit(item,"_") [[1]][1]
dataset <- strsplit(item,"_") [[1]][2]


# read diff expression
if(item%in%c("IL11_Widjaja2024","NTN1_GSE225691"))
{
	f3 <- read.csv(paste0(inputPath,item,".diff"),as.is=T,sep="\t",check.names=F)
	
}else if(item%in%c("TGFB_GSE174686","TGFB3_GSE174686")){
	f3 <- as.matrix(read.csv(gzfile(paste0(inputPath,"aTGFB_XOMA.diff.human.gz")),as.is=T,sep="\t",check.names=F))
	
	if(item%in%c("TGFB_GSE174686")) f3 <- f3[,"mouse.Pan",drop=F]
	if(item=="TGFB3_GSE174686") f3 <- f3[,"Pan-1+2",drop=F]
}else{
	allFiles <- list.files(inputPath)
	
	flag1 <- grepl(dataset,allFiles)
	flag2 <- grepl(".diff.1.cntmap",allFiles)
	
	fileName_full <- allFiles[flag1&flag2]
	fileName_full <- gsub(".diff.1.cntmap","",fileName_full)
		
	f3 <- as.matrix(read.csv(gzfile(paste0(inputPath,fileName_full,".diff.1")),as.is=T,sep="\t"))
}

# calculate z
if(sigName%in%sigNames)
{
	library(SecAct)
	
	inputProfile <- f3
	ref <- paste0(finalSignaturesPath,sigName,".tsv.gz")
	is.group.sig <- ifelse(sigName=="SecAct",TRUE,FALSE)
	lambda <- as.numeric(readLines(paste0(finalSignaturesPath,sigName,"_lambda.txt")))

	res <- SecAct.activity.inference(
		inputProfile = inputProfile,
		is.differential = TRUE,
		sigMatrix = ref, 
		is.group.sig = is.group.sig,
		lambda = lambda,
		nrand = 1000
	)
	
	f3 <- res$zscore
}else{ # expression
	f3 <- scale(f3)
}


# comb samples
f3 <- rbind(f3,TGFB=colMeans(f3[rownames(f3)%in%c("TGFB1","TGFB2","TGFB3"),,drop=F]))
		
f3 <- f3[,!grepl("Placebo",colnames(f3),fixed=T) & !is.null(colnames(f3)),drop=F]
f3 <- f3[,!grepl("non.responder",colnames(f3),fixed=T) & !is.null(colnames(f3)),drop=F]

if(item=="IFN1_GSE72754")
{
	f3 <- f3[,"IFNA.high.whole.blood",drop=F]
}
if(item=="IL1B_GSE80060")
{
	f3 <- f3[,c("IL1B.Day3.Canakinumab_100.0","IL1B.Day3.Canakinumab_90.0","IL1B.Day3.Canakinumab_70.0"),drop=F]
}
if(item=="TNFSF12_GSE42048")
{
	f3 <- f3[,c("TWEAK.ACHN.xenograft.72h","TWEAK.ACHN.xenograft.24h"),drop=F]
}
if(item=="TNF_GSE48498")
{
	f3 <- f3[,"TNFA.whole.blood.cells",drop=F]
}

f3 <- apply(f3,1,function(x) mean(x,na.rm=T))


# extract z
LRdb <- extract_experimental_varified_LR( paste0(dataValPath,"LR/") )
LRdb <- LRdb[!is.na(LRdb[,2]),]
LRdb <- rbind(LRdb, NTN1_DCC=c("NTN1","DCC"))
LRdb <- rbind(LRdb, NTN1_UNC5C=c("NTN1","UNC5C"))

if(sigName%in%c("ReceptorExp"))	
{
	if(gene=="TGFB") gene <- c("TGFB1","TGFB2","TGFB3")
	z <- mean(f3[unique(LRdb[LRdb[,1]%in%gene,2])], na.rm=TRUE)
}else if(sigName%in%c("LRsumExp")){
	if(gene=="TGFB") gene <- c("TGFB1","TGFB2","TGFB3")
	z <- mean(f3[c(gene,unique(LRdb[LRdb[,1]%in%gene,2]))], na.rm=TRUE)
}else{
	z <- f3[gene]
}


save(z, file = paste0(outputPath,"/",item,"/",sigName,".RData"))

