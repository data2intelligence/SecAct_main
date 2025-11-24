source("Secretome_s0_path.R")

inputPath <- paste0(dataValPath,"VEGF/")
outputPath <- paste0(validationPath,"VEGF/")

#SWARM -t 2 -g 50 --time 3:00:00
args = commandArgs(trailingOnly=TRUE)
sigName <- args[1]
item <- args[2]

dir.create(paste0(outputPath,item))

gene <- strsplit(item,"_") [[1]][1]
dataset <- strsplit(item,"_") [[1]][2]


# read diff expression
if(item=="VEGFA_GSE72951")
{
	f3 <- as.matrix(read.csv(gzfile(paste0(inputPath,"GSE72951.self_subtract.gz")),as.is=T,sep="\t",check.names=F))
	f33 <- read.csv(paste0(inputPath,"GSE72951.OS.Bevacizumab"),as.is=T,sep="\t")
}else{
	f3 <- as.matrix(read.csv(gzfile(paste0(inputPath,"E-MTAB-3267.norm_subtract.gz")),as.is=T,sep="\t",check.names=F))
	f33 <- read.csv(paste0(inputPath,"E-MTAB-3267.PFS"),as.is=T,sep="\t")
}

rownames(f3) <- transferSymbol(rownames(f3))
f3 <- rm_duplicates(f3)
	

# calculate activity
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
	
	if(sigName=="SecAct") save(res, file = paste0(outputPath,"/",item,"/",sigName,"_res.RData"))
	
	f3 <- res$zscore
}


# calculate risk score
library(survival)

if(sigName=="ReceptorExp")
{
	gene <- c("FLT1", "KDR", "NRP1")
}else if(sigName=="LRsumExp"){
	gene <- c("VEGFA","FLT1", "KDR", "NRP1")
}else{
	gene <- c("VEGFA")
}

olp <- intersect(rownames(f33),colnames(f3))
neww <- cbind(f33[olp,], value=colMeans(f3[rownames(f3)%in%gene,olp,drop=F]))

if(item%in%c("VEGFA_GSE72951"))
{
	coxmodel_fit <- coxph(Surv(OS, Event) ~ value, data = neww)
}else{
	coxmodel_fit <- coxph(Surv(PFS, Event) ~ value, data = neww)
}
coxmodel_obj <- summary(coxmodel_fit)
z <- coxmodel_obj$coefficients[1,"z"]
save(z, file = paste0(outputPath,"/",item,"/",sigName,".RData"))
