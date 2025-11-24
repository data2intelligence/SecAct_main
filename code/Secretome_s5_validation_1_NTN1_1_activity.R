source("Secretome_s0_path.R")

inputPath <- paste0(dataValPath,"NTN1/")
outputPath <- paste0(validationPath,"NTN1/")

#SWARM -t 2 -g 50 --time 3:00:00
args = commandArgs(trailingOnly=TRUE)
sampleName <- args[1]


library(SpaCET)
library(SecAct)

visiumPath <- paste0(inputPath,"GSE225690_RAW/",sampleName,"/") 
SpaCET_obj <- create.SpaCET.object.10X(visiumPath = visiumPath)

SpaCET_obj <- SpaCET.quality.control(SpaCET_obj)
SpaCET_obj <- SpaCET.deconvolution(SpaCET_obj, cancerType="UCEC", coreNo=15)

ref <- paste0(finalSignaturesPath,"SecAct.tsv.gz")

if(sampleName%in%c("01-034_C1D1","01-039_C1D1"))
{
	SpaCET_obj <- SecAct.activity.inference.ST(SpaCET_obj, sigMatrix=ref, lambda=5e05, is.group.sig=TRUE)
}else{
	if(sampleName=="01-034_C3D1")
	{
		SpaCET_obj_CTL <- readRDS(paste0(outputPath,"01-034_C1D1.rds"))
		SpaCET_obj <- SecAct.activity.inference.ST(inputProfile=SpaCET_obj, inputProfile_control=SpaCET_obj_CTL, sigMatrix=ref, lambda=5e05, is.group.sig=TRUE)
	}
	
	if(sampleName=="01-039_C3D1") 
	{
		SpaCET_obj_CTL <- readRDS(paste0(outputPath,"01-039_C1D1.rds"))
		SpaCET_obj <- SecAct.activity.inference.ST(inputProfile=SpaCET_obj, inputProfile_control=SpaCET_obj_CTL, sigMatrix=ref, lambda=5e05, is.group.sig=TRUE)
	}
}
saveRDS(SpaCET_obj, file = paste0(outputPath,sampleName,".rds"))


