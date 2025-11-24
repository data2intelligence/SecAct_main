source("Secretome_s0_path.R")

#SWARM -t 8 -g 100 --time 03:00:00
args = commandArgs(trailingOnly=TRUE)
dataset_cancer <- args[1]

#for(dataset_cancer in dataset_cancers)
#{
	dataset <- strsplit(dataset_cancer,"_")[[1]][1]
	cancer <- strsplit(dataset_cancer,"_")[[1]][2]
	
	if(dataset=="TCGA")
	{
		cdata <- read.Xena(cancer)
		
		rownames(cdata) <- transferSymbol(rownames(cdata))
		cdata <- rm_duplicates(cdata)
		
		cdata <- filter.counts(cdata)
		
		cdata_T <- extract.samples(cdata,"T")
		cdata_N <- extract.samples(cdata,"N")
		
		if(ncol(cdata_N)>5)
		{
			cdata_T_minusBG <- cdata_T-rowMeans(cdata_N)
		}else{
			cdata_T_minusBG <- cdata_T-rowMeans(cdata_T)
		}
		
		cdata_T_minusBG <- as.matrix(cdata_T_minusBG)
	}else{
		cdata_T <- as.matrix(read.csv(paste0("../data/1_Development/ICGC/",cancer,".seq.expression.gz"),sep="\t",row.names=1))	
					
		rownames(cdata_T) <- transferSymbol(rownames(cdata_T))
		cdata_T <- rm_duplicates(cdata_T)
		
		cdata_T <- filter.counts(cdata_T)
		
		cdata_T_minusBG <- cdata_T-rowMeans(cdata_T)
		cdata_T_minusBG <- as.matrix(cdata_T_minusBG)
	}	
	
	
	library(SecAct)
	for(sigName in sigNames)
	{
		#for(lambda in lambdas)
		runRidge <- function(lambda)
		{
			ref <- paste0(finalSignaturesPath,sigName,".tsv.gz")
						
			res <- SecAct.activity.inference(
				inputProfile=cdata_T_minusBG, 
				is.differential=TRUE, 
				sigMatrix=ref,
				is.group.sig=FALSE,
				lambda=lambda
			)
			
			save(res, file = paste0(lambdaPath, dataset,"_",cancer,"_",sigName,"_",lambda,".RData"))
		}
		
		parallel::mclapply(lambdas, runRidge, mc.cores=7) 
	}
	
#}
