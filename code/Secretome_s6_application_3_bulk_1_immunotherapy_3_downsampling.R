source("Secretome_s0_path.R")

#SWARM -t 10 -g 200 --time 5:30:00
args = commandArgs(trailingOnly=TRUE)
cancer <- args[1]
datatype <- args[2]

inputPath <- paste0(dataAppPath,"Immunotherapy/")
outputPath <- paste0(applicationPath,"Immunotherapy/downsampling/")
dir.create(outputPath)

expr <- as.matrix(read.csv(gzfile(paste0(inputPath,cancer,".",datatype,".gz")),as.is=T,sep="\t",check.names=F))
expr <- expr[apply(expr,1,function(x) sum(x>0)>5),]

rownames(expr) <- transferSymbol(rownames(expr))
expr <- rm_duplicates(expr)
#write.csv(expr, gzfile(paste0(outputPath,cancer,".csv.gz")),quote=F)	# for expr risk score


cdata_T_minusBG <- expr-rowMeans(expr)

library(SecAct)
#res <- SecAct.activity.inference(inputProfile = cdata_T_minusBG, is.differential = TRUE)
#save(res, file = paste0(outputPath,cancer,".RData"))


# down-sampling
if(!cancer%in%notGenomeWide)
{
	#ref <- file.path(system.file(package = "SecAct"), "extdata/SecAct.tsv.gz")
	#f1 <- read.table(ref)
	#
	#cdata_T_minusBG <- cdata_T_minusBG[rownames(cdata_T_minusBG)%in%rownames(f1),]
	
	for(n in n_downsampling)
	{
		dir.create(paste0(outputPath,"/",n))
		
		set.seed(123)
		#for(i in 1:n_downsampling_rep)
		runSecAct <- function(i)
		{
			if(n<nrow(cdata_T_minusBG))
			{
				geneSub <- sample(1:nrow(cdata_T_minusBG), n)
			}else{
				geneSub <- 1:nrow(cdata_T_minusBG)
			}
			
			cdata_T_minusBG_sub <- cdata_T_minusBG[geneSub,]
			res <- SecAct.activity.inference(inputProfile = cdata_T_minusBG_sub, is.differential = TRUE)
		
			save(res, file = paste0(outputPath,"/",n,"/rep_",i,"_",cancer,".RData"))
		}
		
		parallel::mclapply(1:n_downsampling_rep, runSecAct, mc.cores=10) 
	}

}
