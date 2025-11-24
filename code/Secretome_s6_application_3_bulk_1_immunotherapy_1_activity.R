source("Secretome_s0_path.R")

#SWARM -t 2 -g 30 --time 00:30:00
args = commandArgs(trailingOnly=TRUE)
cancer <- args[1]
datatype <- args[2]

inputPath <- paste0(dataAppPath,"Immunotherapy/")
outputPath <- paste0(applicationPath,"Immunotherapy/")
dir.create(outputPath, recursive = TRUE)

expr <- as.matrix(read.csv(gzfile(paste0(inputPath,cancer,".",datatype,".gz")),as.is=T,sep="\t",check.names=F))
expr <- expr[apply(expr,1,function(x) sum(x>0)>5),]

rownames(expr) <- transferSymbol(rownames(expr))
expr <- rm_duplicates(expr)
write.csv(expr, gzfile(paste0(outputPath,cancer,".csv.gz")),quote=F)	# for expr risk score


cdata_T_minusBG <- expr-rowMeans(expr)

library(SecAct)
res <- SecAct.activity.inference(inputProfile = cdata_T_minusBG, is.differential = TRUE)
save(res, file = paste0(outputPath,cancer,".RData"))
