source("Secretome_s0_path.R")

#SWARM -t 2 -g 20 --time 00:30:00
args = commandArgs(trailingOnly=TRUE)
sigName <- args[1]
cancer <- args[2]


inputPath <- paste0(dataValPath,cancer,"/")
outputPath <- paste0(validationPath,cancer,"/")


# read proteomics
Protein <- read.csv(paste0(inputPath,"Proteomics_20250211/Protein_matrix_averaged_20250211.tsv"),check.names=F,sep="\t",comment.char = "#")
Protein_anno <- Protein[,1:2]
Protein <- Protein[,3:ncol(Protein)]
rownames(Protein) <- Protein_anno[,1]
Protein <- t(Protein)

Protein <- Protein[rowSums(!is.na(Protein))>ncol(Protein)*0.2,,drop=F]

rownames(Protein) <- transferSymbol(rownames(Protein))
Protein <- rm_duplicates(Protein)


# read expression
cdata_T <- read.csv(paste0(inputPath,"rnaseq_merged_20250117/rnaseq_merged_rsem_tpm_20250117.csv"),check.names=F,comment.char = "#")
cdata_T_anno <- cdata_T[,1:3]
cdata_T <- cdata_T[,4:ncol(cdata_T)]
cdata_T <- cdata_T[,colnames(cdata_T)%in%colnames(Protein)]


cdata_T_anno <- cdata_T_anno[apply(cdata_T,1,function(x) sum(x=="nan") )==0,]
cdata_T <- cdata_T[apply(cdata_T,1,function(x) sum(x=="nan") )==0,]

cdata_T_anno <- cdata_T_anno[apply(cdata_T,1,function(x) sum(x=="") )==0,]
cdata_T <- cdata_T[apply(cdata_T,1,function(x) sum(x=="") )==0,]

cdata_T_anno <- cdata_T_anno[apply(cdata_T,1,function(x) sum(grepl(" ",x)) )==0,]
cdata_T <- cdata_T[apply(cdata_T,1,function(x) sum(grepl(" ",x)) )==0,]


cdata_T <- apply(cdata_T, 2, as.numeric)
cdata_T <- as.matrix(cdata_T)


rownames(cdata_T) <- cdata_T_anno[,1]
rownames(cdata_T) <- transferSymbol(rownames(cdata_T))
cdata_T <- rm_duplicates(cdata_T)
	

# calculate secreted protein activity
if(sigName=="Pathway")
{
	library(GSVA)
	
	load(paste0(dataValPath,"Reactome/secreted_to_pathway_gmt.RData"))
	gsva_scores <- gsva(cdata_T, gmt, method = "gsva", kcdf = "Gaussian")

	save(gsva_scores, file = paste0(outputPath,"/",cancer,"_",sigName,".RData"))

}else{
	library(SecAct)

	inputProfile <- cdata_T
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
	
	save(res, file = paste0(outputPath,"/",cancer,"_",sigName,".RData"))
}
