source("Secretome_s0_path.R")

#SWARM -t 2 -g 20 --time 00:30:00
args = commandArgs(trailingOnly=TRUE)
sigName <- args[1]
cancer <- args[2]


inputPath <- paste0(dataValPath,"CPTAC/")
outputPath <- paste0(validationPath,"CPTAC/")


# Ensemble ID -> gene symbol
geneAnno <- read.csv(paste0(inputPath,"Proteome_BCM_GENCODE_v34_harmonized_v1/README/Gene_annotation_and_representable_isoform_mapping_table.txt"),sep="\t",check.names=F)
geneAnno[geneAnno[,1]=="ENSG00000284024.2",4] <- "MSANTD7"
geneAnno <- geneAnno[,c(1,4)]
geneAnno <- geneAnno[!duplicated(geneAnno[,1]),]
rownames(geneAnno) <- geneAnno[,1]


# read expression
cdata_T <- as.matrix(read.csv(paste0(inputPath,"RNA_BCM_v1/",cancer,"_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt"),sep="\t",row.names=1,check.names=F))
rownames(cdata_T) <- geneAnno[rownames(cdata_T),2]
rownames(cdata_T) <- transferSymbol(rownames(cdata_T))
cdata_T <- rm_duplicates(cdata_T)


if(file.exists(paste0(inputPath,"RNA_BCM_v1/",cancer,"_RNAseq_gene_RSEM_coding_UQ_1500_log2_Normal.txt")))
{
	cdata_N <- as.matrix(read.csv(paste0(inputPath,"RNA_BCM_v1/",cancer,"_RNAseq_gene_RSEM_coding_UQ_1500_log2_Normal.txt"),sep="\t",row.names=1,check.names=F))
	
	rownames(cdata_N) <- geneAnno[rownames(cdata_N),2]
	rownames(cdata_N) <- transferSymbol(rownames(cdata_N))
	cdata_N <- rm_duplicates(cdata_N)
	
	cdata_T <- cdata_T-rowMeans(cdata_N)
}else{
	cdata_T <- cdata_T-rowMeans(cdata_T)
}


# calculate secreted protein activity
if(sigName=="Pathway")
{
	library(GSVA)
	
	load(paste0(dataValPath,"Reactome/secreted_to_pathway_gmt.RData"))
	gsva_scores <- gsva(cdata_T, gmt, method = "gsva", kcdf = "Gaussian")

	save(gsva_scores, file = paste0(outputPath,"/",cancer,"_",sigName,".RData"))

}else if(sigName=="SecAct.LeaveOneOut"){

	library(SecAct)
	
	inputProfile <- cdata_T
	
	if(cancer=="COAD")
	{
		ref <- paste0(signatureCombPath,"AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_non-CRC_filterByPan_ds3_vst.tsv")
	}else{
		ref <- paste0(signatureCombPath,"AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_non-",cancer,"_filterByPan_ds3_vst.tsv")
	}
	
	is.group.sig <- TRUE
	lambda <- 5e+05

	res <- SecAct.activity.inference(
		inputProfile = inputProfile,
		is.differential = TRUE,
		sigMatrix = ref, 
		is.group.sig = is.group.sig,
		lambda = lambda,
		nrand = 1000
	)
	
	save(res, file = paste0(outputPath,"/",cancer,"_",sigName,".RData"))

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
