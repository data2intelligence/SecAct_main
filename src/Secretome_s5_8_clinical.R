source("Secretome_s0_path.R")
source("../../SpaCE/code/ifun.R")

#SWARM -t 2 -g 50 --time 3:00:00
args = commandArgs(trailingOnly=TRUE)
item <- args[1]

#for(item in items)
#{
	gene <- strsplit(item,"_") [[1]][1]
	dataset <- strsplit(item,"_") [[1]][2]
	
	#if(item%in%c("TGFB1_GSE174686.12","TGFB1_GSE174686.pan","TGFB1_GSE174686.mousepan",
	#	"TGFB2_GSE174686.12","TGFB2_GSE174686.pan","TGFB2_GSE174686.mousepan",
	#	"TGFB3_GSE174686","TGFB3_GSE174686.pan","TGFB3_GSE174686.mousepan",
	#	"TGFB_GSE174686.12","TGFB_GSE174686.pan","TGFB_GSE174686.mousepan"
	#	))
	#{
	#	f3 <- as.matrix(read.csv(gzfile("../../SpaCE/data/aTGFB_XOMA.diff.human.gz"),as.is=T,sep="\t",check.names=F))
	#	
	#	if(item%in%c("TGFB1_GSE174686.12","TGFB2_GSE174686.12","TGFB_GSE174686.12"))
	#	{	
	#		f3 <- f3[,"1+2",drop=F]
	#	}
	#	if(item%in%c("TGFB1_GSE174686.pan","TGFB2_GSE174686.pan","TGFB3_GSE174686.pan","TGFB_GSE174686.pan"))
	#	{	
	#		f3 <- f3[,"Pan",drop=F]
	#	}
	#	if(item%in%c("TGFB1_GSE174686.mousepan","TGFB2_GSE174686.mousepan","TGFB3_GSE174686.mousepan","TGFB_GSE174686.mousepan"))
	#	{	
	#		f3 <- f3[,"mouse.Pan",drop=F]
	#	}
	#	if(item=="TGFB3_GSE174686")
	#	{	
	#		f3 <- f3[,"Pan-1+2",drop=F]
	#	}
	if(item%in%c("TGFB_GSE174686","TGFB3_GSE174686"))
	{
		dataPath <- "../../SpaCE/data/"
		f3 <- as.matrix(read.csv(gzfile(paste0(dataPath,"aTGFB_XOMA.diff.human.gz")),as.is=T,sep="\t",check.names=F))
		
		if(item%in%c("TGFB_GSE174686"))
		{	
			f3 <- f3[,"mouse.Pan",drop=F]
		}
		if(item=="TGFB3_GSE174686")
		{	
			f3 <- f3[,"Pan-1+2",drop=F]
		}
	}else if(item%in%c("VEGFA_GSE72951","VEGFA_E-MTAB-3267")){
		dataPath <- "../../SpaCE/code/CytoSig_prediction-master/data/bulk/tumor/"
		
		if(item=="VEGFA_GSE72951")
		{
			f3 <- as.matrix(read.csv(gzfile(paste0(dataPath,"GSE72951.self_subtract.gz")),as.is=T,sep="\t"))
		}else if(item=="VEGFA_E-MTAB-3267"){
			f3 <- as.matrix(read.csv(gzfile(paste0(dataPath,"E-MTAB-3267.norm_subtract.gz")),as.is=T,sep="\t",check.names=F))
		}
	
	}else if(item%in%c("NTN1_GSE225691")){
		
		NTN_bulk <- as.matrix(read.csv("/data/rub2/project/Secretome/data/GSE225691_NTN1/GSE225687_Patient_Raw_Count.csv",sep=";",row.names=1,header=T))
		rownames(NTN_bulk) <- substr(rownames(NTN_bulk),1,15)
		NTN_bulk <- rm_duplicates(NTN_bulk)
		
		gene_anno <- as.matrix(read.csv(gzfile(paste0("../_raw/BRCA_10x_Datasets/Version1.0.0_Breast.Cancer_rep1/filtered_feature_bc_matrix/features.tsv.gz")),as.is=T,header=F,sep="\t"))
		
		olp <- intersect(rownames(NTN_bulk),gene_anno[,1])
		
		NTN_bulk <- NTN_bulk[olp,]
		rownames(NTN_bulk) <- gene_anno[match(rownames(NTN_bulk),gene_anno[,1]),2]
		rownames(NTN_bulk) <- transferSymbol(rownames(NTN_bulk))
		NTN_bulk <- rm_duplicates(NTN_bulk)
		
		expr <- NTN_bulk
		
		expr.scaled <- t(t(expr)*1e6/colSums(expr))
		expr.scaled.log <- log2(expr.scaled + 1 )
				
		f3 <- expr.scaled.log[,(1:12)*2] - expr.scaled.log[,(1:12)*2 - 1]
	
	}else if(item%in%c("IL11_Widjaja2024")){
		
		library(rio) 
		xlsx <- import_list(paste0(validationPath,"/IL11/IL11_Widjaja2024_TableS5.xlsx")) 
		
		s0 <- as.data.frame(xlsx[[1]]) 
		s0 <- s0[s0[,2]>0,]
		s0[,3] <- as.numeric(s0[,3])
		s0 <- s0[order(s0[,2],decreasing=T),]
		s0 <- s0[!duplicated(s0[,8]),]
		rownames(s0) <- s0[,8]
		s1 <- s0[,3,drop=F]
		
		s0 <- as.data.frame(xlsx[[2]]) 
		s0 <- s0[s0[,2]>0,]
		s0[,3] <- as.numeric(s0[,3])
		s0 <- s0[order(s0[,2],decreasing=T),]
		s0 <- s0[!duplicated(s0[,8]),]
		rownames(s0) <- s0[,8]
		s2 <- s0[,3,drop=F]
		
		s0 <- as.data.frame(xlsx[[3]]) 
		s0 <- s0[s0[,2]>0,]
		s0[,3] <- as.numeric(s0[,3])
		s0 <- s0[order(s0[,2],decreasing=T),]
		s0 <- s0[!duplicated(s0[,8]),]
		rownames(s0) <- s0[,8]
		s3 <- s0[,3,drop=F]
		
		olp <- intersect(rownames(s1), intersect(rownames(s2),rownames(s3)))
		
		f3 <- cbind(vWat=s1[olp,], cbind(liver=s2[olp,], gastro=s3[olp,]) )
		rownames(f3) <- olp
		
		f3 <- cbind(f3, meanLogFC=apply(f3,1,mean))
		
		f3 <- transferMouseToHuman(f3)
		
		f3 <- f3[,4,drop=F]
	
	}else if(item%in%c("IL11_Widjaja2024_vWat")){
		
		library(rio) 
		xlsx <- import_list(paste0(validationPath,"/IL11/IL11_Widjaja2024_TableS5.xlsx")) 
		
		s0 <- as.data.frame(xlsx[[1]]) 
		s0 <- s0[s0[,2]>0,]
		s0[,3] <- as.numeric(s0[,3])
		s0 <- s0[order(s0[,2],decreasing=T),]
		s0 <- s0[!duplicated(s0[,8]),]
		rownames(s0) <- s0[,8]
		s1 <- s0[,3,drop=F]
		
		s0 <- as.data.frame(xlsx[[2]]) 
		s0 <- s0[s0[,2]>0,]
		s0[,3] <- as.numeric(s0[,3])
		s0 <- s0[order(s0[,2],decreasing=T),]
		s0 <- s0[!duplicated(s0[,8]),]
		rownames(s0) <- s0[,8]
		s2 <- s0[,3,drop=F]
		
		s0 <- as.data.frame(xlsx[[3]]) 
		s0 <- s0[s0[,2]>0,]
		s0[,3] <- as.numeric(s0[,3])
		s0 <- s0[order(s0[,2],decreasing=T),]
		s0 <- s0[!duplicated(s0[,8]),]
		rownames(s0) <- s0[,8]
		s3 <- s0[,3,drop=F]
		
		olp <- intersect(rownames(s1), intersect(rownames(s2),rownames(s3)))
		
		f3 <- cbind(vWat=s1[olp,], cbind(liver=s2[olp,], gastro=s3[olp,]) )
		rownames(f3) <- olp
		
		f3 <- cbind(f3, meanLogFC=apply(f3,1,mean))
		
		f3 <- transferMouseToHuman(f3)
		
		f3 <- f3[,1,drop=F]
		
	}else if(item%in%c("IL11_Widjaja2024_liver")){
		
		library(rio) 
		xlsx <- import_list(paste0(validationPath,"/IL11/IL11_Widjaja2024_TableS5.xlsx")) 
		
		s0 <- as.data.frame(xlsx[[1]]) 
		s0 <- s0[s0[,2]>0,]
		s0[,3] <- as.numeric(s0[,3])
		s0 <- s0[order(s0[,2],decreasing=T),]
		s0 <- s0[!duplicated(s0[,8]),]
		rownames(s0) <- s0[,8]
		s1 <- s0[,3,drop=F]
		
		s0 <- as.data.frame(xlsx[[2]]) 
		s0 <- s0[s0[,2]>0,]
		s0[,3] <- as.numeric(s0[,3])
		s0 <- s0[order(s0[,2],decreasing=T),]
		s0 <- s0[!duplicated(s0[,8]),]
		rownames(s0) <- s0[,8]
		s2 <- s0[,3,drop=F]
		
		s0 <- as.data.frame(xlsx[[3]]) 
		s0 <- s0[s0[,2]>0,]
		s0[,3] <- as.numeric(s0[,3])
		s0 <- s0[order(s0[,2],decreasing=T),]
		s0 <- s0[!duplicated(s0[,8]),]
		rownames(s0) <- s0[,8]
		s3 <- s0[,3,drop=F]
		
		olp <- intersect(rownames(s1), intersect(rownames(s2),rownames(s3)))
		
		f3 <- cbind(vWat=s1[olp,], cbind(liver=s2[olp,], gastro=s3[olp,]) )
		rownames(f3) <- olp
		
		f3 <- cbind(f3, meanLogFC=apply(f3,1,mean))
		
		f3 <- transferMouseToHuman(f3)
		
		f3 <- f3[,2,drop=F]
	
	}else if(item%in%c("IL11_Widjaja2024_gastro")){
		
		library(rio) 
		xlsx <- import_list(paste0(validationPath,"/IL11/IL11_Widjaja2024_TableS5.xlsx")) 
		
		s0 <- as.data.frame(xlsx[[1]]) 
		s0 <- s0[s0[,2]>0,]
		s0[,3] <- as.numeric(s0[,3])
		s0 <- s0[order(s0[,2],decreasing=T),]
		s0 <- s0[!duplicated(s0[,8]),]
		rownames(s0) <- s0[,8]
		s1 <- s0[,3,drop=F]
		
		s0 <- as.data.frame(xlsx[[2]]) 
		s0 <- s0[s0[,2]>0,]
		s0[,3] <- as.numeric(s0[,3])
		s0 <- s0[order(s0[,2],decreasing=T),]
		s0 <- s0[!duplicated(s0[,8]),]
		rownames(s0) <- s0[,8]
		s2 <- s0[,3,drop=F]
		
		s0 <- as.data.frame(xlsx[[3]]) 
		s0 <- s0[s0[,2]>0,]
		s0[,3] <- as.numeric(s0[,3])
		s0 <- s0[order(s0[,2],decreasing=T),]
		s0 <- s0[!duplicated(s0[,8]),]
		rownames(s0) <- s0[,8]
		s3 <- s0[,3,drop=F]
		
		olp <- intersect(rownames(s1), intersect(rownames(s2),rownames(s3)))
		
		f3 <- cbind(vWat=s1[olp,], cbind(liver=s2[olp,], gastro=s3[olp,]) )
		rownames(f3) <- olp
		
		f3 <- cbind(f3, meanLogFC=apply(f3,1,mean))
		
		f3 <- transferMouseToHuman(f3)
		
		f3 <- f3[,3,drop=F]
		
	}else{
		dataPath <- "../../SpaCE/code/CytoSig_prediction-master/data/bulk/inflam/"
		allFiles <- list.files(dataPath)
		
		flag1 <- grepl(dataset,allFiles)
		flag2 <- grepl(".diff.1.cntmap",allFiles)
		
		fileName_full <- allFiles[flag1&flag2]
		fileName_full <- gsub(".diff.1.cntmap","",fileName_full)
			
		f3 <- as.matrix(read.csv(gzfile(paste0(dataPath,fileName_full,".diff.1")),as.is=T,sep="\t"))

	}
	
	rownames(f3) <- transferSymbol(rownames(f3))
	f3 <- rm_duplicates(f3)
	
	library(SecAct)
	#for(a in c(10000,50000,100000,500000,1000000,5000000,10000000,50000000))
	#{
		for(sig in c(100011,100012,100021,100022,100031,100032,100041,100042,100051,100052))
		#for(sig in c(61811,61822,61833,61844,61855,61866,61877))
		#for(sig in c(0,11,112,36,618))
		#for(sig in c(7777))
		{
			if(sig==0)
			{
				ref <- paste0("/data/rub2/project/Secretome/code/CytoSig-master/CytoSig/signature.centroid")
			}else if(sig==618){
				ref <- paste0("/data/rub2/project/Secretome/results/signatureComb/AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3.tsv")
			}else if(sig==6182){
				ref <- paste0("/data/rub2/project/Secretome/results/signatureComb/AllSigFilteredBy_MoranI_TCGA_0.25_ds3_ICGC_highConfidence.tsv")
			}else if(sig==61811){
				ref <- paste0("/data/rub2/project/Secretome/results/signatureComb/AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3_weighted1.tsv")
			}else if(sig==61822){
				ref <- paste0("/data/rub2/project/Secretome/results/signatureComb/AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3_weighted2.tsv")
			}else if(sig==61833){
				ref <- paste0("/data/rub2/project/Secretome/results/signatureComb/SecAct.signature.adjusted_v1")
			}else if(sig==61844){
				ref <- paste0("/data/rub2/project/Secretome/results/signatureComb/SecAct.signature.adjusted_v2")
			}else if(sig==61855){
				ref <- paste0("/data/rub2/project/Secretome/results/signatureComb/SecAct.signature.adjusted_v2_0.1")
			}else if(sig==61866){
				ref <- paste0("/data/rub2/project/Secretome/results/signatureComb/SecAct.signature.adjusted_v2_0.2")
			}else if(sig==61877){
				ref <- paste0("/data/rub2/project/Secretome/results/signatureComb/AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3_weighted2_0.1")
			}else if(sig==61888){
				ref <- paste0("/data/rub2/project/Secretome/results/signatureComb/AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3_1k_0.1")
			}else if(sig==61899){
				ref <- paste0("/data/rub2/project/Secretome/results/signatureComb/AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3_1k_residual_0.1.tsv")
			}else if(sig==61800){
				ref <- paste0("/data/rub2/project/Secretome/results/signatureComb/AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3_1k_weighted2.tsv")
			}else if(sig==100011){
				ref <- paste0(signatureCombPath,"AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3_1k_vst.tsv")
			}else if(sig==100012){
				ref <- paste0(signatureCombPath,"AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3_1k_vst_0.9.tsv")
			}else if(sig==100021){
				ref <- paste0(signatureCombPath,"AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3_1k_weighted.tsv")
			}else if(sig==100022){
				ref <- paste0(signatureCombPath,"AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3_1k_weighted_0.9.tsv")
			}else if(sig==100031){
				ref <- paste0(signatureCombPath,"AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3_1k_weighted2.tsv")
			}else if(sig==100032){
				ref <- paste0(signatureCombPath,"AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3_1k_weighted2_0.9.tsv")
			}else if(sig==100041){
				ref <- paste0(signatureCombPath,"AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3_1k_vst_condition_logUMI_cellType.tsv")
			}else if(sig==100042){
				ref <- paste0(signatureCombPath,"AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3_1k_vst_condition_logUMI_cellType_0.9.tsv")
			}else if(sig==100051){
				ref <- paste0(signatureCombPath,"AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3_1k_vst_condition_cellType.tsv")
			}else if(sig==100052){
				ref <- paste0(signatureCombPath,"AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3_1k_vst_condition_cellType_0.9.tsv")
			}else if(sig==7777){
				ref <- paste0("/data/rub2/project/Secretome/results/signatureComb/AllSigFilteredBy_Membrane_MoranI_TCGA_ICGC_0.25_ds3.tsv")
			}else{
				ref <- paste0("/data/rub2/project/Secretome/code/CytoSig-master/CytoSig/signature.centroid.",sig)
			}
			
			#aMat <- read.csv(paste0(QCSummaryPath,"lambda_summary/lambda_compare_",sig,".csv"))
			#a <- aMat[aMat[,2]==max(aMat[,2]),1]
			a <- 5e+05
			
			res <- SecAct.inference(f3, SigMat=ref, lambda=a, nrand=1000)
			save(res, file = paste0(validationPath,"clinical/",item,"_",sig,"_",a,".RData"))
			
			print(paste0(a,"_",sig))
			
		}
	#}
#}




if(FALSE)
{


items <- c(
#"TGFB1_GSE174686.12","TGFB1_GSE174686.pan","TGFB1_GSE174686.mousepan",
#"TGFB2_GSE174686.12","TGFB2_GSE174686.pan","TGFB2_GSE174686.mousepan",
#"TGFB3_GSE174686","TGFB3_GSE174686.pan","TGFB3_GSE174686.mousepan",
#"TGFB_GSE174686.12","TGFB_GSE174686.pan","TGFB_GSE174686.mousepan",
"TGFB_GSE174686","TGFB3_GSE174686",
"IL1B_GSE80060","IL1B_GSE57253","IL1B_GSE70019",
#"IL17A_GSE31652","IL17A_GSE11903","IL17A_GSE55201",
#"IL22_GSE99802",
"IFNG_GSE100093","IFNG_GSE78193",
"TNFSF12_GSE42049","TNFSF12_GSE42048", # tumor
"IL6_GSE62941","IL6_GSE45867","IL6_GSE61201", # tumor x x 
"TNF_GSE48498","TNF_GSE11903",
"IL1A_GSE57253","IL1A_GSE70019",
#"IL4_GSE130588","IL4_GSE59294",
#"IL13_GSE130588","IL13_GSE59294",
"NTN1_GSE225691", # tumor
#"IL11_Widjaja2024_vWat",
#"IL11_Widjaja2024_liver",
#"IL11_Widjaja2024_gastro",
"IL11_Widjaja2024"

)


compSigs <- c(
	0,11,"11_abs",112,"112_abs",
	36,37,618,6182)
compSigNames <- c(
	"CytoSig","NicheNet1","NicheNet1_abs","NicheNet","NicheNet_abs",
	"ImmuneDicFilter","SecAct398","SecAct","SecActv2")


compSigs <- c(
	0,11,112,
	36,618,61811,61822,61833,61844,61855,61866,61877,61888,61899,61800,
	100011,100012,100021,100022,100031,100032,100041,100042,100051,100052
	#,7777
	)
compSigNames <- c(
	"CytoSig","NicheNet.v1","NicheNet.v2",
	"ImmuneDic","SecAct","SecAct.w1","SecAct.w2","SecAct.v1","SecAct.v2","SecAct.v2_0.9","SecAct.v2_0.8","SecAct.w2_0.9","SecAct.1k_0.9","SecAct.1k_residual","SecAct.1k",
	"SecAct.1k_11","SecAct.1k_12","SecAct.1k_21","SecAct.1k_22","SecAct.1k_31","SecAct.1k_32","SecAct.1k_41","SecAct.1k_42","SecAct.1k_51","SecAct.1k_52"
	#,"MemAct"
	)
CompSigAlphas <- c(
	1e+07,1e+07,1e+07,
	50000,5e+05,5e+05,5e+05,5e+05,5e+05,5e+05,5e+05,5e+05,5e+05,1e+05,5e+05,
	5e+05,5e+05,5e+05,5e+05,5e+05,5e+05,5e+05,5e+05,5e+05,5e+05
	#,5e+06
)


compSigs <- c(
	0,11,112,
	36,100042
	)
compSigNames <- c(
	"CytoSig","NicheNet.v1","NicheNet.v2",
	"ImmuneDic","SecAct"
	)
CompSigAlphas <- c(
	1e+07,1e+07,1e+07,
	50000,5e+05
)


source("createLRdb.R")
LRdb <- rbind(LRdb,c("NTN1","DCC"))
LRdb <- rbind(LRdb,c("NTN1","UNC5C"))


for(compSig in compSigs)
{	
	clinical <- data.frame()
	
	for(item in items)
	{
		gene <- strsplit(item,"_") [[1]][1]
		dataset <- strsplit(item,"_") [[1]][2]
	
		#for(a in c(10000,50000,1e+05,5e+05,1e+06,5e+06,1e+07,5e+07) )
		#{
			#if(!file.exists(paste0(validationPath,"clinical/",item,"_",compSig,"_",a,".RData"))) next
			a <- CompSigAlphas[compSigs==compSig]

			load(paste0(validationPath,"clinical/",item,"_",compSig,"_",a,".RData"))
			f3 <- as.matrix(res$zscore)
			f3 <- expand_rows(f3)
			
			if(compSig==0)
			{
				rownames(f3) <- transferSymbol(rownames(f3))
				f3 <- rm_duplicates(f3)
			}
			
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
			
			clinical[paste0(item,"_",compSig,"_",a),"item"] <- item
			clinical[paste0(item,"_",compSig,"_",a),"compSig"] <- paste0(compSig,"_",a)
			
			if(compSig==7777)
			{
				if(gene=="TGFB")
				{
					if(length(unique(LRdb[LRdb[,1]%in%c("TGFB1","TGFB2","TGFB3"),2]))==0)
					{
						clinical[paste0(item,"_",compSig,"_",a),"value"] <- NA
					}else{
						clinical[paste0(item,"_",compSig,"_",a),"value"] <- mean( f3[unique(LRdb[LRdb[,1]%in%c("TGFB1","TGFB2","TGFB3"),2])], na.rm=T)
					}
				}else{
					if(length(unique(LRdb[LRdb[,1]%in%gene,2]))==0)
					{
						clinical[paste0(item,"_",compSig,"_",a),"value"] <- NA
					}else{
						
						clinical[paste0(item,"_",compSig,"_",a),"value"] <- mean( f3[unique(LRdb[LRdb[,1]%in%gene,2])], na.rm=T)
					}
				}
			}else{
				clinical[paste0(item,"_",compSig,"_",a),"value"] <- f3[gene]
			}

		#} #a
	}
	
	
	mat = reshape2::dcast( clinical , compSig~item )
	rownames(mat) <- mat[,1]
	mat <- t(mat[,-1])
	
	mat <- mat[order(apply(mat,1,function(x) median(x,na.rm=T))),,drop=F]
	mat <- mat[,order(apply(mat,2,function(x) median(x,na.rm=T))),drop=F]
	
	
	#library(ComplexHeatmap)
	#
	#png(paste0(validationPath,"validation_clinical_blockade_",compSigNames[which(compSigs==compSig)],".png"), width = 20, height = 20, res=200, units = "cm")
	#
	#row_ha <- rowAnnotation(
	#	z = anno_boxplot(as.matrix(mat), height = unit(3, "cm") ) 
	#)
	#	column_ha <- columnAnnotation(
	#	z = anno_boxplot(as.matrix(mat), height = unit(3, "cm") ) 
	#)
	#		
	#ht <- Heatmap(as.matrix(mat),
	#	column_title=compSigNames[which(compSigs==compSig)],
	#	name = "Risk z",
	#	col = circlize::colorRamp2(c(-3, 0,3), c("green", "white", "red")),
	#	row_names_max_width = max_text_width(
	#        rownames(mat), 
	#        gp = gpar(fontsize = 12)
	#        ),
	#    column_names_max_height = max_text_width(
	#        colnames(mat), 
	#        gp = gpar(fontsize = 12)
	#        ),
	#    top_annotation = column_ha,
	#    right_annotation = row_ha,
  	#	cluster_rows = FALSE,
	#	cluster_columns = FALSE,
	#	cell_fun = function(j, i, x, y, width, height, fill) {grid.text(sprintf("%.3f", mat[i, j]), x, y, gp = gpar(fontsize = 10))}
	#)
	#
	#draw(ht)
	#dev.off()
	
	#clinicalBest <- clinical[clinical[,2]==colnames(mat)[1],,drop=F]
	clinicalBest <- clinical
	
	if(compSig==0)
	{
		clinicalComb <- clinicalBest
	}else{
		clinicalComb <- rbind(clinicalComb,clinicalBest)
	}
}


for(compSig in compSigs)
{	
	clinicalComb[clinicalComb[,2]==paste0(compSig,"_",CompSigAlphas[compSigs==compSig]),2] <- compSigNames[compSigs==compSig]	
}


#clinical[clinical[,2]=="0_1e+06",2] <- "CytoSig"		
#clinical[clinical[,2]=="11_1e+06",2] <- "NicheNetSig"		
#clinical[clinical[,2]=="12_1e+06",2] <- "SecAct"		
#clinical[clinical[,2]=="13_1e+06",2] <- "SecAct_normalSPs_all"		
#clinical[clinical[,2]=="14_1e+06",2] <- "SecAct_normalSPs_MoranI"		
#clinical[clinical[,2]=="15_1e+06",2] <- "SecAct_normalSPs_Both"		
#clinical[clinical[,2]=="16_1e+06",2] <- "SecAct_normalSPs_z_Both"		
#clinical[clinical[,2]=="17_1e+06",2] <- "SecAct_normalFilterSPs_MoranI"		
#clinical[clinical[,2]=="18_1e+06",2] <- "SecAct_normalFilterSPs_Both"		
#		
#clinical <- clinical[clinical[,2]%in%c("CytoSig","SecAct","NicheNetSig","SecAct_normalSPs_all",
#"SecAct_normalSPs_MoranI","SecAct_normalSPs_Both","SecAct_normalSPs_z_Both",
#"SecAct_normalFilterSPs_MoranI","SecAct_normalFilterSPs_Both"
#),]

#clinicalComb[clinicalComb[,2]=="0_5e+07",2] <- "CytoSig"		
#clinicalComb[clinicalComb[,2]=="11_5e+07",2] <- "NicheNet"		
#clinicalComb[clinicalComb[,2]=="11_abs_1e+07",2] <- "NicheNet_abs"		
#clinicalComb[clinicalComb[,2]=="112_1e+07",2] <- "NicheNet2"		
#clinicalComb[clinicalComb[,2]=="112_abs_1e+07",2] <- "NicheNet2_abs"	
#clinicalComb[clinicalComb[,2]=="36_50000",2] <- "ImmuneDic"
#clinicalComb[clinicalComb[,2]=="37_5e+05",2] <- "SecAct398"
#clinicalComb[clinicalComb[,2]=="618_5e+05",2] <- "SecAct"
#clinicalComb[clinicalComb[,2]=="6182_5e+05",2] <- "SecActv2"
#clinical[clinical[,2]=="20_1e+06",2] <- "ImmuneDicSig"
#clinical[clinical[,2]=="21_1e+06",2] <- "ImmuneDicSig.B_cell"
#clinical[clinical[,2]=="22_1e+06",2] <- "ImmuneDicSig.cDC1"
#clinical[clinical[,2]=="23_1e+06",2] <- "ImmuneDicSig.cDC2"
#clinical[clinical[,2]=="24_1e+06",2] <- "ImmuneDicSig.eTAC"
#clinical[clinical[,2]=="25_1e+06",2] <- "ImmuneDicSig.ILC"
#clinical[clinical[,2]=="26_1e+06",2] <- "ImmuneDicSig.Macrophage"
#clinical[clinical[,2]=="27_1e+06",2] <- "ImmuneDicSig.MigDC"
#clinical[clinical[,2]=="28_1e+06",2] <- "ImmuneDicSig.Monocyte"
#clinical[clinical[,2]=="29_1e+06",2] <- "ImmuneDicSig.Neutrophil"
#clinical[clinical[,2]=="30_1e+06",2] <- "ImmuneDicSig.NK_cell"
#clinical[clinical[,2]=="31_1e+06",2] <- "ImmuneDicSig.pDC"
#clinical[clinical[,2]=="32_1e+06",2] <- "ImmuneDicSig.T_cell_CD4"
#clinical[clinical[,2]=="33_1e+06",2] <- "ImmuneDicSig.T_cell_CD8"
#clinical[clinical[,2]=="34_1e+06",2] <- "ImmuneDicSig.T_cell_gd"
#clinical[clinical[,2]=="35_1e+06",2] <- "ImmuneDicSig.Treg"
#clinical[clinical[,2]=="36_1e+06",2] <- "ImmuneDicQCSig"



for(item in items)
{
	gene <- strsplit(item,"_") [[1]][1]
	dataset <- strsplit(item,"_") [[1]][2]
	
	if(item%in%c("TGFB_GSE174686","TGFB3_GSE174686"))
	{
		f3 <- as.matrix(read.csv(gzfile("../../SpaCE/data/aTGFB_XOMA.diff.human.gz"),as.is=T,sep="\t",check.names=F))
		rownames(f3) <- transferSymbol(rownames(f3))
		f3 <- rm_duplicates(f3)
	
		if(item%in%c("TGFB_GSE174686"))
		{	
			f3 <- f3[,"mouse.Pan",drop=F]
			f3 <- rbind(f3,TGFB=colMeans(f3[rownames(f3)%in%c("TGFB1","TGFB2","TGFB3"),,drop=F]))
		}
		if(item=="TGFB3_GSE174686")
		{	
			f3 <- f3[,"Pan-1+2",drop=F]
		}
		
	}else if(item%in%c("NTN1_GSE225691")){	
		gene_anno <- as.matrix(read.csv(gzfile(paste0("../_raw/BRCA_10x_Datasets/Version1.0.0_Breast.Cancer_rep1/filtered_feature_bc_matrix/features.tsv.gz")),as.is=T,header=F,sep="\t"))

		NTN_bulk <- as.matrix(read.csv("/data/rub2/project/Secretome/data/GSE225691_NTN1/GSE225687_Patient_Raw_Count.csv",sep=";",row.names=1,header=T))
		rownames(NTN_bulk) <- substr(rownames(NTN_bulk),1,15)
		NTN_bulk <- rm_duplicates(NTN_bulk)
		
		olp <- intersect(rownames(NTN_bulk),gene_anno[,1])
		
		NTN_bulk <- NTN_bulk[olp,]
		rownames(NTN_bulk) <- gene_anno[match(rownames(NTN_bulk),gene_anno[,1]),2]
		rownames(NTN_bulk) <- transferSymbol(rownames(NTN_bulk))
		NTN_bulk <- rm_duplicates(NTN_bulk)
		
		expr <- NTN_bulk
		
		expr.scaled <- t(t(expr)*1e6/colSums(expr))
		expr.scaled.log <- log2(expr.scaled + 1 )
				
		
		f3 <- expr.scaled.log[,(1:12)*2] - expr.scaled.log[,(1:12)*2 - 1]
		
	}else if(item%in%c("IL11_Widjaja2024")){
		
		library(rio) 
		xlsx <- import_list(paste0(validationPath,"/IL11/IL11_Widjaja2024_TableS5.xlsx")) 
		
		s0 <- as.data.frame(xlsx[[1]]) 
		s0 <- s0[s0[,2]>0,]
		s0[,3] <- as.numeric(s0[,3])
		s0 <- s0[order(s0[,2],decreasing=T),]
		s0 <- s0[!duplicated(s0[,8]),]
		rownames(s0) <- s0[,8]
		s1 <- s0[,3,drop=F]
		
		s0 <- as.data.frame(xlsx[[2]]) 
		s0 <- s0[s0[,2]>0,]
		s0[,3] <- as.numeric(s0[,3])
		s0 <- s0[order(s0[,2],decreasing=T),]
		s0 <- s0[!duplicated(s0[,8]),]
		rownames(s0) <- s0[,8]
		s2 <- s0[,3,drop=F]
		
		s0 <- as.data.frame(xlsx[[3]]) 
		s0 <- s0[s0[,2]>0,]
		s0[,3] <- as.numeric(s0[,3])
		s0 <- s0[order(s0[,2],decreasing=T),]
		s0 <- s0[!duplicated(s0[,8]),]
		rownames(s0) <- s0[,8]
		s3 <- s0[,3,drop=F]
		
		olp <- intersect(rownames(s1), intersect(rownames(s2),rownames(s3)))
		
		f3 <- cbind(s1[olp,], cbind(s2[olp,], s3[olp,]) )
		rownames(f3) <- olp
		
		f3 <- cbind(f3, meanLogFC=apply(f3,1,mean))
		
		f3 <- transferMouseToHuman(f3)
		
		f3 <- f3[,4,drop=F]

		rownames(f3) <- transferSymbol(rownames(f3))
		f3 <- rm_duplicates(f3)
	}else if(item%in%c("IL11_Widjaja2024_vWat")){
		
		library(rio) 
		xlsx <- import_list(paste0(validationPath,"/IL11/IL11_Widjaja2024_TableS5.xlsx")) 
		
		s0 <- as.data.frame(xlsx[[1]]) 
		s0 <- s0[s0[,2]>0,]
		s0[,3] <- as.numeric(s0[,3])
		s0 <- s0[order(s0[,2],decreasing=T),]
		s0 <- s0[!duplicated(s0[,8]),]
		rownames(s0) <- s0[,8]
		s1 <- s0[,3,drop=F]
		
		s0 <- as.data.frame(xlsx[[2]]) 
		s0 <- s0[s0[,2]>0,]
		s0[,3] <- as.numeric(s0[,3])
		s0 <- s0[order(s0[,2],decreasing=T),]
		s0 <- s0[!duplicated(s0[,8]),]
		rownames(s0) <- s0[,8]
		s2 <- s0[,3,drop=F]
		
		s0 <- as.data.frame(xlsx[[3]]) 
		s0 <- s0[s0[,2]>0,]
		s0[,3] <- as.numeric(s0[,3])
		s0 <- s0[order(s0[,2],decreasing=T),]
		s0 <- s0[!duplicated(s0[,8]),]
		rownames(s0) <- s0[,8]
		s3 <- s0[,3,drop=F]
		
		olp <- intersect(rownames(s1), intersect(rownames(s2),rownames(s3)))
		
		f3 <- cbind(vWat=s1[olp,], cbind(liver=s2[olp,], gastro=s3[olp,]) )
		rownames(f3) <- olp
		
		f3 <- cbind(f3, meanLogFC=apply(f3,1,mean))
		
		f3 <- transferMouseToHuman(f3)
		
		f3 <- f3[,1,drop=F]
		rownames(f3) <- transferSymbol(rownames(f3))
		f3 <- rm_duplicates(f3)
	}else if(item%in%c("IL11_Widjaja2024_liver")){
		
		library(rio) 
		xlsx <- import_list(paste0(validationPath,"/IL11/IL11_Widjaja2024_TableS5.xlsx")) 
		
		s0 <- as.data.frame(xlsx[[1]]) 
		s0 <- s0[s0[,2]>0,]
		s0[,3] <- as.numeric(s0[,3])
		s0 <- s0[order(s0[,2],decreasing=T),]
		s0 <- s0[!duplicated(s0[,8]),]
		rownames(s0) <- s0[,8]
		s1 <- s0[,3,drop=F]
		
		s0 <- as.data.frame(xlsx[[2]]) 
		s0 <- s0[s0[,2]>0,]
		s0[,3] <- as.numeric(s0[,3])
		s0 <- s0[order(s0[,2],decreasing=T),]
		s0 <- s0[!duplicated(s0[,8]),]
		rownames(s0) <- s0[,8]
		s2 <- s0[,3,drop=F]
		
		s0 <- as.data.frame(xlsx[[3]]) 
		s0 <- s0[s0[,2]>0,]
		s0[,3] <- as.numeric(s0[,3])
		s0 <- s0[order(s0[,2],decreasing=T),]
		s0 <- s0[!duplicated(s0[,8]),]
		rownames(s0) <- s0[,8]
		s3 <- s0[,3,drop=F]
		
		olp <- intersect(rownames(s1), intersect(rownames(s2),rownames(s3)))
		
		f3 <- cbind(vWat=s1[olp,], cbind(liver=s2[olp,], gastro=s3[olp,]) )
		rownames(f3) <- olp
		
		f3 <- cbind(f3, meanLogFC=apply(f3,1,mean))
		
		f3 <- transferMouseToHuman(f3)
		
		f3 <- f3[,2,drop=F]
		rownames(f3) <- transferSymbol(rownames(f3))
		f3 <- rm_duplicates(f3)
	}else if(item%in%c("IL11_Widjaja2024_gastro")){
		
		library(rio) 
		xlsx <- import_list(paste0(validationPath,"/IL11/IL11_Widjaja2024_TableS5.xlsx")) 
		
		s0 <- as.data.frame(xlsx[[1]]) 
		s0 <- s0[s0[,2]>0,]
		s0[,3] <- as.numeric(s0[,3])
		s0 <- s0[order(s0[,2],decreasing=T),]
		s0 <- s0[!duplicated(s0[,8]),]
		rownames(s0) <- s0[,8]
		s1 <- s0[,3,drop=F]
		
		s0 <- as.data.frame(xlsx[[2]]) 
		s0 <- s0[s0[,2]>0,]
		s0[,3] <- as.numeric(s0[,3])
		s0 <- s0[order(s0[,2],decreasing=T),]
		s0 <- s0[!duplicated(s0[,8]),]
		rownames(s0) <- s0[,8]
		s2 <- s0[,3,drop=F]
		
		s0 <- as.data.frame(xlsx[[3]]) 
		s0 <- s0[s0[,2]>0,]
		s0[,3] <- as.numeric(s0[,3])
		s0 <- s0[order(s0[,2],decreasing=T),]
		s0 <- s0[!duplicated(s0[,8]),]
		rownames(s0) <- s0[,8]
		s3 <- s0[,3,drop=F]
		
		olp <- intersect(rownames(s1), intersect(rownames(s2),rownames(s3)))
		
		f3 <- cbind(vWat=s1[olp,], cbind(liver=s2[olp,], gastro=s3[olp,]) )
		rownames(f3) <- olp
		
		f3 <- cbind(f3, meanLogFC=apply(f3,1,mean))
		
		f3 <- transferMouseToHuman(f3)
		
		f3 <- f3[,3,drop=F]
		rownames(f3) <- transferSymbol(rownames(f3))
		f3 <- rm_duplicates(f3)
	}else{
		dataPath <- "../../SpaCE/code/CytoSig_prediction-master/data/bulk/inflam/"
		allFiles <- list.files(dataPath)
		
		flag1 <- grepl(dataset,allFiles)
		flag2 <- grepl(".diff.1.cntmap",allFiles)
		
		fileName_full <- allFiles[flag1&flag2]
		fileName_full <- gsub(".diff.1.cntmap","",fileName_full)
			
		f3 <- as.matrix(read.csv(gzfile(paste0(dataPath,fileName_full,".diff.1")),as.is=T,sep="\t"))
		rownames(f3) <- transferSymbol(rownames(f3))
		f3 <- rm_duplicates(f3)
	}
	
	
	f3 <- f3[,!grepl("Placebo",colnames(f3),fixed=T),drop=F]
	f3 <- f3[,!grepl("non.responder",colnames(f3),fixed=T),drop=F]
	
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
		f3 <- f3[,c("TWEAK.ACHN.xenograft.4h","TWEAK.ACHN.xenograft.8h","TWEAK.ACHN.xenograft.72h","TWEAK.ACHN.xenograft.24h"),drop=F]
	}
	if(item=="TNF_GSE48498")
	{
		f3 <- f3[,"TNFA.whole.blood.cells",drop=F]
	}
		
	f3 <- apply(f3,1,function(x) mean(x,na.rm=T))
	f3 <- scale(f3)[,1]
	
	clinicalComb[paste0(item,"LigandExp"),"item"] <- item
	clinicalComb[paste0(item,"LigandExp"),"compSig"] <- "LigandExp"
	clinicalComb[paste0(item,"LigandExp"),"value"] <- f3[gene]
	
	
	if(gene=="TGFB") gene <- c("TGFB1","TGFB2","TGFB3")
	clinicalComb[paste0(item,"ReceptorExp"),"item"] <- item
	clinicalComb[paste0(item,"ReceptorExp"),"compSig"] <- "ReceptorExp"
	clinicalComb[paste0(item,"ReceptorExp"),"value"] <- mean(f3[unique(LRdb[LRdb[,1]%in%gene,2])], na.rm=TRUE)

	clinicalComb[paste0(item,"LRsumExp"),"item"] <- item
	clinicalComb[paste0(item,"LRsumExp"),"compSig"] <- "LRsumExp"
	clinicalComb[paste0(item,"LRsumExp"),"value"] <- mean(f3[c(gene,unique(LRdb[LRdb[,1]%in%gene,2]))], na.rm=TRUE)
}
		






	mat = reshape2::dcast( clinicalComb , compSig~item )
	rownames(mat) <- mat[,1]
	mat <- t(mat[,-1])
	mat <- mat[!is.na(mat[,"SecAct"]),]
	

	TNFSF12_GSE42048 <- read.csv("/data/rub2/project/Secretome/data/GSE42048.HG-U133_Plus_2.rma",sep="\t",row.names=1)
	TNFSF12_GSE42049 <- read.csv("/data/rub2/project/Secretome/data/GSE42049.HG-U133_Plus_2.rma",sep="\t",row.names=1)
	
	TNFSF12_GSE42048_logFC <- rowMeans(TNFSF12_GSE42048[,22:41])-rowMeans(TNFSF12_GSE42048[,42:62])
	TNFSF12_GSE42049_logFC <- rowMeans(TNFSF12_GSE42049[,1:4])-rowMeans(TNFSF12_GSE42049[,5:8])
	
	TNFSF12_GSE42048_logFC <- scale(TNFSF12_GSE42048_logFC)
	TNFSF12_GSE42049_logFC <- scale(TNFSF12_GSE42049_logFC)
		
	mat["TNFSF12_GSE42048","LigandExp"] <- TNFSF12_GSE42048_logFC["205611_at@TNFSF12",1] 
	mat["TNFSF12_GSE42048","ReceptorExp"]<- TNFSF12_GSE42048_logFC["218368_s_at@TNFRSF12A",1]
	mat["TNFSF12_GSE42048","LRsumExp"] <- mean(TNFSF12_GSE42048_logFC[c("205611_at@TNFSF12","218368_s_at@TNFRSF12A"),1])
	
	mat["TNFSF12_GSE42049","LigandExp"] <- TNFSF12_GSE42049_logFC["205611_at@TNFSF12",1] 
	mat["TNFSF12_GSE42049","ReceptorExp"] <- TNFSF12_GSE42049_logFC["205611_at@TNFSF12",1] 
	mat["TNFSF12_GSE42049","LRsumExp"] <- mean(TNFSF12_GSE42049_logFC[c("205611_at@TNFSF12","218368_s_at@TNFRSF12A"),1])


	mat <- mat[order(apply(mat,1,function(x) mean(x,na.rm=T))),]
	mat <- mat[,order(apply(mat,2,function(x) median(x,na.rm=T)))]

	
	write.csv(mat,paste0(validationPath,"validation_clinical_final.csv"),quote=F)
	mat <- read.csv(paste0(validationPath,"validation_clinical_final.csv"), row.names=1)
	
	
	meta <- read.csv(paste0(validationPath,"cytokine_disease_GEOID.csv"),as.is=T,row.names=1)
	for(i in 1:nrow(meta))
	{
		if(meta[i,"Organism"]=="Mouse")
		(
			meta[i,"Dataset_alt"] <- paste0(meta[i,"Dataset_alt"], "^")
		)
	}
	meta[meta[,"Dataset_alt"]=="TGFB1&2&3_Breast cancer^","Dataset_alt"] <- "TGFB1&2&3_Breast cancer^#"
	
	
	
	meta_sorted <- meta[rownames(mat),c(3,2,1,4)]
	write.csv(meta_sorted,paste0(validationPath,"cytokine_disease_GEOID_sorted.csv"))

	
	rownames(mat) <- meta[rownames(mat),"Dataset_alt"]
	rownames(mat) <- paste0("Anti-",rownames(mat))
	
	
	
	
	#png(paste0(validationPath,"validation_clinical_final.png"), width = 17.5, height = 19, res=200, units = "cm")
	pdf(paste0(validationPath,"validation_clinical_final.pdf"), width = 8.1, height = 7.8)
	#pdf(paste0(validationPath,"validation_clinical_final_2506.pdf"), width = 18, height = 7.8)

	library(ComplexHeatmap)
	disease_vec <- meta_sorted[,"Disease"]
	dataset_vec <- rownames(mat)
	
	#row_ha <- rowAnnotation(
	#	z = anno_boxplot(as.matrix(mat), height = unit(3, "cm") )
	#)
	row_ha <- rowAnnotation(
		Disease = disease_vec,
		col = list(Disease = c( "Cancer" = "#88cdbc", "Inflammatory" = "#e0d579", "Aging"="#eca680")),
		labels = anno_text(dataset_vec, which = "row", gp = gpar(fontsize = 12) ), 
    	width = max(grobWidth(textGrob(dataset_vec)))
	)
	
	column_ha <- columnAnnotation(
		"Activity Change" = anno_boxplot(as.matrix(mat), height = unit(3, "cm"), gp = gpar(fill = sigColors[colnames(mat)]) ),
		annotation_name_side = "left"
	)
			
	ht <- Heatmap(as.matrix(mat),
		name = "Activity Change",
		rect_gp = gpar(col = "white", lwd = 2),
		col = circlize::colorRamp2(c(-5, 0,5), c("#91bfdb", "white", "#fc8d59")),
		row_names_max_width = max_text_width(
	        rownames(mat), 
	        gp = gpar(fontsize = 12)
	        ),
	    column_names_max_height = max_text_width(
	        colnames(mat), 
	        gp = gpar(fontsize = 12)
	        ),
	    show_row_names = FALSE,
	    show_column_names = TRUE,
	    column_names_rot = 38,
	    top_annotation = column_ha,
	    right_annotation = row_ha,
    	cluster_rows = FALSE,
		cluster_columns = FALSE,
		cell_fun = function(j, i, x, y, width, height, fill) {grid.text(sprintf("%.1f", mat[i, j]), x, y, gp = gpar(fontsize = 10))}
	)
	draw(ht)

	dev.off()
	
	
	
	
	
	mat <- read.csv(paste0(validationPath,"validation_clinical_final.csv"),row.names=1)
	
	mat["IL1B_GSE57253",1] <- -25.328
	mat["TGFB_GSE174686",1] <- -13.515
	
	fg.df <- data.frame(gene=rownames(mat),value=mat[,1])
	fg.df[10,1] <- "TGFBpan_GSE174686"
	fg.df[,2] <- as.numeric(fg.df[,2])
	fg.df <- fg.df[order(fg.df[,2]),]
	fg.df[,1] <- factor(fg.df[,1], levels=fg.df[,1])
	
	library(ggplot2)
	p1 <- ggplot(fg.df,aes(x = gene, y = value))+
		geom_bar(stat="identity", width=0.01, color="grey66")+
		geom_hline(yintercept=0, color = "grey", linewidth=0.8)+
		geom_point(color="#619CFF",size=0.8)+
		xlab(" ")+
		ylab(paste0("Activity Change"))+
		theme_classic()+ 
		theme(
		  panel.grid = element_blank(),
		  panel.background = element_blank(),
		  axis.text = element_text(colour = "black"),
		  axis.text.x = element_text(size=7.5, angle = 90, hjust = 1, vjust = 0.5),
		  axis.title = element_text(colour = "black")
		)

	ggsave(paste0(validationPath,"validation_clinical_final_SecAct.png"), p1, width = 8.2, height =7, dpi=400, units = "cm",)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	



items <- c("VEGFA_GSE72951","VEGFA_E-MTAB-3267")

compSigs <- c(
	0,11,112,
	36,618,
	66,77,88,
	61811,61822,61833,61844,61855,61866,61877,61888,61899,61800,
	100011,100012,100021,100022,100031,100032,100041,100042,100051,100052
	)
compSigNames <- c(
	"CytoSig","NicheNet.v1","NicheNet.v2",
	"ImmuneDic","SecAct",
	"LigandExp","ReceptorExp","LRsumExp",
	"SecAct.w1","SecAct.w2","SecAct.v1","SecAct.v2","SecAct.v2_0.9","SecAct.v2_0.8","SecAct.w2_0.9","SecAct.1k_0.9","SecAct.1k_residual","SecAct.1k",
	"SecAct.1k_11","SecAct.1k_12","SecAct.1k_21","SecAct.1k_22","SecAct.1k_31","SecAct.1k_32","SecAct.1k_41","SecAct.1k_42","SecAct.1k_51","SecAct.1k_52"
	)
CompSigAlphas <- c(
	1e+07,1e+07,1e+07,
	50000,5e+05,
	10000,10000,10000,
	5e+05,5e+05,5e+05,5e+05,5e+05,5e+05,5e+05,5e+05,5e+05,5e+05,
	5e+05,5e+05,5e+05,5e+05,5e+05,5e+05,5e+05,5e+05,5e+05,5e+05
)


compSigs <- c(
	0,11,112,
	36,100042,
	66,77,88
	)
compSigNames <- c(
	"CytoSig","NicheNet.v1","NicheNet.v2",
	"ImmuneDic","SecAct",
	"LigandExp","ReceptorExp","LRsumExp"
	)
CompSigAlphas <- c(
	1e+07,1e+07,1e+07,
	50000,5e+05,
	10000,10000,10000
)



for(compSig in compSigs)
{	
	clinical <- data.frame()
	
	for(item in items)
	{
		gene <- strsplit(item,"_") [[1]][1]
		dataset <- strsplit(item,"_") [[1]][2]

		#for(a in c(10000,50000,1e+05,5e+05,1e+06,5e+06,1e+07,5e+07) )
		#{
			a <- CompSigAlphas[compSigs==compSig]
			dataPath <- "../../SpaCE/code/CytoSig_prediction-master/data/bulk/tumor/"
			
			if(compSig%in%c(66,77,88))
			{
				if(item=="VEGFA_GSE72951")
				{
					f3 <- as.matrix(read.csv(gzfile(paste0(dataPath,"GSE72951.self_subtract.gz")),as.is=T,sep="\t"))
					f33 <- read.csv(paste0(dataPath,"GSE72951.OS.Bevacizumab"),as.is=T,sep="\t")
				}else{
					f3 <- as.matrix(read.csv(gzfile(paste0(dataPath,"E-MTAB-3267.norm_subtract.gz")),as.is=T,sep="\t",check.names=F))
					f33 <- read.csv(paste0(dataPath,"E-MTAB-3267.PFS"),as.is=T,sep="\t")
				}
				
				if(compSig==66)
				{
					gene <- c("VEGFA")
				}else if(compSig==77){
					gene <- c("FLT1", "KDR", "NRP1")
				}else{
					gene <- c("VEGFA","FLT1", "KDR", "NRP1")
				}
				
			}else{
				if(!file.exists(paste0(validationPath,"clinical/",item,"_",compSig,"_",a,".RData"))) next
				load(paste0(validationPath,"clinical/",item,"_",compSig,"_",a,".RData"))

				f3 <- as.matrix(res$zscore)
				f3 <- expand_rows(f3)
			
				if(item=="VEGFA_GSE72951")
				{
					f33 <- read.csv(paste0(dataPath,"GSE72951.OS.Bevacizumab"),as.is=T,sep="\t")
				}else{
					f33 <- read.csv(paste0(dataPath,"E-MTAB-3267.PFS"),as.is=T,sep="\t")
				}
			
				gene <- c("VEGFA")
			}
			
			
			library(survival)
			
			olp <- intersect(rownames(f33),colnames(f3))
			neww <- cbind(f33[olp,], value=colMeans(f3[rownames(f3)%in%gene,olp,drop=F]))
			
			if(item%in%c("VEGFA_GSE72951"))
			{
				coxmodel_fit <- coxph(Surv(OS, Event) ~ value, data = neww)
			}else{
				coxmodel_fit <- coxph(Surv(PFS, Event) ~ value, data = neww)
			}
			
			coxmodel_obj <- summary(coxmodel_fit)
			
			clinical[paste0(item,"_",compSig,"_",a),"item"] <- item
			clinical[paste0(item,"_",compSig,"_",a),"compSig"] <- paste0(compSig,"_",a)
			clinical[paste0(item,"_",compSig,"_",a),"value"] <- coxmodel_obj$coefficients[1,"z"]
			
		#} # a 
		
	}
	
	
	mat = reshape2::dcast( clinical , compSig~item )
	rownames(mat) <- mat[,1]
	mat <- t(mat[,-1])
	
	mat <- mat[order(apply(mat,1,function(x) median(x,na.rm=T))),,drop=F]
	mat <- mat[,order(apply(mat,2,function(x) median(x,na.rm=T))),drop=F]
	
	
	library(ComplexHeatmap)
	
	#png(paste0(validationPath,"validation_clinical_VEGF_",compSigNames[which(compSigs==compSig)],".png"), width = 20, height = 10, res=200, units = "cm")
	#
	#row_ha <- rowAnnotation(
	#	z = anno_boxplot(as.matrix(mat), height = unit(3, "cm") ) 
	#)
	#	column_ha <- columnAnnotation(
	#	z = anno_boxplot(as.matrix(mat), height = unit(3, "cm") ) 
	#)
	#		
	#ht <- Heatmap(as.matrix(mat),
	#	column_title=compSigNames[which(compSigs==compSig)],
	#	name = "Risk z",
	#	col = circlize::colorRamp2(c(-3, 0,3), c("green", "white", "red")),
	#	row_names_max_width = max_text_width(
	#        rownames(mat), 
	#        gp = gpar(fontsize = 12)
	#        ),
	#    column_names_max_height = max_text_width(
	#        colnames(mat), 
	#        gp = gpar(fontsize = 12)
	#        ),
	#    top_annotation = column_ha,
	#    right_annotation = row_ha,
   #	cluster_rows = FALSE,
	#	cluster_columns = FALSE,
	#	cell_fun = function(j, i, x, y, width, height, fill) {grid.text(sprintf("%.3f", mat[i, j]), x, y, gp = gpar(fontsize = 10))}
	#)
	#
	#draw(ht)
	#dev.off()
	
	#clinicalBest <- clinical[clinical[,2]==colnames(mat)[1],,drop=F]
	clinicalBest <- clinical
	
	if(compSig==0)
	{
		clinicalComb <- clinicalBest
	}else{
		clinicalComb <- rbind(clinicalComb,clinicalBest)
	}
	
}


for(compSig in compSigs)
{	
	clinicalComb[clinicalComb[,2]==paste0(compSig,"_",CompSigAlphas[compSigs==compSig]),2] <- compSigNames[compSigs==compSig]	
}


mat = reshape2::dcast( clinicalComb , compSig~item )
rownames(mat) <- mat[,1]
mat <- t(mat[,-1])
mat <- mat[!is.na(mat[,"SecAct"]),]

mat <- mat[order(apply(mat,1,function(x) median(x,na.rm=T))),]
mat <- mat[,order(apply(mat,2,function(x) median(x,na.rm=T)))]

write.csv(mat,paste0(validationPath,"validation_clinical_VEGF_summary.csv"),quote=F)


png(paste0(validationPath,"validation_clinical_VEGF_summary_heatmap.png"), width = 30, height = 10, res=200, units = "cm")

library(ComplexHeatmap)
row_ha <- rowAnnotation(
	z = anno_boxplot(as.matrix(mat), height = unit(3, "cm") ) 
)
	column_ha <- columnAnnotation(
	z = anno_boxplot(as.matrix(mat), height = unit(3, "cm") ) 
)		
Heatmap(as.matrix(mat),
	name = "Risk z",
	col = circlize::colorRamp2(c(-3, 0,3), c("green", "white", "red")),
	row_names_max_width = max_text_width(
        rownames(mat), 
        gp = gpar(fontsize = 12)
        ),
    column_names_max_height = max_text_width(
        colnames(mat), 
        gp = gpar(fontsize = 12)
        ),
    top_annotation = column_ha,
    right_annotation = row_ha,
	cluster_rows = FALSE,
	cluster_columns = FALSE,
	cell_fun = function(j, i, x, y, width, height, fill) {grid.text(sprintf("%.3f", mat[i, j]), x, y, gp = gpar(fontsize = 10))}
)

dev.off()
	


clinicalComb <- clinicalComb[!clinicalComb[,2]%in%c("ImmuneDicAll","SecAct398"),]
clinicalComb[clinicalComb[,2]=="ImmuneDicFilter",2] <- "ImmuneDic"

sigOrder <- data.frame()
for(i in unique(clinicalComb[,"compSig"]))
{
	sigOrder[i,"r_mean"] <- mean(clinicalComb[clinicalComb[,"compSig"]==i,"value"])
}
sigOrder <- sigOrder[order(sigOrder[,"r_mean"]),,drop=F]

clinicalComb[,2] <- factor(clinicalComb[,2], levels=rownames(sigOrder))


library(ggplot2)
p1 <- ggplot(clinicalComb, aes(x=compSig, y=value, fill=compSig)) +
  geom_bar(stat="identity", color="white", width=1, alpha=0.7) +
  scale_fill_manual( values=sigColors[rownames(sigOrder)] )+
  ylab("Risk score z")+
	xlab(" ")+
	theme_classic()+ 
	theme(
		panel.background = element_blank(),
		panel.grid = element_blank(),
		axis.title = element_text(colour = "black"),
		axis.text = element_text(colour = "black"),
		axis.text.x = element_text(angle = 38,hjust = 1),
		strip.background = element_rect(colour="white", fill="white"),
		legend.position = "none"
	)+ facet_wrap(~ item, ncol =2) #scales = "free"
ggsave(paste0(validationPath,"validation_clinical_VEGF_summary_bar.png"), p1, width = 12.5, height = 8, dpi=500, units = "cm", limitsize =FALSE)


clinicalCombMean <- data.frame(
	compSig=colnames(mat),
	value=colMeans(mat))

library(ggplot2)
p1 <- ggplot() +
	geom_bar(aes(x=compSig, y=value, fill=compSig), data=clinicalCombMean, stat="identity", color="white", width=1, alpha=0.88) +
	geom_jitter(aes(x=compSig, y=value), data=clinicalComb, width=0.01, color="grey60") +
	scale_fill_manual( values=sigColors[rownames(sigOrder)] )+
	ylab("Risk score z")+
	xlab(" ")+
	theme_classic()+ 
	theme(
		panel.background = element_blank(),
		panel.grid = element_blank(),
		axis.title = element_text(colour = "black"),
		axis.text = element_text(colour = "black"),
		axis.text.x = element_text(angle = 38,hjust = 1),
		strip.background = element_rect(colour="white", fill="white"),
		legend.position = "none"
	)
ggsave(paste0(validationPath,"validation_clinical_VEGF_summary_bar2.png"), p1, width = 7.2, height = 7.8, dpi=500, units = "cm", limitsize =FALSE)










items <- c("VEGFA_GSE72951","VEGFA_E-MTAB-3267")
dataPath <- "../../SpaCE/code/CytoSig_prediction-master/data/bulk/tumor/"
source(file.path("./survival_util.R"))
sig = 100042
a = 5e+05

for(item in items)
{
	load(paste0(validationPath,"clinical/",item,"_",sig,"_",a,".RData"))
	Act <- as.matrix(res$zscore)
	Act <- expand_rows(Act)
	data <- t(as.matrix(Act)) [,c("LY86","VEGFA")]
	
	
	if(item=="VEGFA_GSE72951")
	{
		survival <- read.csv(paste0(dataPath,"GSE72951.OS.Bevacizumab"),as.is=T,sep="\t")
		xtext <- "Overall (Months)"
		widthValue <- 9.7
		addtext <- 20
	}else{
		survival <- read.csv(paste0(dataPath,"E-MTAB-3267.PFS"),as.is=T,sep="\t")
		xtext <- "Progression-Free (Months)"
		widthValue <- 9.5
		addtext <- 50
	}
	
	out <- Beibei_revised(data, survival, margin=5)
	data <- out[[1]]
	survival <- out[[2]]
	result <- out[[3]]
	cutoff <- round(result["VEGFA","thres.opt"],3)

	
		
	library(survival)
	olp <- intersect(rownames(survival),rownames(data))
	X_olp <- survival[olp,,drop=F]
	Y_olp <- data[olp,,drop=F]
	
	comb <- cbind(X_olp, Act=Y_olp[,"VEGFA"])
	
	if(item=="VEGFA_GSE72951")
	{
		coxmodel_fit <- coxph(Surv(OS, Event) ~ ., data = comb)
	}else{
		coxmodel_fit <- coxph(Surv(PFS, Event) ~ ., data = comb)
	}
	coxmodel_obj <- summary(coxmodel_fit)
	zs <- coxmodel_obj$coefficients["Act","z"]
	pv <- coxmodel_obj$coefficients["Act","Pr(>|z|)"]
	pv <- pv/2 # two sided -> one sided
	
	
	library(survival)
	surv_o <- Surv(survival[,1],survival[,2])
	
	groups <- as.character(data[,"VEGFA"]>cutoff)
	surv_d <- survdiff(surv_o ~ groups)
	#pv <- 1 - pchisq(surv_d$chisq, length(surv_d$n) - 1)
				
	coxmodel_fit <- coxph(surv_o ~ groups)
	coxmodel_obj <- summary(coxmodel_fit)
	#zs <- coxmodel_obj$coefficients["groupsTRUE","z"]

	surv_c  <- survfit(surv_o ~ groups)
	
	
	jpeg(file=paste0(validationPath,item,"_",sig,".jpg"),width=widthValue, height=11, units="cm",res=300)
	plot(surv_c,mark.time=T,col=colors_surv,
		#main=item, # a gap exists between plot and text
		#xlab=xtext, # a gap exists between plot and text
		ylab="Percentage",
		frame=F,lwd=2,las=1
		)
	title(item, adj = 0.5, line = 0.4, font.main= 1)
	mtext(side=1, text=xtext, line=2.2)
	legend("topright",1,
				legend=c(
					paste0("High (n=",sum(groups=="TRUE"),")"),
					paste0("Low (n=",sum(groups=="FALSE"),")") ),
				bty="n",cex=1.1,lwd=3,col=rev(colors_surv))
	legend("right",1,
				legend=paste0("z = ", round(zs,2), "\np = ", signif(pv,2) ),
				bty="n",cex=1.1)
	dev.off()
	
	
	
	library("survminer")
	surv.df <- cbind(survival,groups)
	
	if(item=="VEGFA_GSE72951")
	{
		fit <- survfit(Surv(OS, Event) ~ groups, data = surv.df)
	}else{
		fit <- survfit(Surv(PFS, Event) ~ groups, data = surv.df)
	}
	
	p2 <- ggsurvplot(fit, data=surv.df, palette = colors_surv, legend.labs = c(paste0("Low (n=",sum(groups=="FALSE"),")"),paste0("High (n=",sum(groups=="TRUE"),")")) )$plot+
		xlab(xtext)+
		ylab("Percentage")+
		annotate("text", x = addtext, y=0.88, label = paste0("z = ", round(zs,2), "\np = ", signif(pv,2)), size=4 )

	ggsave(paste0(validationPath,item,"_",sig,".pdf"), p2, width = 3.2, height = 3.5)
	
}








item <- "VEGFA_E-MTAB-3267"

load(paste0(validationPath,"clinical/",item,"_37_5e+05.RData"))

f3 <- as.matrix(res$zscore)
f33 <- read.csv("E-MTAB-3267.sdrf.txt",as.is=T,sep="\t")
f33 <- f33[f33[,"Characteristics.disease."]=="Tumor",]
f33 <- f33[,c("Characteristics.individual.","Characteristics.sunitinib.response.")]
f33[f33[,2]=="CLINICAL BENEFIT",2] <- "CR"
f33[f33[,2]%in%c("CR","PR"),2] <- "CR/PR"
f33[f33[,2]%in%c("SD","PD"),2] <- "SD/PD"
f33 <- f33[order(f33[,1]),]
f333 <- cbind(f33,f3["VEGFA",])
colnames(f333) <- c("sample","Response","Activity")

library(ggplot2)
p2 <- ggplot(f333,aes(x = Response, y = Activity, fill=Response))+
	geom_boxplot(outlier.shape=NA, alpha=0.3, width=0.5)+
	geom_jitter(size=0.8, alpha=0.3)+
	xlab(" ")+
	ylab("Activity")+
	theme_bw()+ 
	theme(
	  panel.grid = element_blank(),
	  panel.background = element_blank(),
	  plot.title = element_text(size=11, hjust = 0.5),
	  axis.text = element_text(size=11,colour = "black"),
	  axis.title = element_text(size=11,colour = "black"),
	  legend.position="none"
	) 
ggsave(paste0(validationPath,"validation_clinical_final3.png"), p2, width = 16, height = 10, dpi=500, units = "cm", limitsize =FALSE)




}


