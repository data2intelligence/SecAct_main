source("Secretome_s0_path.R")
source("../../SpaCE/code/ifun.R")
source("importHPA.R")

group <- "Membrane_"
group <- ""

CTL.corrected <- "CTL.corrected_"
CTL.corrected <- ""


if(TRUE)
{

#SWARM -t 10 -g 20 --time 00:30:00
args = commandArgs(trailingOnly=TRUE)
cancer <- args[1]
datatype <- args[2]


dataPath <- "/data/rub2/data/Immunotherapy/"
expr <- as.matrix(read.csv(gzfile(paste0(dataPath,cancer,".",datatype,".gz")),as.is=T,sep="\t",check.names=F))
expr <- expr[apply(expr,1,function(x) sum(x>0)>5),]

rownames(expr) <- transferSymbol(rownames(expr))
expr <- rm_duplicates(expr)

write.table(expr, paste0(predictionPath,"immunotherapy/",cancer,".txt"),sep="\t",quote=F)	
stop

cdata_T_minusBG <- expr-rowMeans(expr)

write.table(cdata_T_minusBG, paste0(predictionPath,"immunotherapy/",cancer,".diff"),sep="\t",quote=F)	


#runCytoSig <- 
#	paste0(
#	"CytoSig_run.py -a 1000000 -s 15 -i ",
#	paste0(predictionPath,"immunotherapy/",cancer,".diff"),
#	" -o ",
#	paste0(predictionPath,"immunotherapy/",cancer)
#	)
#
#system("module load python") # run it in node first
#system(runCytoSig)


#cdata_T_minusBG_sub <- cdata_T_minusBG[rownames(cdata_T_minusBG)%in%c("CD8A","CD8B","GZMA","GZMB","PRF1"),,drop=F]
#cdata_T_minusBG_sub <- rbind(cdata_T_minusBG_sub, CTL=colMeans(cdata_T_minusBG_sub))
#write.table(cdata_T_minusBG_sub, paste0(predictionPath,"immunotherapy/",cancer,".CTL"),sep="\t",quote=F)	


library(SecAct)
if(group=="Membrane_")
{
	ref <- paste0("/data/rub2/project/Secretome/results/signatureComb/AllSigFilteredBy_Membrane_MoranI_TCGA_ICGC_0.25_ds3.tsv")
	res <- SecAct.inference(Y=cdata_T_minusBG, SigMat=ref, lambda=5e+6, nrand=1000)
	save(res, file = paste0(predictionPath,"immunotherapy/Membrane_",cancer,".RData"))
}else{
	ref <- paste0("/data/rub2/project/Secretome/results/signatureComb/AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3_1k_vst_condition_logUMI_cellType_0.9.tsv")
	res <- SecAct.inference(Y=cdata_T_minusBG, SigMat=ref, lambda=5e+5, nrand=1000)
	save(res, file = paste0(predictionPath,"immunotherapy/",cancer,".RData"))
}


}















n_downsampling <- c(4,5,6,8,10,12,14,16,18,20,40,50,60,70,80,90,100,200,400,600,800,1000,2000,4000,6000,8000,10000,12000,14000,16000,18000)
n_downsampling <- c(5,10,15,20,40,60,80,100,200,400,600,800,1000,2000,4000,6000,8000)

notGenomeWide <- c(
"HNSCC_ICB_Foy2022",
"NSCLC_PD1_Foy2022",
"PanCancer_PD1_Prat2017"
)

if(FALSE)
{

#SWARM -t 60 -g 20 --time 40:20:00
args = commandArgs(trailingOnly=TRUE)
cancer <- args[1]
datatype <- args[2]


dataPath <- "/data/rub2/data/Immunotherapy/"
expr <- as.matrix(read.csv(gzfile(paste0(dataPath,cancer,".",datatype,".gz")),as.is=T,sep="\t",check.names=F))
expr <- expr[apply(expr,1,function(x) sum(x>0)>5),]

rownames(expr) <- transferSymbol(rownames(expr))
expr <- rm_duplicates(expr)

cdata_T_minusBG <- expr-rowMeans(expr)

#write.table(cdata_T_minusBG, paste0(predictionPath,"immunotherapy/",cancer,".diff"),sep="\t",quote=F)	


#runCytoSig <- 
#	paste0(
#	"CytoSig_run.py -a 1000000 -s 15 -i ",
#	paste0(predictionPath,"immunotherapy/",cancer,".diff"),
#	" -o ",
#	paste0(predictionPath,"immunotherapy/",cancer)
#	)
#
#system("module load python") # run it in node first
#system(runCytoSig)


#cdata_T_minusBG_sub <- cdata_T_minusBG[rownames(cdata_T_minusBG)%in%c("CD8A","CD8B","GZMA","GZMB","PRF1"),,drop=F]
#cdata_T_minusBG_sub <- rbind(cdata_T_minusBG_sub, CTL=colMeans(cdata_T_minusBG_sub))
#write.table(cdata_T_minusBG_sub, paste0(predictionPath,"immunotherapy/",cancer,".CTL"),sep="\t",quote=F)	


library(SecAct)
#if(group=="Membrane_")
#{
#	ref <- paste0("/data/rub2/project/Secretome/results/signatureComb/AllSigFilteredBy_Membrane_MoranI_TCGA_ICGC_0.25_ds3.tsv")
#	res <- SecAct.inference(Y=cdata_T_minusBG, SigMat=ref, lambda=5e+6, nrand=1000)
#	save(res, file = paste0(predictionPath,"immunotherapy/Membrane_",cancer,".RData"))
#}else{
#	ref <- paste0("/data/rub2/project/Secretome/results/signatureComb/AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3.tsv")
#	res <- SecAct.inference(Y=cdata_T_minusBG, SigMat=ref, lambda=5e+5, nrand=1000)
#	save(res, file = paste0(predictionPath,"immunotherapy/",cancer,".RData"))
#}


ref <- paste0("/data/rub2/project/Secretome/results/signatureComb/AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3.tsv")
f1 <- read.table(ref)

cdata_T_minusBG <- cdata_T_minusBG[rownames(cdata_T_minusBG)%in%rownames(f1),]

for(n in n_downsampling)
{
	if(cancer%in%notGenomeWide) next
	
	dir.create(paste0(predictionPath,"immunotherapy/",n))
	
	set.seed(123)
	for(i in 1:10)
	{
		if(n<nrow(cdata_T_minusBG))
		{
			geneSub <- sample(1:nrow(cdata_T_minusBG), n)
		}else{
			geneSub <- 1:nrow(cdata_T_minusBG)
		}
		
		cdata_T_minusBG_sub <- cdata_T_minusBG[geneSub,]
		res <- SecAct.inference(Y=cdata_T_minusBG_sub, SigMat=ref, lambda=5e+5, nrand=1000)
	
		save(res, file = paste0(predictionPath,"immunotherapy/",n,"/",i,"_",cancer,".RData"))
	}
}


}


#stop("xx")





			
if(FALSE)
{


notICB <- c(
"CCRCC_Everolimus_Braun2020",
"Hepatocellular_Sorafenib_Finn2020",
"mRCC_Sunitinib_McDermott2018",
"Breast_aHER2+Chemo_Sammut2022",
"Breast_Chemo_Sammut2022",
"Mouse_PD1_Chen2020",
"Mouse_PDL1_Meskini2021",
"SCLC_Placebo_Nabet2024",
"RCC_Sunitinib_Motzer2020"
)

notGenomeWide <- c(
"HNSCC_ICB_Foy2022",
"NSCLC_PD1_Foy2022",
"PanCancer_PD1_Prat2017"
)
notGenomeWide_alt <- c(
"HNSCC_ICB_Foy2022.OS",
"NSCLC_PD1_Foy2022.OS",
"PanCancer_PD1_Prat2017.PFS"
)

notGenomeWide_alt2 <- c(
"HeadNeck_ICB_Foy 2022",
"Lung (Non-Small Cell)_PD1_Foy 2022",
"PanCancer_PD1_Prat 2017"
)

notSolidTumor <- c(
"Melanoma_ICB_Lozano2022" # from blood
)

notExpr <- c(
"HeadNeck_Pembrolizumab_Cristescu2018", # with CNV and mutation
"Melanoma_Pembrolizumab_Cristescu2018", # with CNV and mutation
"PanCancer_Pembrolizumab_Cristescu2018", # with CNV and mutation
"Melanoma_CTLA4_Snyder2014", # with CNV and mutation
"NSCLC_Pembrolizumab_Rizvi2015" # with CNV and mutation
)

notResonse <- c(
"PanCancer_ICB_Sung2023" # only irAE
)


dataPath <- "/data/rub2/data/Immunotherapy/"
allFiles <- list.files(dataPath)
allFiles_OS <- allFiles[grepl("OS",allFiles)]
allFiles_PFS <- allFiles[grepl("PFS",allFiles)]
allFiles_RECIST <- allFiles[grepl("RECIST",allFiles)]
allFiles_response <- allFiles[grepl("response",allFiles)]
allFiles <- sort(c(allFiles_OS,allFiles_PFS,allFiles_RECIST,allFiles_response))
allFiles <- allFiles[!grepl("NSCLC_ICB_Ravi2023",allFiles)]

x1 <- sapply(strsplit(allFiles,".",fixed=T),function(x) return(x[1]))
x2 <- sapply(strsplit(allFiles,".",fixed=T),function(x) return(x[2]))
x1 <- c(x1,"NSCLC_ICB_Ravi2023.Adeno","NSCLC_ICB_Ravi2023.Squamous")
x2 <- c(x2,"OS","OS")


mat <- data.frame(study=x1, meta=x2)
rownames(mat) <- c(allFiles,"NSCLC_ICB_Ravi2023.Adeno.OS","NSCLC_ICB_Ravi2023.Squamous.OS")

mat <- mat[!mat[,1]%in%notICB,]
mat <- mat[!mat[,1]%in%notSolidTumor,]
mat <- mat[!mat[,1]%in%notExpr,]
#mat <- mat[!mat[,1]%in%notGenomeWide,]

datasets <- c()
for(study in unique(mat[,1]))
{
	mat_sub <- mat[mat[,1]==study,]
	
	if("OS"%in%mat_sub[,2] | "OS_Naive"%in%mat_sub[,2] | "OS_Prog"%in%mat_sub[,2] |
		"OS_Nivo+Chemo"%in%mat_sub[,2] | "OS_Nivo+Sotiga+Chemo"%in%mat_sub[,2] | "OS_Sotiga+Chemo"%in%mat_sub[,2]
	)
	{
		datasets <- c(datasets, rownames(mat_sub)[mat_sub[,2]%in%c("OS","OS_Naive","OS_Prog","OS_Nivo+Chemo","OS_Nivo+Sotiga+Chemo","OS_Sotiga+Chemo")])
		next
	}else if("PFS"%in%mat_sub[,2]){
		datasets <- c(datasets, rownames(mat_sub)[mat_sub[,2]%in%"PFS"])
		next
	}else if("RECIST"%in%mat_sub[,2]){
		datasets <- c(datasets, rownames(mat_sub)[mat_sub[,2]%in%"RECIST"])
		next
	}else{
		datasets <- c(datasets, rownames(mat_sub)[mat_sub[,2]%in%c("response","response_Basal","response_Luminal")])
	}
}

# already have subtype
datasets <- datasets[!datasets%in%"Breast_Pembrolizumab_Wolf2022.response"]

#low exon coverage
datasets <- datasets[!datasets%in%"Colorectal_ICB_Thibaudin2023.RECIST"]
datasets <- datasets[!datasets%in%"Esophageal_Atezolizumab_VanDenEnde2021.response"]

#low sample number 
datasets <- datasets[!datasets%in%"Melanoma_CTLA4_Campbell2023.RECIST"]
datasets <- datasets[!datasets%in%"Melanoma_PD1-to-CTLA4_Campbell2023.RECIST"]
datasets <- datasets[!datasets%in%"Melanoma_CTLA4_Roh2017.response"] # targeted
datasets <- datasets[!datasets%in%"Melanoma_PD1_Roh2017.response"] # targeted

# new
datasets <- datasets[!datasets%in%"Urothelial_Atezolizumab_Hamidi2024-IMvigor010.OS"]
datasets <- datasets[!datasets%in%"Urothelial_Atezolizumab_Hamidi2024-IMvigor210.OS"]
datasets <- datasets[!datasets%in%"Urothelial_Chemotherapy_Hamidi2024-IMvigor130.OS"]
datasets <- datasets[!datasets%in%"Urothelial_Chemotherapy_Hamidi2024-IMvigor211.OS"]
datasets <- datasets[!datasets%in%"Urothelial_Surveillance_Hamidi2024-IMvigor010.OS"]
datasets <- datasets[!datasets%in%"Urothelial_Observation_Powles2024.OS"]

datasets <- datasets[!datasets%in%"NSCLC_Docetaxel_Patil2022-OAK.OS"]
datasets <- datasets[!datasets%in%"NSCLC_Docetaxel_Patil2022-POPLAR.OS"]
datasets <- datasets[!datasets%in%"RCC_Sunitinib_Motzer2020-IMmotion151.PFS"]

datasets <- datasets[!datasets%in%"NSCLC_PD1_Yan2025.RECIST"]

datasets <- datasets[!datasets%in%"Colorectal_Regorafenib_Eng2019.OS"]


#this treatment failed
datasets <- datasets[!datasets%in%"Pancreatic_Nivolumab_Padron2022.OS_Sotiga+Chemo"]



datasets <- sort(datasets)



xx <- sapply(strsplit(datasets,".",fixed=T),function(x) return(x[1]))
xx <- sapply(strsplit(xx,"_",fixed=T),function(x) return(x[3]))
xx <- sort(unique(xx))


CIDE <- c("Braun2020","Budhu2023","Campbell2023","Cui2021","Finn2020","Gide2019","Hu2023","Hugo2016","Jung2019","Kim2018","Lauss2017","Lee2021","Li2023","Li2024","Liu2019","Mariathasan2018","McDermott2018","Miao2018","Montoya2013","Motzer2020","Nabet2024","Padron2022","Pusztai2021","Ravi2023","Riaz2017","Roper2021","Uppaluri2020","VanAllen2015","Vos2023","Wolf2022","Yang2021","Zhao2019","Hamidi2024-IMvigor130","Hamidi2024-IMvigor211","Patil2022-OAK","Patil2022-POPLAR","Macarulla2024","Powles2024","Motzer2020-IMmotion151","Mamdani2021")
setdiff(xx,CIDE)
setdiff(CIDE,xx)


###########################
### sample count
###########################

ctype_vec <- c()
patient_count <- c()
coverage_vec <- c()
for(dataset in datasets)
{
	if(!dataset%in%c("NSCLC_ICB_Ravi2023.Adeno.OS","NSCLC_ICB_Ravi2023.Squamous.OS" ))
	{
		study <- strsplit(dataset,".",fixed=T) [[1]][1]
		clini <- strsplit(dataset,".",fixed=T) [[1]][2]
		ctype <- strsplit(clini,"_",fixed=T) [[1]][1]
		surv <- read.csv(gzfile(paste0(dataPath,dataset)),as.is=T,sep="\t",check.names=F)
	}else{
		study <- strsplit(dataset,".O",fixed=T) [[1]][1]
		ctype <- strsplit(dataset,".",fixed=T) [[1]][3]
		surv <- read.csv(gzfile(paste0(dataPath,"NSCLC_ICB_Ravi2023.OS")),as.is=T,sep="\t",check.names=F)
	}
	
	
	if(dataset%in%c(
		"Breast_Durvalumab+Olaparib_Pusztai2021.response",
		"Breast_Pembrolizumab_Wolf2022.response_Basal",
		"Breast_Pembrolizumab_Wolf2022.response_Luminal",
		"HeadNeck_Pembrolizumab_Uppaluri2020.response"))
	{
		colnames(surv)[1] <- "Response"
	}
	if(dataset=="NSCLC_Pembrolizumab_Lee2021.response")
	{
		rownames(surv) <- surv[,1]
		surv <- surv[,-1,drop=F]
	}
	if(ctype%in%c("OS"))
	{
		colnames(surv)[1] <- "OS"
		colnames(surv)[2] <- "OS.Event"	
	}
	if(ctype%in%c("PFS"))
	{
		colnames(surv)[1] <- "PFS"
		colnames(surv)[2] <- "PFS.Event"	
	}
	if(ctype%in%c("RECIST"))
	{
		colnames(surv)[1] <- "RECIST"
	}
	
	ctype_vec <- c(ctype_vec, ctype)
	
	#Exp <- read.table(paste0(predictionPath,"immunotherapy/",study,".diff"),sep="\t",check.names=F)
	#print(nrow(Exp))
	#patient_count <- c(patient_count, dim(res$zscore)[2] )

	load(paste0(predictionPath,"immunotherapy/",study,".RData"))
	
	Act <- as.matrix(res$zscore)	
	Act <- expand_rows(Act)

	matx <- t(Act)
	olp <- intersect(rownames(surv),rownames(matx))
	
	patient_count <- c(patient_count, length(olp))
	
	coverage_vec <- c(coverage_vec, ifelse(dataset%in%notGenomeWide_alt, "Targeted", "Genome wide") )
}	

rename <- read.csv(paste0(inputPath,"cohort_rename.txt"),row.names=1,header=T,sep="\t")

stat <- data.frame(
	ID=datasets, 
	Name=rename[match(datasets,rename[,1]),2],
	clinical=ctype_vec, 
	expr_count=patient_count, 
	Coverage=coverage_vec)

stat <- stat[order(stat[,"Coverage"]),]

sum(patient_count)
# 2835

write.table(stat,paste0(predictionPath,"rename.txt"),quote=F,sep="\t")





downsampling <- FALSE


#SWARM -t 2 -g 20 --time 1:00:00
args = commandArgs(trailingOnly=TRUE)
n <- args[1]


for(i in 1:10)
{

smy <- data.frame()
ctype_vec <- c()
for(dataset in datasets)
{
	if(!dataset%in%c("NSCLC_ICB_Ravi2023.Adeno.OS","NSCLC_ICB_Ravi2023.Squamous.OS" ))
	{
		study <- strsplit(dataset,".",fixed=T) [[1]][1]
		clini <- strsplit(dataset,".",fixed=T) [[1]][2]
		ctype <- strsplit(clini,"_",fixed=T) [[1]][1]
		surv <- read.csv(gzfile(paste0(dataPath,dataset)),as.is=T,sep="\t",check.names=F)
	}else{
		study <- strsplit(dataset,".O",fixed=T) [[1]][1]
		ctype <- strsplit(dataset,".",fixed=T) [[1]][3]
		surv <- read.csv(gzfile(paste0(dataPath,"NSCLC_ICB_Ravi2023.OS")),as.is=T,sep="\t",check.names=F)
	}
	
	
	if(dataset%in%c(
		"Breast_Durvalumab+Olaparib_Pusztai2021.response",
		"Breast_Pembrolizumab_Wolf2022.response_Basal",
		"Breast_Pembrolizumab_Wolf2022.response_Luminal",
		"HeadNeck_Pembrolizumab_Uppaluri2020.response"))
	{
		colnames(surv)[1] <- "Response"
	}
	if(dataset=="NSCLC_Pembrolizumab_Lee2021.response")
	{
		rownames(surv) <- surv[,1]
		surv <- surv[,-1,drop=F]
	}
	if(ctype%in%c("OS"))
	{
		colnames(surv)[1] <- "OS"
		colnames(surv)[2] <- "OS.Event"	
	}
	if(ctype%in%c("PFS"))
	{
		colnames(surv)[1] <- "PFS"
		colnames(surv)[2] <- "PFS.Event"	
	}
	if(ctype%in%c("RECIST"))
	{
		colnames(surv)[1] <- "RECIST"
	}
	
	ctype_vec <- c(ctype_vec, ctype)
	
	#print(dataset)
	#print(colnames(surv))
	#}
	#print(paste0(colnames(surv)[1]," ",dataset)) 
	#
	#if(colnames(surv)[1]%in%c("Response"))
	#{
	#	print(surv)
	#}
	#}
	
	#Act <- t(as.matrix(read.table(paste0(predictionPath,"immunotherapy/",study,".Zscore"),sep="\t",check.names=F)	))
	
	
	if(downsampling==TRUE)
	{
		if(!file.exists(paste0(predictionPath,"immunotherapy/",group,n,"/",i,"_",study,".RData"))) next
		load(paste0(predictionPath,"immunotherapy/",group,n,"/",i,"_",study,".RData"))
	}else{
		load(paste0(predictionPath,"immunotherapy/",group,study,".RData"))
	}
	
	Act <- as.matrix(res$zscore)
	Act <- expand_rows(Act)
	Act <- t(Act)
	
	
	#print(cor(Act[,"SPP1"],Act[,"CXCL9"]))
	
	matx <- Act
	
	olp <- intersect(rownames(surv),rownames(matx))
	X_olp <- surv[olp,,drop=F]
	Y_olp <- matx[olp,,drop=F]
	
	if(CTL.corrected=="CTL.corrected_")
	{
		CTL <- read.table(paste0(predictionPath,"immunotherapy/",study,".CTL"),sep="\t",check.names=F)	
		CTL <- t(as.matrix(CTL))
		X_olp <- cbind(X_olp, CTL=CTL[olp,"CTL"])
	}
	
	
	library(survival)

	if(!grepl("response",dataset))
	{
		genes <- colnames(Y_olp)
		for(gene in genes)
		{
			comb <- cbind(X_olp, Act=Y_olp[,gene])
			
			if(ctype=="OS")
			{
				errflag <- F
				coxmodel_fit <- tryCatch(
					coxph(Surv(OS, OS.Event) ~ ., data = comb),
					
					error = function(e){
					  errflag <<- T
					  },
					
					warning = function(w){
					  errflag <<- T
					  }
				)
				
				smy[gene,dataset] <- ifelse(errflag, NA, summary(coxmodel_fit)$coefficients["Act","z"])
				
			}else if(ctype=="PFS"){
				errflag <- F
				coxmodel_fit <- tryCatch(
					coxph(Surv(PFS, PFS.Event) ~ ., data = comb),
					
					error = function(e){
					  errflag <<- T
					  },
					
					warning = function(w){
					  errflag <<- T
					  }
				)
				
				smy[gene,dataset] <- ifelse(errflag, NA, summary(coxmodel_fit)$coefficients["Act","z"])
				
			}else if(ctype=="RECIST"){
				lm_fit <- lm(RECIST ~ ., data = comb)
				lm_obj <- summary(lm_fit)
				
				smy[gene,dataset] <- -lm_obj$coefficients["Act","t value"]
	
			}else{
				glm_fit <- glm(Response ~ ., data = comb, family = "binomial")
				glm_obj <- summary(glm_fit)
				
				smy[gene,dataset] <- -glm_obj$coefficients["Act","z value"]
			}
		}
	}else{
		
		if(ncol(X_olp)>1)
		{
			if(downsampling==TRUE)
			{
				predictionPath_alt <- paste0(predictionPath,"immunotherapy/",n,"/")
			}else{
				predictionPath_alt <- predictionPath
			}
			
			
			write.table(Y_olp,paste0(predictionPath_alt,"X"),quote=F,sep="\t")
			write.table(X_olp[,"Response",drop=F],paste0(predictionPath_alt,"Y"),quote=F,sep="\t")
			write.table(X_olp[,2:ncol(X_olp),drop=F],paste0(predictionPath_alt,"B"),quote=F,sep="\t")
			system(paste0("/data/rub2/software/logis/bin/logis_batch -B ",predictionPath_alt,"B -X ",predictionPath_alt,"X -Y ",predictionPath_alt,"Y -cntthres 3 -out ",predictionPath_alt,"output"))
			
			zfile <- read.table(paste0(predictionPath_alt,"output.zscore"),sep="\t",check.names=F)	
			
			smy[rownames(zfile),dataset] <- - zfile[,1]
			
			system(paste0("rm ",predictionPath_alt,"X"))
			system(paste0("rm ",predictionPath_alt,"Y"))
			system(paste0("rm ",predictionPath_alt,"B"))
			system(paste0("rm ",predictionPath_alt,"output.coef"))
			system(paste0("rm ",predictionPath_alt,"output.zscore"))
			system(paste0("rm ",predictionPath_alt,"output.pvalue"))
		}else{
			genes <- colnames(Y_olp)
			for(gene in genes)
			{
				comb <- data.frame(Response=as.factor(X_olp[,"Response"]), Act=Y_olp[,gene])
				wilcox_res <- coin::wilcox_test(Act ~ Response, data = comb)
				smy[gene,dataset] = wilcox_res @ statistic @ standardizedlinearstatistic
			}
		}
	
	}
	
}

if(downsampling==TRUE)
{
	write.csv(smy,paste0(predictionPath,"immunotherapy/",n,"/","prediction_act_",i,".csv"),quote=F)
}else{
	write.csv(smy,paste0(predictionPath,group,CTL.corrected,"prediction_act.csv"),quote=F)
}


}

stop("end")



smy <- read.csv(paste0(predictionPath,group,CTL.corrected,"prediction_act.csv"),as.is=T,row.names=1,header=T,check.names=F)

for(vers in c("genomeWide","all"))
{

	if(vers == "genomeWide")
	{
		smy_v <- smy[, setdiff(colnames(smy),notGenomeWide_alt) ]
	}else{
		smy_v <- smy
	}	


	smy.median <- sort(apply(smy_v,1,function(x) median(x,na.rm=T)))
	
	smy_v <- smy_v[names(smy.median),]
	write.csv(smy_v,paste0(predictionPath,group,CTL.corrected,"prediction_act_",vers,".csv"),quote=F)
	
	p.value <- apply(smy_v, 1, function(x) wilcox.test(x,mu=0)$p.value)
	p.adjusted <- p.adjust(p.value,method="BH")
	
	ratio.neg <- apply(smy_v, 1, function(x) {x<-x[!is.na(x)];sum(x< -2)})
	ratio.pos <- apply(smy_v, 1, function(x) {x<-x[!is.na(x)];sum(x>  2)})
	
	is.signif <- (smy.median > 0.5 | smy.median < -0.5) &
		p.value < 0.01 & 
		p.adjusted < 0.01 & 
		( 
			smy.median < 0 & ratio.neg-ratio.pos >=3 | 
			smy.median > 0 & ratio.pos-ratio.neg >=3
		)
		
	smy_test <- cbind(smy.median, p.value=p.value, p.adjusted=p.adjusted, ratio.neg, ratio.pos, is.signif)
	
	write.csv(smy_test,paste0(predictionPath,group,CTL.corrected,"prediction_act_",vers,"_test.csv"),quote=F)
	
	
	smy_top <- as.matrix(smy_v[smy_test[,"is.signif"]==TRUE,])
	
	
	dataCol <- ctype_vec
	geneCol <- ifelse(rownames(smy_top)%in%c("CR1L","AOAH","LY86","COLQ","ADAMTS7"),"green","black")
	dataCol[dataCol=="OS"] <- "black"
	dataCol[dataCol=="PFS"] <- "red"
	dataCol[dataCol=="RECIST"] <- "cyan"
	dataCol[dataCol=="response"] <- "blue"
	
	png(paste0(predictionPath,group,CTL.corrected,"prediction_act_",vers,"_test_signif.png"), width = 25, height = 120, res=200, units = "cm")
	
	library(ComplexHeatmap)
	
	row_ha <- rowAnnotation(
		z = anno_boxplot(smy_top, height = unit(3, "cm") ) 
	)
			
	ht <- Heatmap(smy_top,
		name = "Risk z",
		col = circlize::colorRamp2(c(-4, 0,4), c("blue", "white", "red")),
		column_names_max_height = max_text_width(
	        colnames(smy_top), 
	        gp = gpar(fontsize = 12)
	        ),
	    column_title = paste0(nrow(smy_top)," SPs"),
	    row_names_gp = grid::gpar(fontsize = 9 , col = geneCol ),
	    column_names_gp = grid::gpar(fontsize = 12 , col = dataCol ),
	    right_annotation = row_ha,
		cluster_rows = FALSE,
		cluster_columns = TRUE
	)
	draw(ht)
	
	dev.off()
	
}










smy <- data.frame()
ctype_vec <- c()
for(dataset in datasets)
{
	if(!dataset%in%c("NSCLC_ICB_Ravi2023.Adeno.OS","NSCLC_ICB_Ravi2023.Squamous.OS"))
	{
		study <- strsplit(dataset,".",fixed=T) [[1]][1]
		clini <- strsplit(dataset,".",fixed=T) [[1]][2]
		ctype <- strsplit(clini,"_",fixed=T) [[1]][1]
		surv <- read.csv(gzfile(paste0(dataPath,dataset)),as.is=T,sep="\t",check.names=F)
	}else{
		study <- strsplit(dataset,".O",fixed=T) [[1]][1]
		ctype <- strsplit(dataset,".",fixed=T) [[1]][3]
		surv <- read.csv(gzfile(paste0(dataPath,"NSCLC_ICB_Ravi2023.OS")),as.is=T,sep="\t",check.names=F)
	}
	
	
	if(dataset%in%c(
		"Breast_Durvalumab+Olaparib_Pusztai2021.response",
		"Breast_Pembrolizumab_Wolf2022.response_Basal",
		"Breast_Pembrolizumab_Wolf2022.response_Luminal",
		"HeadNeck_Pembrolizumab_Uppaluri2020.response"))
	{
		colnames(surv)[1] <- "Response"
	}
	if(dataset=="NSCLC_Pembrolizumab_Lee2021.response")
	{
		rownames(surv) <- surv[,1]
		surv <- surv[,-1,drop=F]
	}
	if(ctype%in%c("OS"))
	{
		colnames(surv)[1] <- "OS"
		colnames(surv)[2] <- "OS.Event"	
	}
	if(ctype%in%c("PFS"))
	{
		colnames(surv)[1] <- "PFS"
		colnames(surv)[2] <- "PFS.Event"	
	}
	if(ctype%in%c("RECIST"))
	{
		colnames(surv)[1] <- "RECIST"
	}
	
	ctype_vec <- c(ctype_vec, ctype)
	
	#print(dataset)
	#print(colnames(surv))
	#}
	#print(paste0(colnames(surv)[1]," ",dataset)) 
	#
	#if(colnames(surv)[1]%in%c("Response"))
	#{
	#	print(surv)
	#}
	#}
	
	#Act <- t(as.matrix(read.table(paste0(predictionPath,"immunotherapy/",study,".Zscore"),sep="\t",check.names=F)	))
	Exp <- t(as.matrix(read.table(paste0(predictionPath,"immunotherapy/",study,".diff"),sep="\t",check.names=F)	))
	
	
	#if("SPP1"%in%colnames(Exp)&"CXCL9"%in%colnames(Exp)){
	#print(cor(Exp[,"SPP1"],Exp[,"CXCL9"]))
	#}
	
	
	if(group=="Membrane_")
	{
		matx <- Exp[,colnames(Exp)%in%MPs]
	}else{
		matx <- Exp[,colnames(Exp)%in%colnames(Act)]
	}
	
	olp <- intersect(rownames(surv),rownames(matx))
	X_olp <- surv[olp,,drop=F]
	Y_olp <- matx[olp,,drop=F]
	
	if(CTL.corrected=="CTL.corrected_")
	{
		CTL <- read.table(paste0(predictionPath,"immunotherapy/",study,".CTL"),sep="\t",check.names=F)	
		CTL <- t(as.matrix(CTL))
		X_olp <- cbind(X_olp, CTL=CTL[olp,"CTL"])
	}
	
	library(survival)

	if(!grepl("response",dataset))
	{
		genes <- colnames(Y_olp)
		for(gene in genes)
		{
			comb <- cbind(X_olp, Act=Y_olp[,gene])
			
			if(length(unique(comb[,"Act"]))==1)
			{
				smy[gene,dataset] <- 0
				next
			}
			
			if(ctype=="OS")
			{
				errflag <- F
				coxmodel_fit <- tryCatch(
					coxph(Surv(OS, OS.Event) ~ ., data = comb),
					
					error = function(e){
					  errflag <<- T
					  },
					
					warning = function(w){
					  errflag <<- T
					  }
				)
				
				smy[gene,dataset] <- ifelse(errflag, NA, summary(coxmodel_fit)$coefficients["Act","z"])
				
			}else if(ctype=="PFS"){
				errflag <- F
				coxmodel_fit <- tryCatch(
					coxph(Surv(PFS, PFS.Event) ~ ., data = comb),
					
					error = function(e){
					  errflag <<- T
					  },
					
					warning = function(w){
					  errflag <<- T
					  }
				)
				
				smy[gene,dataset] <- ifelse(errflag, NA, summary(coxmodel_fit)$coefficients["Act","z"])
				
			}else if(ctype=="RECIST"){
				lm_fit <- lm(RECIST ~ ., data = comb)
				lm_obj <- summary(lm_fit)
				
				smy[gene,dataset] <- -lm_obj$coefficients["Act","t value"]
	
			}else{
				glm_fit <- glm(Response ~ ., data = comb, family = "binomial")
				glm_obj <- summary(glm_fit)
				
				smy[gene,dataset] <- -glm_obj$coefficients["Act","z value"]
			}
		}
	}else{
		
		if(ncol(X_olp)>1)
		{
			write.table(Y_olp,paste0("X"),quote=F,sep="\t")
			write.table(X_olp[,"Response",drop=F],paste0("Y"),quote=F,sep="\t")
			write.table(X_olp[,2:ncol(X_olp),drop=F],paste0("B"),quote=F,sep="\t")
			system("/data/rub2/software/logis/bin/logis_batch -B B -X X -Y Y -cntthres 3 -out output")
			
			zfile <- read.table("output.zscore",sep="\t",check.names=F)	
			
			smy[rownames(zfile),dataset] <- - zfile[,1]
			
			system("rm /data/rub2/project/Secretome/code/X")
			system("rm /data/rub2/project/Secretome/code/Y")
			system("rm /data/rub2/project/Secretome/code/B")
			system("rm /data/rub2/project/Secretome/code/output.coef")
			system("rm /data/rub2/project/Secretome/code/output.zscore")
			system("rm /data/rub2/project/Secretome/code/output.pvalue")
		}else{
			genes <- colnames(Y_olp)
			for(gene in genes)
			{
				comb <- data.frame(Response=as.factor(X_olp[,"Response"]), Act=Y_olp[,gene])
				
				if(length(unique(comb[,"Act"]))==1)
				{
					smy[gene,dataset] <- 0
					next
				}
				
				wilcox_res <- coin::wilcox_test(Act ~ Response, data = comb)
				smy[gene,dataset] = wilcox_res @ statistic @ standardizedlinearstatistic
			}
		}
	
	}
	
}

write.csv(smy,paste0(predictionPath,group,CTL.corrected,"prediction_exp.csv"),quote=F)


smy <- read.csv(paste0(predictionPath,group,CTL.corrected,"prediction_exp.csv"),as.is=T,row.names=1,header=T,check.names=F)

for(vers in c("genomeWide","all"))
{

	if(vers == "genomeWide")
	{
		smy_v <- smy[, setdiff(colnames(smy),notGenomeWide_alt) ]
		smy_v <- smy_v[apply(smy_v, 1, function(x) sum(is.na(x)))<2,]
	}else{
		smy_v <- smy
	}	


	smy.median <- sort(apply(smy_v,1,function(x) median(x,na.rm=T)))
	
	smy_v <- smy_v[names(smy.median),]
	write.csv(smy_v,paste0(predictionPath,group,CTL.corrected,"prediction_exp_",vers,".csv"),quote=F)
	
	p.value <- apply(smy_v, 1, function(x) wilcox.test(x,mu=0)$p.value)
	p.adjusted <- p.adjust(p.value,method="BH")
	
	ratio.neg <- apply(smy_v, 1, function(x) {x<-x[!is.na(x)];sum(x< -2)})
	ratio.pos <- apply(smy_v, 1, function(x) {x<-x[!is.na(x)];sum(x>  2)})
	
	is.signif <- (smy.median > 0.5 | smy.median < -0.5) &
		p.value < 0.01 & 
		p.adjusted < 0.01 & 
		( 
			smy.median < 0 & ratio.neg-ratio.pos >=3 | 
			smy.median > 0 & ratio.pos-ratio.neg >=3
			
		)
		
	smy_test <- cbind(smy.median, p.value=p.value, p.adjusted=p.adjusted, ratio.neg, ratio.pos, is.signif)
	
	write.csv(smy_test,paste0(predictionPath,group,CTL.corrected,"prediction_exp_",vers,"_test.csv"),quote=F)
	
	
	smy_top <- as.matrix(smy_v[smy_test[,"is.signif"]==TRUE,])
	
	
	dataCol <- ctype_vec
	geneCol <- ifelse(rownames(smy_top)%in%c("CR1L","AOAH","LY86","COLQ","ADAMTS7"),"green","black")
	dataCol[dataCol=="OS"] <- "black"
	dataCol[dataCol=="PFS"] <- "red"
	dataCol[dataCol=="RECIST"] <- "cyan"
	dataCol[dataCol=="response"] <- "blue"
	
	png(paste0(predictionPath,group,CTL.corrected,"prediction_exp_",vers,"_test_signif.png"), width = 25, height = 120, res=200, units = "cm")
	
	library(ComplexHeatmap)
	
	row_ha <- rowAnnotation(
		z = anno_boxplot(smy_top, height = unit(3, "cm") ) 
	)
			
	ht <- Heatmap(smy_top,
		name = "Risk z",
		col = circlize::colorRamp2(c(-4, 0,4), c("blue", "white", "red")),
		column_names_max_height = max_text_width(
	        colnames(smy_top), 
	        gp = gpar(fontsize = 12)
	        ),
	    column_title = paste0(nrow(smy_top)," SPs"),
	    row_names_gp = grid::gpar(fontsize = 9 , col = geneCol ),
	    column_names_gp = grid::gpar(fontsize = 12 , col = dataCol ),
	    right_annotation = row_ha,
		cluster_rows = FALSE,
		cluster_columns = TRUE
	)
	draw(ht)
	
	dev.off()
	
}










































smy <- data.frame()
smy2 <- data.frame()
ctype_vec <- c()
for(dataset in datasets)
{
	if(!dataset%in%c("NSCLC_ICB_Ravi2023.Adeno.OS","NSCLC_ICB_Ravi2023.Squamous.OS"))
	{
		study <- strsplit(dataset,".",fixed=T) [[1]][1]
		clini <- strsplit(dataset,".",fixed=T) [[1]][2]
		ctype <- strsplit(clini,"_",fixed=T) [[1]][1]
		surv <- read.csv(gzfile(paste0(dataPath,dataset)),as.is=T,sep="\t",check.names=F)
	}else{
		study <- strsplit(dataset,".O",fixed=T) [[1]][1]
		ctype <- strsplit(dataset,".",fixed=T) [[1]][3]
		surv <- read.csv(gzfile(paste0(dataPath,"NSCLC_ICB_Ravi2023.OS")),as.is=T,sep="\t",check.names=F)
	}
	
	
	if(dataset%in%c(
		"Breast_Durvalumab+Olaparib_Pusztai2021.response",
		"Breast_Pembrolizumab_Wolf2022.response_Basal",
		"Breast_Pembrolizumab_Wolf2022.response_Luminal",
		"HeadNeck_Pembrolizumab_Uppaluri2020.response"))
	{
		colnames(surv)[1] <- "Response"
	}
	if(dataset=="NSCLC_Pembrolizumab_Lee2021.response")
	{
		rownames(surv) <- surv[,1]
		surv <- surv[,-1,drop=F]
	}
	if(ctype%in%c("OS"))
	{
		colnames(surv)[1] <- "OS"
		colnames(surv)[2] <- "OS.Event"	
	}
	if(ctype%in%c("PFS"))
	{
		colnames(surv)[1] <- "PFS"
		colnames(surv)[2] <- "PFS.Event"	
	}
	if(ctype%in%c("RECIST"))
	{
		colnames(surv)[1] <- "RECIST"
	}
	
	ctype_vec <- c(ctype_vec, ctype)
	
	#print(dataset)
	#print(colnames(surv))
	#}
	#print(paste0(colnames(surv)[1]," ",dataset)) 
	#
	#if(colnames(surv)[1]%in%c("Response"))
	#{
	#	print(surv)
	#}
	#}
	
	load(paste0(predictionPath,"immunotherapy/",study,".RData"))
	Act <- as.matrix(res$zscore)
	Act <- expand_rows(Act)
	Act <- t(Act)
	
	
	Exp <- t(as.matrix(read.table(paste0(predictionPath,"immunotherapy/",study,".txt"),sep="\t",check.names=F)	))
	
	
	data_MSigDB_path <- "/data/rub2/data/MSigDB/"
	gmtName <- "h.all.v2023.2.Hs.symbols.gmt"
	gmt <- read.gmt(paste0(data_MSigDB_path,gmtName))
	gmt <- gmt[c("HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_INTERFERON_ALPHA_RESPONSE","HALLMARK_INTERFERON_GAMMA_RESPONSE")]
	
	gsva_scores <- GSVA::gsva(t(Exp), gmt, method = "gsva", kcdf = "Gaussian")
	rownames(gsva_scores) <- c("NFKB_Sig","IFNA_Sig","IFNG_Sig")
	
	
	# https://amigo.geneontology.org/amigo/term/GO:0034142
	#TLR4Sig <- c("TRAF3","TLR4","MYD88","TICAM1","ECSIT","PIK3AP1","TNIP3","NMI","RAB13","SQSTM1","TRIM32","CHUK","IRAK4","OAS1","LYN","NAGLU","TIRAP","CAV1","RNF115","TRAF6","ZNRF1","LETMD1","TMEM126A","TBK1","LY96","GDI1","IRAK1","CD14","AKT1","RELA","YWHAE","MAP3K7","IRAK2","RIPK2","IRF3","NFKBIA","TICAM2","S100A14","TRIL","SCIMP","TAX1BP1")
	#Exp_TLR4Sig <- Exp[,colnames(Exp)%in%TLR4Sig, drop=F]
	#Exp <- cbind(Exp,TLR4Sig=rowMeans(Exp_TLR4Sig))
	
	colnames(Exp) <- paste0(colnames(Exp),"_exp")
	Exp <- cbind(Exp,t(gsva_scores) )
	Exp <- Exp[,colnames(Exp)%in%c("LY86_exp","CD180_exp","NFKB_Sig","IFNA_Sig","IFNG_Sig"),drop=F]	
	
	matx <- cbind(LY86_act=Act[,"LY86"], Exp)
	
	smy2[colnames(Exp),dataset] <- cor(Act[,"LY86",drop=F], Exp)[1,]
	
	
	olp <- intersect(rownames(surv),rownames(matx))
	X_olp <- surv[olp,,drop=F]
	Y_olp <- matx[olp,,drop=F]
	
	if(CTL.corrected=="CTL.corrected_")
	{
		CTL <- read.table(paste0(predictionPath,"immunotherapy/",study,".CTL"),sep="\t",check.names=F)	
		CTL <- t(as.matrix(CTL))
		X_olp <- cbind(X_olp, CTL=CTL[olp,"CTL"])
	}
	
	library(survival)

	if(!grepl("response",dataset))
	{
		genes <- colnames(Y_olp)
		for(gene in genes)
		{
			comb <- cbind(X_olp, Act=Y_olp[,gene])
			
			if(length(unique(comb[,"Act"]))==1)
			{
				smy[gene,dataset] <- 0
				next
			}
			
			if(ctype=="OS")
			{
				errflag <- F
				coxmodel_fit <- tryCatch(
					coxph(Surv(OS, OS.Event) ~ ., data = comb),
					
					error = function(e){
					  errflag <<- T
					  },
					
					warning = function(w){
					  errflag <<- T
					  }
				)
				
				smy[gene,dataset] <- ifelse(errflag, NA, summary(coxmodel_fit)$coefficients["Act","z"])
				
			}else if(ctype=="PFS"){
				errflag <- F
				coxmodel_fit <- tryCatch(
					coxph(Surv(PFS, PFS.Event) ~ ., data = comb),
					
					error = function(e){
					  errflag <<- T
					  },
					
					warning = function(w){
					  errflag <<- T
					  }
				)
				
				smy[gene,dataset] <- ifelse(errflag, NA, summary(coxmodel_fit)$coefficients["Act","z"])
				
			}else if(ctype=="RECIST"){
				lm_fit <- lm(RECIST ~ ., data = comb)
				lm_obj <- summary(lm_fit)
				
				smy[gene,dataset] <- -lm_obj$coefficients["Act","t value"]
	
			}else{
				glm_fit <- glm(Response ~ ., data = comb, family = "binomial")
				glm_obj <- summary(glm_fit)
				
				smy[gene,dataset] <- -glm_obj$coefficients["Act","z value"]
			}
		}
	}else{
		
		if(ncol(X_olp)>1)
		{
			write.table(Y_olp,paste0("X"),quote=F,sep="\t")
			write.table(X_olp[,"Response",drop=F],paste0("Y"),quote=F,sep="\t")
			write.table(X_olp[,2:ncol(X_olp),drop=F],paste0("B"),quote=F,sep="\t")
			system("/data/rub2/software/logis/bin/logis_batch -B B -X X -Y Y -cntthres 3 -out output")
			
			zfile <- read.table("output.zscore",sep="\t",check.names=F)	
			
			smy[rownames(zfile),dataset] <- - zfile[,1]
			
			system("rm /data/rub2/project/Secretome/code/X")
			system("rm /data/rub2/project/Secretome/code/Y")
			system("rm /data/rub2/project/Secretome/code/B")
			system("rm /data/rub2/project/Secretome/code/output.coef")
			system("rm /data/rub2/project/Secretome/code/output.zscore")
			system("rm /data/rub2/project/Secretome/code/output.pvalue")
		}else{
			genes <- colnames(Y_olp)
			for(gene in genes)
			{
				comb <- data.frame(Response=as.factor(X_olp[,"Response"]), Act=Y_olp[,gene])
				
				if(length(unique(comb[,"Act"]))==1)
				{
					smy[gene,dataset] <- 0
					next
				}
				
				wilcox_res <- coin::wilcox_test(Act ~ Response, data = comb)
				smy[gene,dataset] = wilcox_res @ statistic @ standardizedlinearstatistic
			}
		}
	
	}
	
}


smy <- as.matrix(smy[, setdiff(colnames(smy),notGenomeWide_alt) ])
smy <- smy[1:4,]

fg.df <- reshape::melt(smy)

sigOrder <- data.frame()
for(i in unique(fg.df[,"X1"]))
{
	sigOrder[i,"r_mean"] <- median(fg.df[fg.df[,"X1"]==i,"value"],na.rm=TRUE)
}
sigOrder <- sigOrder[order(sigOrder[,"r_mean"],decreasing=T),,drop=F]

fg.df[,1] <- factor(fg.df[,1], levels=rownames(sigOrder))


library(ggplot2)
p1 <- ggplot(fg.df,aes(x=X1,y=value)) + 
	geom_hline(yintercept=0, color = "grey", linewidth=0.6, linetype="dashed")+
	#geom_violin(aes(group=X1,fill=X1), trim=FALSE, alpha=0.3, width=0.8)+
	geom_jitter(color="darkgrey",alpha=0.3, size=1, width=0.1)+
	geom_boxplot(color="black", alpha=0, width=0.3, outlier.shape = NA)+
	ylab("Risk score z")+
	xlab("")+
	theme_classic()+ 
	theme(
		panel.background = element_blank(),
		panel.grid = element_blank(),
		axis.title = element_text(colour = "black"),
		axis.text = element_text(colour = "black"),
				  axis.text.x = element_text(size=8, angle = 90, hjust = 1, vjust = 0.5),
		#axis.ticks = element_line(colour = "darkgrey"),
		#axis.line = element_line(colour = "darkgrey"),
		axis.title.x = element_blank(),
		legend.position="none"
	)
	
ggsave(paste0(predictionPath,"prediction_LY86_NFKB_z.png"), p1, width = 5, height = 4.8, dpi=500, units = "cm")




smy2 <- as.matrix(smy2[, setdiff(colnames(smy2),notGenomeWide_alt) ])
smy2 <- smy2[1:3,]

fg.df <- reshape::melt(smy2)

sigOrder <- data.frame()
for(i in unique(fg.df[,"X1"]))
{
	sigOrder[i,"r_mean"] <- median(fg.df[fg.df[,"X1"]==i,"value"],na.rm=TRUE)
}
sigOrder <- sigOrder[order(sigOrder[,"r_mean"],decreasing=T),,drop=F]

fg.df[,1] <- factor(fg.df[,1], levels=rownames(sigOrder))



library(ggplot2)
p1 <- ggplot(fg.df,aes(x=X1,y=value)) + 
	geom_hline(yintercept=0, color = "grey", linewidth=0.6, linetype="dashed")+
	geom_jitter(color="darkgrey",alpha=0.3, size=1, width=0.1)+
	geom_boxplot(color="black", alpha=0, width=0.5, outlier.shape = NA)+
	ylab("Cor")+
	xlab("")+
	theme_classic()+ 
	theme(
		panel.background = element_blank(),
		panel.grid = element_blank(),
		axis.title = element_text(colour = "black"),
		axis.text = element_text(colour = "black"),
		#axis.ticks = element_line(colour = "darkgrey"),
		#axis.line = element_line(colour = "darkgrey"),
		axis.title.x = element_blank(),
		legend.position="none"
	)
	
ggsave(paste0(predictionPath,"prediction_LY86_NFKB_cor.png"), p1, width = 10, height = 4.8, dpi=500, units = "cm")














write.csv(smy,paste0(predictionPath,group,CTL.corrected,"prediction_exp_LY86.csv"),quote=F)


smy <- read.csv(paste0(predictionPath,group,CTL.corrected,"prediction_exp.csv"),as.is=T,row.names=1,header=T,check.names=F)

for(vers in c("genomeWide","all"))
{

	if(vers == "genomeWide")
	{
		smy_v <- smy[, setdiff(colnames(smy),notGenomeWide_alt) ]
		smy_v <- smy_v[apply(smy_v, 1, function(x) sum(is.na(x)))<2,]
	}else{
		smy_v <- smy
	}	


	smy.median <- sort(apply(smy_v,1,function(x) median(x,na.rm=T)))
	
	smy_v <- smy_v[names(smy.median),]
	write.csv(smy_v,paste0(predictionPath,group,CTL.corrected,"prediction_exp_",vers,".csv"),quote=F)
	
	p.value <- apply(smy_v, 1, function(x) wilcox.test(x,mu=0)$p.value)
	p.adjusted <- p.adjust(p.value,method="BH")
	
	ratio.neg <- apply(smy_v, 1, function(x) {x<-x[!is.na(x)];sum(x< -2)})
	ratio.pos <- apply(smy_v, 1, function(x) {x<-x[!is.na(x)];sum(x>  2)})
	
	is.signif <- (smy.median > 0.5 | smy.median < -0.5) &
		p.value < 0.01 & 
		p.adjusted < 0.01 & 
		( 
			smy.median < 0 & ratio.neg-ratio.pos >=3 | 
			smy.median > 0 & ratio.pos-ratio.neg >=3
			
		)
		
	smy_test <- cbind(smy.median, p.value=p.value, p.adjusted=p.adjusted, ratio.neg, ratio.pos, is.signif)
	
	write.csv(smy_test,paste0(predictionPath,group,CTL.corrected,"prediction_exp_",vers,"_test.csv"),quote=F)
	
	
	smy_top <- as.matrix(smy_v[smy_test[,"is.signif"]==TRUE,])
	
	
	dataCol <- ctype_vec
	geneCol <- ifelse(rownames(smy_top)%in%c("CR1L","AOAH","LY86","COLQ","ADAMTS7"),"green","black")
	dataCol[dataCol=="OS"] <- "black"
	dataCol[dataCol=="PFS"] <- "red"
	dataCol[dataCol=="RECIST"] <- "cyan"
	dataCol[dataCol=="response"] <- "blue"
	
	png(paste0(predictionPath,group,CTL.corrected,"prediction_exp_",vers,"_test_signif.png"), width = 25, height = 120, res=200, units = "cm")
	
	library(ComplexHeatmap)
	
	row_ha <- rowAnnotation(
		z = anno_boxplot(smy_top, height = unit(3, "cm") ) 
	)
			
	ht <- Heatmap(smy_top,
		name = "Risk z",
		col = circlize::colorRamp2(c(-4, 0,4), c("blue", "white", "red")),
		column_names_max_height = max_text_width(
	        colnames(smy_top), 
	        gp = gpar(fontsize = 12)
	        ),
	    column_title = paste0(nrow(smy_top)," SPs"),
	    row_names_gp = grid::gpar(fontsize = 9 , col = geneCol ),
	    column_names_gp = grid::gpar(fontsize = 12 , col = dataCol ),
	    right_annotation = row_ha,
		cluster_rows = FALSE,
		cluster_columns = TRUE
	)
	draw(ht)
	
	dev.off()
	
}






















######
# compare significant SP count from act and expr
######

smy_test1 <- read.csv(paste0(predictionPath,"prediction_act_genomeWide_test.csv"),row.names=1,header=T)
smy_test2 <- read.csv(paste0(predictionPath,"prediction_exp_genomeWide_test.csv"),row.names=1,header=T)

genes_act <- rownames(smy_test1)[smy_test1[,"is.signif"]==1]
genes_exp <- rownames(smy_test2)[smy_test2[,"is.signif"]==1]

olp <- intersect(genes_act,genes_exp)

length(genes_act)
length(genes_exp)
length(olp)


write.csv(smy_test1[genes_act,], paste0(predictionPath,"prediction_act_genomeWide_test_signif.csv"),quote=F)
write.csv(smy_test2[genes_exp,], paste0(predictionPath,"prediction_exp_genomeWide_test_signif.csv"),quote=F)

write.csv(smy_test1[rownames(smy_test1)%in%olp,], paste0(predictionPath,"prediction_act_genomeWide_test_signif_overlap.csv"),quote=F)
write.csv(smy_test2[rownames(smy_test2)%in%olp,], paste0(predictionPath,"prediction_exp_genomeWide_test_signif_overlap.csv"),quote=F)



genes_act_neg <- rownames(smy_test1)[smy_test1[,"is.signif"]==1&smy_test1[,"smy.median"]<0]
genes_act_pos <- rev(rownames(smy_test1)[smy_test1[,"is.signif"]==1&smy_test1[,"smy.median"]>0])

genes_exp_neg <- rownames(smy_test2)[smy_test2[,"is.signif"]==1&smy_test2[,"smy.median"]<0]
genes_exp_pos <- rev(rownames(smy_test2)[smy_test2[,"is.signif"]==1&smy_test2[,"smy.median"]>0])

genes_olp_pos <- intersect(genes_act_pos,genes_exp_pos)
genes_olp_neg <- intersect(genes_act_neg,genes_exp_neg)


length(genes_act_neg)
length(genes_act_pos)

length(genes_exp_neg)
length(genes_exp_pos)

length(intersect(genes_act_neg,genes_exp_neg))
length(intersect(genes_act_pos,genes_exp_pos))
length(intersect(genes_act_neg,genes_exp_pos))
length(intersect(genes_act_pos,genes_exp_neg))



fg.df <- data.frame(
	group=c("Activity","Activity","Expression","Expression"),
	mod=c("Anti","Pro","Anti","Pro"),
	Count=c(
		length(genes_act_neg),
		length(genes_act_pos),
		length(genes_exp_neg),
		length(genes_exp_pos)
	)
)

library(ggplot2)
p1 <- ggplot(fg.df, aes(x=group, y=Count, fill=mod)) +
  geom_bar(stat="identity", position="dodge", color="white", width=0.4, alpha=0.7) +
  scale_fill_manual( values=c("#4d9221","#c51b7d") )+
  ylab("# Secreted Proteins")+
  theme_classic()+
  theme(
	  panel.grid = element_blank(),
	  panel.background = element_blank(),
	  axis.text = element_text(colour = "black"),
	  axis.title = element_text(colour = "black"),
	  axis.title.x = element_blank(),
	  legend.key.size = unit(0.4, 'cm'),
	  legend.position = c(0.75,0.8),
	  legend.title = element_blank()
	)

ggsave(paste0(predictionPath,"prediction_act_exp_count_compare.pdf"), p1, width = 5.3, height = 4.8, units = "cm")






# same anti and pro number
genes1_top <- c(genes_act_neg[1:length(genes_exp_neg)],rev(genes_act_pos[1:length(genes_exp_pos)]))


# same total number
#smy_test1_alt <- smy_test1
#smy_test1_alt[,1] <- abs(smy_test1_alt[,1])
#smy_test1_alt <- smy_test1_alt[order(smy_test1_alt[,1],decreasing=T),]
#genes1_top_alt <- rownames(smy_test1_alt)[1:length(genes_exp)]


to_be_annotated <- unique(c(genes1_top, genes_exp, genes_olp_pos, genes_olp_neg))




library(rio) 
xlsx <- import_list(paste0(predictionPath,"soluble_factors.doubleblind.Lanqi.xlsx")) 

for(n in 1:length(xlsx))
{
	if(n==1)
	{
		comb <- xlsx[[n]]
	}else{
		comb <- rbind(comb,xlsx[[n]])
	}
}

rownames(comb) <- comb[,1]
comb[is.na(comb[,2]),2] <- "Unknown"


setdiff(to_be_annotated,comb[,1])




write.csv(comb[genes1_top, 2:3], paste0(predictionPath,"prediction_act_genomeWide_test_signif_top_annotation.csv"))
write.csv(comb[genes_exp, 2:3], paste0(predictionPath,"prediction_exp_genomeWide_test_signif_annotation.csv"))


smy_test1 <- smy_test1[rev(rownames(smy_test1)),]
smy_test2 <- smy_test2[rev(rownames(smy_test2)),]
genes_olp <- c( genes_olp_neg, rev(genes_olp_pos) )


lqlist <- read.csv(paste0(predictionPath,"list_annotation.csv"),row.names=1,header=T)
lqlist[lqlist[,1]=="Anti-tumor (verified in-house)",1] <- "Anti-tumor"
lqlist[lqlist[,1]=="Pro-tumor (verified in-house)",1] <- "Pro-tumor"

is.conflict1 <- cbind(annotation=comb[genes_olp_neg, 2, drop=F], prediction="Anti-tumor")
is.conflict2 <- cbind(annotation=comb[genes_olp_pos, 2, drop=F], prediction="Pro-tumor")
is.conflict <- rbind(is.conflict1, is.conflict2)
#is.conflict <- is.conflict[!is.conflict[,1]==is.conflict[,2],]

is.conflict <- is.conflict[!rownames(is.conflict)%in%rownames(lqlist),]


olp_annotation <- comb[genes_olp, 2:3]
olp_annotation <- cbind(olp_annotation, followup="")

for(gene in rownames(olp_annotation))
{
	if(!gene%in%rownames(lqlist))
	{
		print(gene)
		next
	}
	
	if(olp_annotation[gene,1]==lqlist[gene,1])
	{
		next
	}else{
		olp_annotation[gene,"followup"] <- lqlist[gene,1]
	}
}
olp_annotation <- olp_annotation[,c(1,3,2)]

write.csv(olp_annotation, paste0(predictionPath,"prediction_olp_genomeWide_test_signif_annotation.csv"))





comb <- comb[comb[,2]%in%c("Anti-tumor","Pro-tumor"),]

smy <- data.frame()

for(i in c("act","exp"))
{
	for(j in c("top","olp","unique"))
	{
		if(i=="act"&j=="top")
		{
			genes <- genes1_top
			value <- smy_test1[genes,"smy.median"]
		}
		if(i=="act"&j=="olp")
		{
			genes <- intersect(genes_act,genes_exp)
			value <- smy_test1[genes,"smy.median"]
		}
		if(i=="act"&j=="unique")
		{
			genes <- setdiff(genes_act,genes_exp)
			value <- smy_test1[genes,"smy.median"]
		}
		if(i=="exp"&j=="top")
		{
			genes <- genes_exp
			value <- smy_test2[genes,"smy.median"]
		}
		if(i=="exp"&j=="olp")
		{
			genes <- intersect(genes_act,genes_exp)
			value <- smy_test2[genes,"smy.median"]
		}
		if(i=="exp"&j=="unique")
		{
			genes <- setdiff(genes_exp,genes_act)
			value <- smy_test1[genes,"smy.median"]
		}
		
		
		
		fg.df <- data.frame(
			value=value,
			event=NA
		)
		rownames(fg.df) <- genes
		
		
		olp <- intersect(genes,rownames(comb))
		
		
		fg.df <- fg.df[olp,]
		comb_olp <- comb[olp,,drop=FALSE]
		
		
		fg.df[comb_olp[,2]%in%c("Pro-tumor"),"event"] <- 1
		fg.df[comb_olp[,2]%in%c("Anti-tumor"),"event"] <- 0
		
		fg.df <- fg.df[!is.na(fg.df[,"event"]),]
		
		smy[paste0(i,"_",j),"n"] <- dim(fg.df)[1]
		smy[paste0(i,"_",j),"F1base"] <- sum(fg.df[,"event"])/dim(fg.df)[1]
		
		
		predicted <- factor(as.numeric(fg.df[,1]>0),levels = c(1,0)) # 1:event; 0:nonEvent
		actual <- factor(as.numeric(fg.df[,2]==1),levels = c(1,0))
		
		xtab <- table(predicted, actual)
 
		library(caret)
		cm <- caret::confusionMatrix(xtab)
		
		TP <- cm$table[1,1]
		TN <- cm$table[2,2]
		FP <- cm$table[1,2]
		FN <- cm$table[2,1]
		
		
		smy[paste0(i,"_",j),"n"] <- dim(fg.df)[1]
		smy[paste0(i,"_",j),"Accuracy"] <- cm[["overall"]]["Accuracy"]
		smy[paste0(i,"_",j),"Sensitivity"] <- cm[["byClass"]]["Sensitivity"]
		smy[paste0(i,"_",j),"Specificity"] <- cm[["byClass"]]["Specificity"]
		smy[paste0(i,"_",j),"Precision"] <- cm[["byClass"]]["Precision"]
		smy[paste0(i,"_",j),"Recall"] <- cm[["byClass"]]["Recall"]
		smy[paste0(i,"_",j),"F1"] <- cm[["byClass"]]["F1"]
		smy[paste0(i,"_",j),"MCC"] <- (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))

		
		library(ROCR)
		predM <- prediction(fg.df[,1], fg.df[,"event"])
		roc = performance(predM, measure = "auc")
		roc_plot = performance(predM, measure = "tpr", x.measure = "fpr")
		
		
		png(paste0(predictionPath,i,"_",j,"_ROC.png"), width = 8, height = 9.5, res=400, units = "cm")
		
		par(mfrow= c(1, 1))
		plot(roc_plot,main="",col=c("blue","red")) 
		abline(a=0, b=1, lty=2)
		
		dev.off()
		
		
		AUC <- round(roc@ y.values [[1]],2)
		
		smy[paste0(i,"_",j),"AUC"] <- AUC

		pro <- sum(fg.df[,"event"])
		anti <- dim(fg.df)[1] - pro
		
		fg.df[fg.df[,"event"]==1,"event"] <- paste0("Pro-tumor (",pro,")")
		fg.df[fg.df[,"event"]=="0","event"] <- paste0("Anti-tumor (",anti,")")
		
		library(ggplot2)
		p2 <- ggplot(fg.df,aes(x=event,y=value)) + 
				geom_hline(yintercept=0, color = "grey", linewidth=0.8, linetype="dashed")+
				geom_boxplot( colour="grey1", alpha=0, outlier.shape = NA)+
				geom_jitter(alpha=0.2,size=1)+
				annotate("text", x = 1, y=0.8, label = paste0("AUC = ",AUC))+
				ggtitle(paste0(i," ",j))+
				ylab("Median Risk z")+
				xlab("")+
				theme_bw()+ 
				theme(
					panel.background = element_blank(),
					panel.grid = element_blank(),
					axis.title = element_text(colour = "black"),
					axis.text = element_text(colour = "black"),
					legend.position="none"
				)
		ggsave(paste0(predictionPath,i,"_",j,"_AUC.png"), p2, width = 8, height = 8, dpi=200, units = "cm",limitsize = FALSE)

	}
}

smy <- smy[c(1,4,2,5,3,6),]

write.csv(smy, paste0(predictionPath,"performance_compare_act_exp.csv"))


fg.df <- data.frame(
	x=c("ActExp","ActExp","ActExp","ActExp","ActExp","ActExp","ActExp","ActExp"),
	group=c("Activity","Activity","Activity","Activity","Expression","Expression","Expression","Expression"),
	mod=c("Accuracy","Sensitivity","Specificity","MCC","Accuracy","Sensitivity","Specificity","MCC"),
	Value=c(unlist(smy[1,c(3:5,9)]),unlist(smy[2,c(3:5,9)]))
)
fg.df[,3] <- factor(fg.df[,3],levels=c("Sensitivity","Specificity","Accuracy","MCC"))


library(ggplot2)
p0 <- ggplot(fg.df, aes(x=group, y=Value, fill=mod)) +
  #geom_hline(yintercept=0.5, color = "grey", linewidth=0.6, linetype="dashed")+
  geom_bar(stat="identity", position="dodge", color="white", width=0.66, alpha=0.7) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE,size=3))+
  scale_y_continuous(breaks = seq(0, 1, by = 0.2))+
  theme_classic()+
  theme(
	  panel.grid = element_blank(),
	  panel.background = element_blank(),
	  axis.text = element_text(colour = "black"),
	  axis.title = element_text(colour = "black"),
	  axis.title.x = element_blank(),
	  strip.background = element_rect(colour="white", fill="white"),
	  legend.key.size = unit(0.3, 'cm'),
	  legend.position = "bottom",
	  legend.title = element_blank(),
	  legend.text=element_text(size=8)
	)
	
ggsave(paste0(predictionPath,"prediction_act_exp_accuracy_compare.png"), p0, width = 5.3, height = 7, dpi=400, units = "cm")




fg.df1 <- fg.df[fg.df[,3]%in%c("Precision","Recall"),]
fg.df2 <- fg.df[fg.df[,3]%in%c("Accuracy","F1"),]

library(ggplot2)
p1 <- ggplot(fg.df1, aes(x=x, y=Value, fill=group)) +
  geom_bar(stat="identity", position="dodge", color="white", width=0.66, alpha=0.7) +
  scale_fill_manual( values=c("#8c510a","#dfc27d") )+
  theme_classic()+
  theme(
	  panel.grid = element_blank(),
	  panel.background = element_blank(),
	  axis.text = element_text(colour = "black"),
	  axis.title = element_text(colour = "black"),
	  axis.text.x = element_blank(),
	  axis.title.x = element_blank(),
	  axis.ticks.x = element_blank(),
	  strip.background = element_rect(colour="white", fill="white"),
	  legend.key.size = unit(0.4, 'cm'),
	  legend.position = c(0.75,0.8),
	  legend.title = element_blank()
	)+ facet_wrap(~ mod, ncol =2)+
	coord_cartesian(ylim = c(0.2, 1))
  
ggsave(paste0(predictionPath,"prediction_act_exp_accuracy_compare_sub1.png"), p1, width = 5, height = 5, dpi=400, units = "cm")


library(ggplot2)
p2 <- ggplot(fg.df2, aes(x=x, y=Value, fill=group)) +
  geom_bar(stat="identity", position="dodge", color="white", width=0.66, alpha=0.7) +
  scale_fill_manual( values=c("#8c510a","#dfc27d") )+
  theme_classic()+
  theme(
	  panel.grid = element_blank(),
	  panel.background = element_blank(),
	  axis.text = element_text(colour = "black"),
	  axis.title = element_text(colour = "black"),
	  axis.text.x = element_blank(),
	  axis.title.x = element_blank(),
	  axis.ticks.x = element_blank(),
	  strip.background = element_rect(colour="white", fill="white"),
	  legend.position = "none"
	)+ facet_wrap(~ mod, ncol =2)+
	coord_cartesian(ylim = c(0.5, 0.8))
  
ggsave(paste0(predictionPath,"prediction_act_exp_accuracy_compare_sub2.png"), p2, width = 5, height = 5, dpi=400, units = "cm")









################ 
# long heatmap
################


smy1 <- read.csv(paste0(predictionPath,"prediction_act.csv"),row.names=1,header=T,check.names=F)
smy2 <- read.csv(paste0(predictionPath,"prediction_exp.csv"),row.names=1,header=T,check.names=F)

smy1_l <- smy1[, setdiff(colnames(smy1),notGenomeWide_alt) ]
smy1_s <- smy1[, notGenomeWide_alt]

smy2_l <- smy2[, setdiff(colnames(smy2),notGenomeWide_alt) ]
smy2_s <- smy2[, notGenomeWide_alt]

smy_l <- cbind(smy1_l, smy2_l[rownames(smy1_l),])
smy_s <- cbind(smy1_s, smy2_s[rownames(smy1_s),])

smy_l.median <- apply(smy_l,1,function(x) median(x,na.rm=T))
smy_s.median <- apply(smy_s,1,function(x) median(x,na.rm=T))

sorted.by.act.exp <- names(sort(smy_l.median))


olp <- intersect(genes_act,genes_exp)
#olp <- sorted.by.act.exp[sorted.by.act.exp%in%olp]
olp_neg <- intersect(olp,genes_act_neg)




smy1 <- read.csv(paste0(predictionPath,"prediction_act.csv"),row.names=1,header=T,check.names=F)
smy2 <- read.csv(paste0(predictionPath,"prediction_exp.csv"),row.names=1,header=T,check.names=F)

rename <- read.csv(paste0(predictionPath,"rename.txt"),row.names=1,header=T,sep="\t")

colnames(smy1) <- rename[match(colnames(smy1),rename[,1]),2]
colnames(smy2) <- rename[match(colnames(smy2),rename[,1]),2]

smy1 <- smy1[rownames(smy1)%in%olp,]
smy2 <- smy2[rownames(smy2)%in%olp,]

smy1 <- smy1[rev(olp),]
smy2 <- smy2[rev(olp),]

smy1.m <- reshape2::melt(as.matrix(smy1))
smy2.m <- reshape2::melt(as.matrix(smy2))

smy <- rbind(cbind(smy1.m,group="Act"),cbind(smy2.m,group="Exp"))

smy_l <- smy[!smy[,2]%in%notGenomeWide_alt2,]
smy_s <- smy[smy[,2]%in%notGenomeWide_alt2,]

smy_l_sorted <- sort(apply(smy1[olp_neg,unique(smy_l[,2])],2,function(x) mean(x,na.rm=T)))
smy_s_sorted <- sort(apply(smy1[olp_neg,unique(smy_s[,2])],2,function(x) mean(x,na.rm=T)))


library(tidyverse)

mydat <- smy_l
colnames(mydat) <- c("Species","month","measurement","flower_att")
mydat[,2] <- as.character(mydat[,2])


make_triangles <- function(x, y, point = "up") {
  x <- as.integer(as.factor((x)))
  y <- as.integer(as.factor((y)))

  if (point == "up") {
    newx <- sapply(x, function(x) {
      c(x - 0.5, x - 0.5, x + 0.5)
    }, simplify = FALSE)
    newy <- sapply(y, function(y) {
      c(y - 0.5, y + 0.5, y + 0.5)
    }, simplify = FALSE)
  } else if (point == "down") {
    newx <- sapply(x, function(x) {
      c(x - 0.5, x + 0.5, x + 0.5)
    }, simplify = FALSE)
    newy <- sapply(y, function(y) {
      c(y - 0.5, y - 0.5, y + 0.5)
    }, simplify = FALSE)
  }
  data.frame(x = unlist(newx), y = unlist(newy))
}

# required, otherwise you cannot use the values as fill
mydat_wide <- mydat %>% pivot_wider(names_from = "flower_att", values_from = "measurement")
# making your ordered months factor
mydat_wide$month <- factor(mydat_wide$month)
# The actual triangle computation
newcoord_up <- make_triangles(mydat_wide$month, mydat_wide$Species)
newcoord_down <- make_triangles(mydat_wide$month, mydat_wide$Species, point = "down")
# just a dirty trick for renaming
newcoord_down <- newcoord_down %>% select(xdown = x, ydown = y)
# you need to repeat each row of your previous data frame 3 times
repdata <- map_df(1:nrow(mydat_wide), function(i) mydat_wide[rep(i, 3), ])
newdata <- bind_cols(repdata, newcoord_up, newcoord_down)

p1 <- ggplot(newdata) +
  geom_polygon(aes(x = x, y = y, fill = Act, group = interaction(Species, month)), color = "white") +
  scale_fill_gradient2(low = "#4d9221", mid="white", high = "#c51b7d", na.value="grey86", limits = c(-3, 3), oob=scales::squish) +
  geom_polygon(aes(x = xdown, y = ydown, fill = Exp, group = interaction(Species, month)), color = "white") +
  scale_fill_gradient2(low = "#4d9221", mid="white", high = "#c51b7d", na.value="grey86", limits = c(-3, 3), oob=scales::squish) +
  scale_x_continuous(breaks = seq_along(levels(mydat_wide$month)), 
                     labels = unique(levels(mydat_wide$month)),
                     expand = c(0, 0), limits = c(0.5, NA))+
  scale_y_continuous(breaks = seq_along(levels(mydat_wide$Species)),
                     labels = unique(levels(mydat_wide$Species)),
                     expand = c(0, 0), limits = c(0.5, NA))+
  #guides(color=guide_legend("Act"), fill = "none")+ 
  #guides(fill=guide_legend(title='Risk z'))+
  ggtitle(paste0(length(olp), " SPs"))+
  theme(
	  axis.title = element_blank(),
	  axis.text = element_text(color="black"),
	  axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),
	  legend.position="bottom"
	)

ggsave(paste0(predictionPath,"prediction_together_heatmap_l.pdf"), p1, width = 16, height = 42, units = "cm")




################ 
# short heatmap
################


mydat <- smy_s
#mydat[mydat[,4]=="Exp",3] <- NA
#mydat[,3] <- NA
#mydat[,4] <- "Exp"
#mydat <- rbind(smy_s,mydat)

colnames(mydat) <- c("Species","month","measurement","flower_att")
mydat[,2] <- as.character(mydat[,2])


make_triangles <- function(x, y, point = "up") {
  x <- as.integer(as.factor((x)))
  y <- as.integer(as.factor((y)))

  if (point == "up") {
    newx <- sapply(x, function(x) {
      c(x - 0.5, x - 0.5, x + 0.5)
    }, simplify = FALSE)
    newy <- sapply(y, function(y) {
      c(y - 0.5, y + 0.5, y + 0.5)
    }, simplify = FALSE)
  } else if (point == "down") {
    newx <- sapply(x, function(x) {
      c(x - 0.5, x + 0.5, x + 0.5)
    }, simplify = FALSE)
    newy <- sapply(y, function(y) {
      c(y - 0.5, y - 0.5, y + 0.5)
    }, simplify = FALSE)
  }
  data.frame(x = unlist(newx), y = unlist(newy))
}

# required, otherwise you cannot use the values as fill
mydat_wide <- mydat %>% pivot_wider(names_from = "flower_att", values_from = "measurement")
# making your ordered months factor
mydat_wide$month <- factor(mydat_wide$month)
# The actual triangle computation
newcoord_up <- make_triangles(mydat_wide$month, mydat_wide$Species)
newcoord_down <- make_triangles(mydat_wide$month, mydat_wide$Species, point = "down")
# just a dirty trick for renaming
newcoord_down <- newcoord_down %>% select(xdown = x, ydown = y)
# you need to repeat each row of your previous data frame 3 times
repdata <- map_df(1:nrow(mydat_wide), function(i) mydat_wide[rep(i, 3), ])
newdata <- bind_cols(repdata, newcoord_up, newcoord_down)

p2 <- ggplot(newdata) +
geom_polygon(aes(x = x, y = y, fill = Act, group = interaction(Species, month)), color = "white") +
  scale_fill_gradient2(low = "#4d9221", mid="white", high = "#c51b7d", na.value="grey86", limits = c(-3, 3), oob=scales::squish) +
  geom_polygon(aes(x = xdown, y = ydown, fill = Exp, group = interaction(Species, month)), color = "white") +
  scale_fill_gradient2(low = "#4d9221", mid="white", high = "#c51b7d", na.value="grey86", limits = c(-3, 3), oob=scales::squish) +
  scale_x_continuous(breaks = seq_along(levels(mydat_wide$month)), 
                     labels = unique(levels(mydat_wide$month)),
                     expand = c(0, 0), limits = c(0.5, NA))+
  scale_y_continuous(breaks = seq_along(levels(mydat_wide$Species)),
                     labels = unique(levels(mydat_wide$Species)),
                     expand = c(0, 0), limits = c(0.5, NA))+
  #guides(color=guide_legend("Act"), fill = "none")+ 
  #guides(fill=guide_legend(title='Risk z'))+
  ggtitle(paste0(length(olp), " SPs"))+
  theme(
	  axis.title = element_blank(),
	  axis.text = element_text(color="black"),
	  axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),
	  legend.position="bottom"
	)

ggsave(paste0(predictionPath,"prediction_together_heatmap_s.pdf"), p2, width = 4.4, height = 42, units = "cm")






smy <- read.csv(paste0(predictionPath,"prediction_act.csv"),row.names=1,header=T,check.names=F)

smy_l <- smy[, setdiff(colnames(smy),notGenomeWide_alt) ]
smy_s <- smy[, notGenomeWide_alt]

smy_l.median <- apply(smy_l,1,function(x) median(x,na.rm=T))
smy_s.median <- apply(smy_s,1,function(x) median(x,na.rm=T))
	
fg.df <- data.frame(
	gene = olp,
	value = smy_l.median[olp]
)

fg.df[,1] <- factor(fg.df[,1],levels=rev(olp))

library(ggplot2)
p11 <- ggplot(fg.df,aes(x = gene, y = value))+
	geom_bar(stat="identity", width=0.01, color="grey66")+
	geom_hline(yintercept=0, color = "darkgrey", linewidth=0.8)+
	geom_point(color="#aa8800")+
	ylab("Median risk z\n(Activity)")+
	theme_classic()+ 
	theme(
	  panel.grid = element_blank(),
	  panel.background = element_rect(fill = "grey95"),
	  axis.text = element_text(colour = "black"),
	  axis.title = element_text(colour = "black"),
	  axis.ticks.y = element_blank(),
	  axis.title.y = element_blank(),
	  axis.line.y.left = element_line(color = 'white'),
	  axis.line.y.right = element_line(color = 'black')
	)+ coord_flip()

ggsave(paste0(predictionPath,"prediction_together_heatmap_rank_l.pdf"), p11, width = 4.5, height =32, units = "cm",)


fg.df <- data.frame(
	gene = olp,
	value = smy_s.median[olp]
)

fg.df[,1] <- factor(fg.df[,1],levels=rev(olp))

library(ggplot2)
p22 <- ggplot(fg.df,aes(x = gene, y = value))+
	geom_bar(stat="identity", width=0.01, color="grey66")+
	geom_hline(yintercept=0, color = "darkgrey", linewidth=0.8)+
	geom_point(color="#aa8800")+
	ylab("Median risk z\n(Activity)")+
	theme_classic()+ 
	theme(
	  panel.grid = element_blank(),
	  panel.background = element_rect(fill = "grey95"),
	  axis.text = element_text(colour = "black"),
	  axis.title = element_text(colour = "black"),
	  axis.ticks.y = element_blank(),
	  axis.title.y = element_blank(),
	  axis.line.y.left = element_line(color = 'white'),
	  axis.line.y.right = element_line(color = 'black')
	)+ coord_flip()

ggsave(paste0(predictionPath,"prediction_together_heatmap_rank_s.pdf"), p22, width = 4.5, height =32, units = "cm",)





################ 
# annotation
################


lqlist <- read.csv(paste0(predictionPath,"list_annotation.csv"),row.names=1,header=T)

lqlist[lqlist[,1]=="Anti-tumor (verified in-house)",1] <- "Anti-tumor"
lqlist[lqlist[,1]=="Pro-tumor (verified in-house)",1] <- "Pro-tumor"


mylist <- cbind(olp,Category="Unknown")
rownames(mylist) <- mylist[,1]
mylist <- mylist[,-1,drop=F]

olp_olp <- intersect(rownames(mylist),rownames(lqlist))

mylist[olp_olp,1] <- lqlist[olp_olp,1]

fg.df <- data.frame(mylist,x="x",y=rownames(mylist))
fg.df[,3] <- factor(fg.df[,3],levels=rev(fg.df[,3]))
#fg.df[,1] <- factor(fg.df[,1],levels=c("Anti-tumor","Pro-tumor","Anti-tumor (verified in-house)","Pro-tumor (verified in-house)","Conflict/No effect","Unknown"))
#fg.df[,1] <- factor(fg.df[,1],levels=c("Anti-tumor","Pro-tumor","Conflict/No effect","Unknown"))
fg.df[,1] <- factor(fg.df[,1],levels=c("Anti-tumor","Pro-tumor","Conflict","No effect","Unknown"))

library(ggplot2)
p3 <- ggplot(fg.df,aes(x = x, y = y, fill= Category))+
	geom_tile(color="white",size = 0.5)+
	scale_fill_manual(values=c("#4d9221","#c51b7d","darkgrey","darkgrey","black"))+
	theme(
	  panel.background = element_blank(),
      panel.grid = element_blank(),
	  axis.title = element_blank(),
	  axis.text.x = element_blank(),
	  axis.text.y = element_text(color = ifelse(rev(olp)%in%olp_neg, "#4d9221", "#c51b7d")),
	  axis.ticks.x = element_blank(),
	  legend.position = "bottom"
	)
ggsave(paste0(predictionPath,"prediction_together_anno.pdf"), p3, width = 8.6, height = 40, units = "cm")



olp_annotation <- read.csv(paste0(predictionPath,"prediction_olp_genomeWide_test_signif_annotation.csv"),row.names=1,as.is=T)
olp_annotation <- cbind(olp_annotation, followup="")
for(gene in rownames(olp_annotation))
{
	if(olp_annotation[gene,1]==lqlist[gene,1])
	{
		next
	}else{
		olp_annotation[gene,"followup"] <- lqlist[gene,1]
	}
}
olp_annotation <- olp_annotation[,c(1,3,2)]



library(patchwork)
p_comb <- p3 | p1+theme(axis.text.y = element_blank()) | 
	p11+theme(axis.text.y = element_blank()) |
	p2+theme(axis.text.y = element_blank()) |
	p22+theme(axis.text.y = element_blank())
	
p_comb <- p_comb + plot_layout(widths = c(0.25, 10, 1, 0.6, 1))

ggsave(paste0(predictionPath,"prediction_together_heatmap_comb.pdf"), p_comb, width = 27, height =45, units = "cm")
ggsave(paste0(predictionPath,"prediction_together_heatmap_comb.png"), p_comb, width = 27, height =45, units = "cm")






smy_test1 <- read.csv(paste0(predictionPath,"prediction_act_genomeWide_test.csv"),row.names=1,header=T)
smy_test2 <- read.csv(paste0(predictionPath,"prediction_exp_genomeWide_test.csv"),row.names=1,header=T)

olp <- intersect(rownames(smy_test1),rownames(smy_test2))

fg.df <- data.frame(Act=smy_test1[olp,1],Exp=smy_test2[olp,1])
rownames(fg.df) <- olp

maxmax <- max(fg.df)
minmin <- min(fg.df)

library(ggplot2)
p1 <- ggplot(fg.df, aes(x = Exp, y = Act, label = rownames(fg.df)))+
	geom_hline(yintercept=0, color = "darkgrey", linewidth=0.6)+
	geom_vline(xintercept=0, color = "darkgrey", linewidth=0.6)+
	xlim(minmin,maxmax)+
	ylim(minmin,maxmax)+
	geom_text()+
	coord_equal()+
	theme_bw()+ 
	theme(
	  panel.grid = element_blank(),
	  panel.background = element_blank(),
	  axis.text = element_text(size=10,colour = "black"),
	  axis.title = element_text(size=10,colour = "black"),
	  legend.position="none"
	)

ggsave(paste0(predictionPath,"prediction_together_scatter.png"), p1, width = 50, height = 50, dpi=200, units = "cm",limitsize = FALSE)




smy_test1 <- read.csv(paste0(predictionPath,"prediction_act_genomeWide_test.csv"),row.names=1,header=T)
smy_test2 <- read.csv(paste0(predictionPath,"prediction_exp_genomeWide_test.csv"),row.names=1,header=T)

smy_test1 <- smy_test1[rev(rownames(smy_test1)),]

genes1 <- rownames(smy_test1)[smy_test1[,"is.signif"]==1]
genes2 <- rownames(smy_test2)[smy_test2[,"is.signif"]==1]

fg.df <- data.frame(Act=smy_test1[genes1,1],Exp=smy_test2[genes1,1])
rownames(fg.df) <- genes1

fg.df <- cbind(fg.df, group=NA)

neg_max <- max(fg.df[fg.df[,1]<0,1])
pos_min <- min(fg.df[fg.df[,1]>0,1])

fg.df[,"group"] <-  fg.df[,1] > pos_min & fg.df[,2] > pos_min | fg.df[,1]< neg_max & fg.df[,2]< neg_max


library(ggplot2)
p1 <- ggplot(fg.df, aes(x = Exp, y = Act, label = rownames(fg.df)))+
	geom_hline(yintercept=0, color = "darkgrey", linewidth=0.6)+
	geom_vline(xintercept=0, color = "darkgrey", linewidth=0.6)+
	geom_hline(yintercept=pos_min, color = "grey", linewidth=0.6, linetype="dashed")+
	geom_hline(yintercept=neg_max, color = "grey", linewidth=0.6, linetype="dashed")+
	geom_vline(xintercept=pos_min, color = "grey", linewidth=0.6, linetype="dashed")+
	geom_vline(xintercept=neg_max, color = "grey", linewidth=0.6, linetype="dashed")+
	geom_text(aes(colour=group))+
	coord_equal()+
	theme_bw()+ 
	theme(
	  panel.grid = element_blank(),
	  panel.background = element_blank(),
	  axis.text = element_text(size=10,colour = "black"),
	  axis.title = element_text(size=10,colour = "black"),
	  legend.position="none"
	)

ggsave(paste0(predictionPath,"prediction_together_scatter2.png"), p1, width = 50, height = 50, dpi=200, units = "cm",limitsize = FALSE)






smy_test1 <- read.csv(paste0(predictionPath,"prediction_act_genomeWide_test.csv"),row.names=1,header=T)
smy_test2 <- read.csv(paste0(predictionPath,"prediction_exp_genomeWide_test.csv"),row.names=1,header=T)

smy_test1 <- smy_test1[rev(rownames(smy_test1)),]

genes1 <- rownames(smy_test1)[smy_test1[,"is.signif"]==1]
genes1 <- rownames(smy_test2)[smy_test2[,"is.signif"]==1]

fg.df <- data.frame(Act=smy_test1[genes1,1],Exp=smy_test2[genes1,1])
rownames(fg.df) <- genes1

fg.df <- cbind(fg.df, group=NA)

neg_max <- max(fg.df[fg.df[,1]<0,1])
pos_min <- min(fg.df[fg.df[,1]>0,1])

fg.df[,"group"] <-  fg.df[,1] > pos_min & fg.df[,2] > pos_min | fg.df[,1]< neg_max & fg.df[,2]< neg_max


library(ggplot2)
p1 <- ggplot(fg.df, aes(x = Exp, y = Act, label = rownames(fg.df)))+
	geom_hline(yintercept=0, color = "darkgrey", linewidth=0.6)+
	geom_vline(xintercept=0, color = "darkgrey", linewidth=0.6)+
	geom_hline(yintercept=pos_min, color = "grey", linewidth=0.6, linetype="dashed")+
	geom_hline(yintercept=neg_max, color = "grey", linewidth=0.6, linetype="dashed")+
	geom_vline(xintercept=pos_min, color = "grey", linewidth=0.6, linetype="dashed")+
	geom_vline(xintercept=neg_max, color = "grey", linewidth=0.6, linetype="dashed")+
	geom_text(aes(colour=group))+
	coord_equal()+
	theme_bw()+ 
	theme(
	  panel.grid = element_blank(),
	  panel.background = element_blank(),
	  axis.text = element_text(size=10,colour = "black"),
	  axis.title = element_text(size=10,colour = "black"),
	  legend.position="none"
	)

ggsave(paste0(predictionPath,"prediction_together_scatter3.png"), p1, width = 50, height = 50, dpi=200, units = "cm",limitsize = FALSE)








smy_test1 <- read.csv(paste0(predictionPath,"prediction_act_genomeWide_test.csv"),row.names=1,header=T)

genes_act_neg <- rownames(smy_test1)[smy_test1[,"is.signif"]==1&smy_test1[,"smy.median"]<0]
genes_act_pos <- rev(rownames(smy_test1)[smy_test1[,"is.signif"]==1&smy_test1[,"smy.median"]>0])


smy <- read.csv(paste0(predictionPath,"prediction_act.csv"),as.is=T,row.names=1,header=T)


smy_l <- smy[, setdiff(colnames(smy),notGenomeWide_alt) ]
smy_s <- smy[, notGenomeWide_alt]
	
smy_l.median <- apply(smy_l,1,function(x) median(x,na.rm=T))
smy_s.median <- apply(smy_s,1,function(x) median(x,na.rm=T))
	

fg.df <- data.frame(
	x=smy_l.median,
	y=smy_s.median,
	group="Rest"
)

fg.df[genes_act_neg,"group"] <- "Anti"
fg.df[genes_act_pos,"group"] <- "Pro"

fg.df[,"group"] <- factor(fg.df[,"group"], levels=c("Anti","Rest","Pro"))

cor_res <- cor.test(fg.df[,1],fg.df[,2])
rv <- round(cor_res$estimate,2)
pv <- signif(cor_res$p.value,2)


library(ggplot2)
p1 <- ggplot(fg.df,aes(x=x, y=y, color=group)) + 
	geom_point( alpha=0.6, size=0.6)+
	scale_color_manual( values=c("#4d9221","#bababa","#c51b7d") )+
	annotate("text", x = 1, y=-1.7, label = paste0("r = ",rv))+
	guides(color = guide_legend(override.aes = list(size=1.6)))+
	xlab("Median z (Genome-wide)")+
	ylab("Median z (Targeted)")+
	theme_classic()+ 
	theme(
		plot.background = element_blank(),
		panel.grid = element_blank(),
		legend.position = "bottom",
		legend.title = element_blank(),
		legend.key.size = unit(0.7, 'lines')
	)			
ggsave(paste0(predictionPath,"prediction_genome_targeted_compare.pdf"), p1, width = 6, height = 7, dpi=400, units = "cm")







n_downsampling <- c(4,5,6,8,10,12,14,16,18,40,50,60,70,80,90,100,200,400,600,800,1000,2000,4000,6000,8000,10000,12000,14000,16000,18000)
n_downsampling <- c(5,10,15,20,40,60,80,100,200,400,600,800,1000,2000,4000,6000,8000)

smy <- read.csv(paste0(predictionPath,"prediction_act.csv"),as.is=T,row.names=1,header=T)
smy_l <- smy[, setdiff(colnames(smy),notGenomeWide_alt) ]
smy_l.median <- apply(smy_l,1,function(x) median(x,na.rm=T))

stat <- data.frame()
for(n in n_downsampling)
{
	for(i in 1:10)
	{
		smy_s <- read.csv(paste0(predictionPath,"immunotherapy/",n,"/prediction_act_",i,".csv"),as.is=T,row.names=1,header=T)
		smy_s.median <- apply(smy_s,1,function(x) median(x,na.rm=T))
		
		fg.df <- data.frame(
			x=smy_l.median,
			y=smy_s.median
		)
		
		cor_res <- cor.test(fg.df[,1],fg.df[,2])
		rv <- round(cor_res$estimate,2)
		pv <- signif(cor_res$p.value,2)
		
		stat[i,as.character(n)] <- rv
		
		if(i==2&n%in%c(5,10,20,100,1000))
		{
			library(ggplot2)
			p1 <- ggplot(fg.df,aes(x=x, y=y)) + 
				geom_point(alpha=0.3, size=0.6, color="skyblue")+
				annotate("text", x = 1, y=-0.7, label = paste0("r = ",rv))+
				xlab("Median z (Genome-wide)")+
				ylab("Median z (Downsampling)")+
				ggtitle(paste0("Gene coverage: ",n))+
				theme_classic()+ 
				theme(
					plot.background = element_blank(),
					panel.grid = element_blank(),
					plot.title = element_text(size = 10),
					legend.position = "none",
					legend.title = element_blank(),
					legend.key.size = unit(0.7, 'lines')
				)			
			ggsave(paste0(predictionPath,"prediction_genome_downsample_compare_",n,".png"), p1, width = 6.2, height = 6.2, dpi=500, units = "cm")
		}
	}
}


fg.df <- reshape2::melt(as.matrix(stat))
fg.df[,2] <- as.character(fg.df[,2])
fg.df[,2] <- factor(fg.df[,2], levels=as.character(n_downsampling) )

library(ggplot2)
p1 <- ggplot(fg.df,aes(x=Var2,y=value)) + 
	geom_hline(yintercept=0, color = "grey", linewidth=0.8, linetype="dashed")+
	geom_boxplot( aes(color=Var2), alpha=0, width=0.5, outlier.shape = NA)+
	geom_jitter(color="darkgrey",alpha=0.3, size=1, width=0.1)+
	scale_y_continuous(breaks = seq(0, 1, by = 0.2))+
	ylab("Correlation r")+
	xlab("Gene coverage")+
	theme_classic()+ 
	theme(
		panel.background = element_blank(),
		panel.grid = element_blank(),
		axis.title = element_text(colour = "black"),
		axis.text = element_text(colour = "black"),
		axis.text.x = element_text(size=10),
		legend.position="none"
	)
	
ggsave(paste0(predictionPath,"prediction_genome_downsample_compare_summary.png"), p1, width = 26, height = 6, dpi=500, units = "cm")









notGenomeWide <- c(
"HNSCC_ICB_Foy2022",
"NSCLC_PD1_Foy2022",
"PanCancer_PD1_Prat2017"
)
notGenomeWide_alt <- c(
"HNSCC_ICB_Foy2022.OS",
"NSCLC_PD1_Foy2022.OS",
"PanCancer_PD1_Prat2017.PFS"
)

dataPath <- "/data/rub2/data/Immunotherapy/"

smy <- data.frame()
for(cancer in notGenomeWide)
{
	expr <- as.matrix(read.csv(gzfile(paste0(dataPath,cancer,".expression.gz")),as.is=T,sep="\t",check.names=F))
	expr <- expr[apply(expr,1,function(x) sum(x>0)>5),]	
	
	rownames(expr) <- transferSymbol(rownames(expr))
	expr <- rm_duplicates(expr)
	
	print("APLN"%in%rownames(expr))
	
	smy[cancer,"Secreted"] <- sum(rownames(expr)%in%SPs)
	smy[cancer,"Other"] <- dim(expr)[1]-sum(rownames(expr)%in%SPs)
}

fg.df <- reshape2::melt(as.matrix(smy))

library(ggplot2)
p1 <- ggplot(fg.df, aes(x=Var1, y=value, fill=Var2)) +
	geom_bar(stat="identity", color="white", width=0.9, alpha=0.8) +
	scale_fill_manual( values=c("#91bfdb","#fc8d59") )+
	xlab(" ")+
	ylab("# Gene")+
	theme_classic()+ 
	theme(
		panel.background = element_blank(),
		panel.grid = element_blank(),
		axis.title = element_text(colour = "black"),
		axis.text = element_text(colour = "black"),
		axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),
		#legend.position = c(0.75, 0.8),
		legend.position = "right",
		legend.title=element_blank()
	)
ggsave(paste0(predictionPath,"Targeted_gene_count.pdf"), p1, width = 7, height = 10, dpi=500, units = "cm", limitsize =FALSE)



smy <- read.csv(paste0(predictionPath,"prediction_act.csv"),as.is=T,row.names=1,header=T)
smy <- smy["PTK7",notGenomeWide_alt]
colnames(smy) <- notGenomeWide

fg.df <- reshape2::melt(as.matrix(smy))
pv <- wilcox.test(fg.df[,"value"], mu=0)$p.value

library(ggplot2)
p1 <- ggplot(fg.df,aes(x = Var2, y = value))+
	geom_bar(stat="identity", width=0.01, color="grey66")+
	geom_point(color="#91bfdb",size=3)+
	#annotate("text", x = 3, y=1, label = paste0("p = ",pv))+
	xlab(" ")+
	ylab("PTK7\nRisk score z")+
	theme_classic()+ 
	theme(
		panel.background = element_blank(),
		panel.grid = element_blank(),
		axis.title = element_text(colour = "black"),
		axis.text = element_text(colour = "black"),
		axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),
		legend.title=element_blank()
	)

ggsave(paste0(predictionPath,"Targeted_gene_PTK7.png"), p1, width = 6, height = 10, dpi=500, units = "cm", limitsize =FALSE)













######
# risk score water fall plot
######

gene <- "LY86"
smy <- read.csv(paste0(predictionPath,"prediction_act.csv"),as.is=T,row.names=1,header=T)
smy <- as.matrix(smy[, setdiff(colnames(smy),notGenomeWide_alt) ])

smy_test1 <- read.csv(paste0(predictionPath,"prediction_act_genomeWide_test.csv"),row.names=1,header=T)
pv <- signif(smy_test1[gene,"p.value"],2)

temp <- unlist(smy[gene,])
temp_sorted <- sort(temp)

fg.df <- data.frame(
	dataset = names(temp_sorted),
	gene = gene,
	value = temp_sorted
)
fg.df[,1] <- factor(fg.df[,1],levels=names(temp_sorted))

library(ggplot2)
p1 <- ggplot(fg.df,aes(x = dataset, y = value))+
	#geom_bar(stat="identity", width=0.01, color="grey66")+
	geom_hline(yintercept=0, color = "grey", linewidth=0.8, linetype="dashed")+
	geom_point(color="#619CFF",size=0.7)+
	annotate("text", x = 12, y=1.0, label = paste0(gene))+
	annotate("text", x = 12, y=0.6, label = paste0("p = ",pv))+
	xlab("Cohort")+
	ylab(paste0("Risk score z"))+
	theme_classic()+ 
	theme(
	  panel.grid = element_blank(),
	  panel.background = element_blank(),
	  axis.text = element_text(colour = "black"),
	  #axis.text.x = element_text(size=8, angle = 90, hjust = 1, vjust = 0.5),
	  axis.text.x = element_blank(),
	  axis.ticks.x = element_blank(),
	  axis.title = element_text(colour = "black")
	)

ggsave(paste0(predictionPath,"prediction_",gene,".png"), p1, width = 5.6, height =5.5, dpi=400, units = "cm",)
ggsave(paste0(predictionPath,"prediction_",gene,".pdf"), p1, width = 6.8, height =6.8, dpi=400, units = "cm",)



######
# risk score box plot
######

library(ggplot2)
p1 <- ggplot(fg.df,aes(x=gene,y=value)) + 
	geom_hline(yintercept=0, color = "black", linewidth=0.8, linetype="dashed")+
	geom_jitter(color="darkgrey",alpha=0.5, size=1, width=0.2)+
	geom_boxplot(color="purple", alpha=0, width=0.5, outlier.shape = NA)+
	scale_color_manual(values=c("#FF2A2A","#008000"))+
	annotate("text", x = 1, y=2, label = paste0("p = "))+
	annotate("text", x = 1, y=1.5, label = pv)+
	ylab("Risk score z")+
	xlab("")+
	theme_classic()+ 
	theme(
		panel.background = element_blank(),
		panel.grid = element_blank(),
		axis.title = element_text(colour = "black"),
		axis.text = element_text(colour = "black"),
		axis.title.x = element_blank(),
		axis.ticks = element_blank(),
		legend.position="none"
	)
	
ggsave(paste0(predictionPath,"prediction_",gene,"_box.pdf"), p1, width = 3.3, height = 5, dpi=400, units = "cm")





genes <- c("GZMA","IL11")

smy <- read.csv(paste0(predictionPath,"prediction_act.csv"),as.is=T,row.names=1,header=T)
smy <- as.matrix(smy[, setdiff(colnames(smy),notGenomeWide_alt) ])

fg.df <- reshape::melt(smy[genes,])
fg.df <- cbind(fg.df, group=fg.df[,1]%in%genes[1:6])

smy_test1 <- read.csv(paste0(predictionPath,"prediction_act_genomeWide_test.csv"),row.names=1,header=T)
pvs <- signif(smy_test1[rownames(smy_test1)%in%genes,"p.value"],2)


library(ggplot2)
p1 <- ggplot(fg.df,aes(x=X1,y=value)) + 
	geom_hline(yintercept=0, color = "grey", linewidth=0.6, linetype="dashed")+
	geom_jitter(color="darkgrey",alpha=0.3, size=1, width=0.1)+
	geom_boxplot(color="black", alpha=0, width=0.5, outlier.shape = NA)+
	scale_color_manual(values=c("#FF2A2A","#008000"))+
	annotate("text", x = 1, y=4, label = paste0("p = ", pvs[1]))+
	annotate("text", x = 2, y=-2.0, label = paste0("p = ", pvs[2]))+
	ylab("Risk score z")+
	xlab("")+
	theme_classic()+ 
	theme(
		panel.background = element_blank(),
		panel.grid = element_blank(),
		axis.title = element_text(colour = "black"),
		axis.text = element_text(colour = "black"),
		#axis.ticks = element_line(colour = "darkgrey"),
		#axis.line = element_line(colour = "darkgrey"),
		axis.text.x = element_text(color = ifelse(genes%in%genes[1], "#4d9221","#c51b7d")),
		axis.title.x = element_blank(),
		legend.position="none"
	)
	
ggsave(paste0(predictionPath,"prediction_GZMA_IL11.pdf"), p1, width = 5.5, height = 4.8, dpi=500, units = "cm")


smy_t <- t(smy[c(genes,"LY86"),])
smy_t[order(smy_t[,3]),]




GZMs <- c("GZMA","GZMB","GZMH","GZMK","GZMM","IFNG","IL11","IGF2")

smy <- read.csv(paste0(predictionPath,"prediction_act.csv"),as.is=T,row.names=1,header=T)
smy <- as.matrix(smy[, setdiff(colnames(smy),notGenomeWide_alt) ])

fg.df <- reshape::melt(smy[GZMs,])
fg.df <- cbind(fg.df, group=fg.df[,1]%in%GZMs[1:6])

smy_test1 <- read.csv(paste0(predictionPath,"prediction_act_genomeWide_test.csv"),row.names=1,header=T)
pvs <- signif(smy_test1[rownames(smy_test1)%in%GZMs,"p.value"],2)

library(ggplot2)
p1 <- ggplot(fg.df,aes(x=X1,y=value)) + 
	geom_hline(yintercept=0, color = "grey", linewidth=0.8, linetype="dashed")+
	geom_jitter(color="darkgrey",alpha=0.3, size=1, width=0.1)+
	geom_boxplot(color="black", alpha=0, width=0.5, outlier.shape = NA)+
	annotate("text", x = 0.8, y=4.5, label = paste0("p = "))+
	annotate("text", x = 1, y=3.2, label = paste0(pvs[1]))+
	annotate("text", x = 2, y=2.4, label = paste0(pvs[2]))+
	annotate("text", x = 3, y=3.2, label = paste0(pvs[3]))+
	annotate("text", x = 4, y=2.4, label = paste0(pvs[4]))+
	annotate("text", x = 5, y=3.2, label = paste0(pvs[5]))+
	annotate("text", x = 6, y=2.4, label = paste0(pvs[6]))+
	annotate("text", x = 7, y=-2, label = paste0(pvs[7]))+
	annotate("text", x = 8, y=-3, label = paste0(pvs[8]))+
	ylab("Risk score z")+
	xlab("")+
	theme_classic()+ 
	theme(
		panel.background = element_blank(),
		panel.grid = element_blank(),
		axis.title = element_text(colour = "black"),
		axis.text = element_text(colour = "black"),
		axis.text.x = element_text(color = ifelse(GZMs%in%GZMs[1:6], "#4d9221","#c51b7d")),
		axis.title.x = element_blank(),
		legend.position="none"
	)
	
ggsave(paste0(predictionPath,"prediction_GZM.png"), p1, width = 12, height = 5.5, dpi=400, units = "cm")






######
# survival plot
######


gene <- "AOAH"
datasets <- c("Melanoma_PD1_Gide2019.OS")

gene <- "GZMA"
datasets <- c("mRCC_Atezo+Bev_McDermott2018.PFS")

gene <- "IGF2"
datasets <- c("Hepatocellular_Atezo+Bev_Finn2020.OS")

gene <- "IFNG"
datasets <- c("Hepatocellular_Atezo+Bev_Finn2020.OS")

gene <- "IL11"
datasets <- c("Hepatocellular_Atezo+Bev_Finn2020.OS")

datasets <- c("RCC_Avelumab+Axitinib_Motzer2020.PFS")

gene <- "LY86"
datasets <- c(
	"Melanoma_Ipilimumab_VanAllen2015.OS",
	#"Pancreatic_Nivolumab_Padron2022.OS_Nivo+Sotiga+Chemo",
	#"Melanoma_PD1_Liu2019.OS_Prog",
	#"mRCC_Atezo+Bev_McDermott2018.PFS"
	#"PanCancer_ICB_Li2023.OS",
	#"NSCLC_PD1orPDL1_Jung2019.PFS",
	#"SCLC_Durvalumab_Roper2021.OS"
	"NSCLC_Atezolizumab_Patil2022-OAK.OS",
	"Urothelial_Atezo+Chemo_Hamidi2024-IMvigor130.OS"
	)

source(file.path("./survival_util.R"))
dataPath <- "/data/rub2/data/Immunotherapy/"

for(dataset in datasets)
{
	study <- strsplit(dataset,".",fixed=T) [[1]][1]
	clini <- strsplit(dataset,".",fixed=T) [[1]][2]
	survival <- read.csv(gzfile(paste0(dataPath,dataset)),as.is=T,sep="\t",check.names=F)
	colnames(survival)[1] <- "OS"
	colnames(survival)[2] <- "OS.Event"
	
	if(dataset%in%c("mRCC_Atezo+Bev_McDermott2018.PFS","RCC_Avelumab+Axitinib_Motzer2020.PFS"))
	{
		xtext = "Progression-Free (Months)"
	}else if(dataset%in%c("NSCLC_PD1orPDL1_Jung2019.PFS")){
		xtext = "Progression-Free (Days)"
	}else if(dataset=="NSCLC_Atezolizumab_Patil2022-OAK.OS"){
		xtext = "Overall (Months)"
	}else{
		xtext = "Overall (Days)"
	}
	
	load(paste0(predictionPath,"immunotherapy/",study,".RData"))
	Act <- as.matrix(res$zscore)
	Act <- expand_rows(Act)
	Act <- t(Act)
	
	
	data <- Act [,c("AOAH","VEGFA","GZMA","IFNG","IGF2","IL11","LY86")]
	
	
	out <- Beibei_revised(data, survival, margin=5)
	data <- out[[1]]
	survival <- out[[2]]
	result <- out[[3]]
	cutoff <- round(result[gene,"thres.opt"],3)
	
	
	library(survival)
	olp <- intersect(rownames(survival),rownames(data))
	X_olp <- survival[olp,,drop=F]
	Y_olp <- data[olp,,drop=F]
	
	comb <- cbind(X_olp, Act=Y_olp[,gene])

	coxmodel_fit <- coxph(Surv(OS, OS.Event) ~ ., data = comb)
	coxmodel_obj <- summary(coxmodel_fit)
	zs <- coxmodel_obj$coefficients["Act","z"]
	pv <- coxmodel_obj$coefficients["Act","Pr(>|z|)"]
	
	
	surv_o <- Surv(survival[,1],survival[,2])
	groups <- as.character(data[,gene]>cutoff)
	surv_c  <- survfit(surv_o ~ groups)
	
	
	library("survminer")
	surv.df <- cbind(survival,groups)
	
	fit <- survfit(Surv(OS, OS.Event) ~ groups, data = surv.df)
	
	p2 <- ggsurvplot(fit, data=surv.df, palette = c("#0000FF","#FF0000"), legend.labs = c(paste0("Low (n=",sum(groups=="FALSE"),")"),paste0("High (n=",sum(groups=="TRUE"),")")) )$plot+
		xlab(xtext)+
		ylab("Fraction")+
		annotate("text", x = max(survival[,1])*0.7, y=0.88, label = paste0("z = ", round(zs,2), "\np = ", signif(pv,2)), size=4 )

	ggsave(paste0(predictionPath,gene,"_",dataset,".pdf"), p2, width = 3.1, height = 3.5)
	
}




gene <- "LY86"
dataset <- "Melanoma_MAGEA3_Montoya2013.response"

dataPath <- "/data/rub2/data/Immunotherapy/"

study <- strsplit(dataset,".",fixed=T) [[1]][1]
clini <- strsplit(dataset,".",fixed=T) [[1]][2]
survival <- read.csv(gzfile(paste0(dataPath,dataset)),as.is=T,sep="\t",check.names=F)

#rownames(survival) <- survival[,1]
#survival <- survival[,-1,drop=F]
		
load(paste0(predictionPath,"immunotherapy/",study,".RData"))
data <- t(as.matrix(res$zscore)) [,c("VEGFA","GZMA","IFNG","WNT7B","LY86")]
	
olp <- intersect(rownames(survival),rownames(data))
X_olp <- survival[olp,,drop=F]
Y_olp <- data[olp,,drop=F]
	
fg.df <- data.frame(x=as.character(X_olp[,1]), y=Y_olp[,1])
fg.df[,1] <- ifelse(fg.df[,1]%in%"1","Responder","Non-Res.")

x1 <- fg.df[fg.df[,1]=="Responder",2]
x2 <- fg.df[fg.df[,1]=="Non-Res.",2]
pv <- signif(wilcox.test(x1, x2)$p.value,2)


library(ggplot2)
p1 <- ggplot(fg.df,aes(x=x,y=y)) + 
	geom_hline(yintercept=0, color = "grey", linewidth=0.8, linetype="dashed")+
	geom_boxplot( aes(color=x), alpha=0, width=0.5, outlier.shape = NA)+
	geom_jitter(color="darkgrey",alpha=0.3, size=1, width=0.1)+
	ylab("Activity (LY86)")+
	xlab("")+
	theme_classic()+ 
	theme(
		panel.background = element_blank(),
		panel.grid = element_blank(),
		axis.title = element_text(colour = "black"),
		axis.text = element_text(colour = "black"),
		axis.title.x = element_blank(),
		legend.position="none"
	)
ggsave(paste0(predictionPath, gene, "_", dataset,".png"), p1, width = 6, height = 5.5, dpi=400, units = "cm")





######
# GSEA
######

gene <- "LY86"

smy <- read.csv(paste0(predictionPath,"prediction_act.csv"),as.is=T,row.names=1,header=T,check.names=F)
smy <- as.matrix(smy[, setdiff(colnames(smy),notGenomeWide_alt) ])
datasets <- colnames(smy)[!is.na(smy[gene,])&smy[gene,]< -1.8]

z_sampleNo <- data.frame(smy[gene,])
z_sampleNo[order(z_sampleNo[,1]),]

for(dataset in datasets)
{
	study <- strsplit(dataset,".",fixed=T) [[1]][1]
	load(paste0(predictionPath,"immunotherapy/",study,".RData"))
	
	Act <- as.matrix(res$zscore)
	Act <- expand_rows(Act)
	
	z_sampleNo[dataset,"No"] <- ncol(Act)
	
	q25 <- quantile(Act[gene,])[2]
	q75 <- quantile(Act[gene,])[4]
	
	if(file.exists(paste0(dataPath,study,".","TPM",".gz")))
	{
		expr <- as.matrix(read.csv(gzfile(paste0(dataPath,study,".","TPM",".gz")),as.is=T,sep="\t",check.names=F))
	}else{
		expr <- as.matrix(read.csv(gzfile(paste0(dataPath,study,".","expression",".gz")),as.is=T,sep="\t",check.names=F))
	}
	
	expr <- expr[apply(expr,1,function(x) sum(x>0)>5),]
	
	rownames(expr) <- transferSymbol(rownames(expr))
	expr <- rm_duplicates(expr)
	
	
	library(limma)
	TT <- as.numeric(Act[gene,] > q75)
	WT <- as.numeric(Act[gene,] < q25)
	design <- cbind(TT,WT)
	fit <- lmFit(expr,design)
	cont.matrix <- makeContrasts(TTvsWT=TT-WT,levels=design)
	fit2 <- contrasts.fit(fit, cont.matrix)
	fit2 <- eBayes(fit2)
	res <- topTable(fit2,coef=1,number=nrow(expr))
	res <- res[order(res[,"t"],decreasing=T),]
	
	rnk <- res[,3]
	names(rnk) <- rownames(res)
	
	
	library(fgsea)
	data_MSigDB_path <- "/data/rub2/data/MSigDB/"
	gmtNames <- list.files(data_MSigDB_path)	
	
	for(gmtName in gmtNames)
	{
		gmt <- read.gmt(paste0(data_MSigDB_path,gmtName))
		gmt_length <- lapply(gmt, length)
		gmt <- gmt[gmt_length>=15&gmt_length<=500]
		
		fgseaRes <- fgsea(pathways = gmt, stats = rnk, nPermSimple = 100000)
		fgseaRes <- fgseaRes[order(fgseaRes[,5],decreasing=T),]
		
		write.csv(as.matrix(fgseaRes[,1:7]),paste0("/data/rub2/project/Secretome/results/LY86_RNAseq/immunotherapy/",dataset,"@",gene,"@",gmtName,".csv"),quote=F)
	}
		
}





fg.df <- data.frame()
p1 <- c(
	"HALLMARK_INFLAMMATORY_RESPONSE",
	"HALLMARK_INTERFERON_ALPHA_RESPONSE",
	"HALLMARK_INTERFERON_GAMMA_RESPONSE",
	"HALLMARK_MYC_TARGETS_V1","HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT")
p2 <- c(
	"GOBP_LEUKOCYTE_MEDIATED_IMMUNITY",
	"GOBP_ANTIGEN_RECEPTOR_MEDIATED_SIGNALING_PATHWAY",
	"GOBP_MITOTIC_SPINDLE_ORGANIZATION","GOBP_MITOCHONDRIAL_TRANSLATION")
p3 <- c(
	"HALLMARK_INFLAMMATORY_RESPONSE",
	"HALLMARK_INTERFERON_ALPHA_RESPONSE",
	"HALLMARK_INTERFERON_GAMMA_RESPONSE",
	"GOBP_LEUKOCYTE_MEDIATED_IMMUNITY",
	"GOBP_ANTIGEN_RECEPTOR_MEDIATED_SIGNALING_PATHWAY",
	"HALLMARK_MYC_TARGETS_V1",
	"HALLMARK_E2F_TARGETS",
	"HALLMARK_G2M_CHECKPOINT",
	"GOBP_MITOTIC_SPINDLE_ORGANIZATION",
	"GOBP_MITOCHONDRIAL_TRANSLATION"
)

for(dataset in datasets)
{
	fgseaRes1 <- read.csv(paste0("/data/rub2/project/Secretome/results/LY86_RNAseq/immunotherapy/",dataset,"@",gene,"@h.all.v2023.2.Hs.symbols.gmt.csv"),row.names=2)
	fgseaRes2 <- read.csv(paste0("/data/rub2/project/Secretome/results/LY86_RNAseq/immunotherapy/",dataset,"@",gene,"@c5.go.bp.v2023.2.Hs.symbols.gmt.csv"),row.names=2)
	
	#fg.df[rownames(fgseaRes1),dataset] <- fgseaRes1[rownames(fgseaRes1),"NES"]
	#fg.df[rownames(fgseaRes2),dataset] <- fgseaRes2[rownames(fgseaRes2),"NES"]
	
	fg.df[p1,dataset] <- fgseaRes1[p1,"NES"]
	fg.df[p2,dataset] <- fgseaRes2[p2,"NES"]
}


library(ComplexHeatmap)

mat <- as.matrix(fg.df[p3,])

#mat <- mat[,order(mat[1,],decreasing=T),drop=F]


#colnames(mat) <- sapply(strsplit(colnames(mat),".",fixed=T),function(x) return(x[1])) 
colnames(mat) <- rename[match(colnames(mat),rename[,1]),2]




rownames(mat) <- gsub("HALLMARK_","",rownames(mat)) 
rownames(mat) <- gsub("GOBP_","",rownames(mat)) 
rownames(mat) <- gsub("_PATHWAY","",rownames(mat)) 
rownames(mat) <- gsub("_"," ",rownames(mat)) 
rownames(mat) <- tolower(rownames(mat)) 

rownames(mat) <- stringr::str_to_title(rownames(mat)) 
rownames(mat) <- gsub("E2f","E2F",rownames(mat)) 
rownames(mat) <- gsub("Myc","MYC",rownames(mat)) 
rownames(mat) <- gsub("G2m","G2M",rownames(mat)) 




#png(paste0("/data/rub2/project/Secretome/results/LY86_RNAseq/immunotherapy/",gene,"_pathway.png"), width = 19, height = 17, res=500, units = "cm")
pdf(paste0("/data/rub2/project/Secretome/results/LY86_RNAseq/immunotherapy/",gene,"_pathway.pdf"), width = 7, height = 7)
	
ht <- Heatmap(mat,
	name = "NES",
	rect_gp = gpar(col = "white", lwd = 2),
	col = circlize::colorRamp2(c(-5, 0,5), c("#66bd63", "#ffffbf", "#f46d43")),
    row_names_side = "left",
    column_title_side = "bottom",
    row_title = " ", 
    column_title = "Cohort",
    column_names_max_height = max_text_width(
	        colnames(mat), 
	        gp = gpar(fontsize = 12)
	        ),
    show_heatmap_legend = TRUE,
	cluster_rows = FALSE,
	cluster_columns = FALSE
)

draw(ht)
dev.off()













for(gene in c("FCN1","APOL4","LY86","COL4A4"))
{

smy <- read.csv(paste0(predictionPath,"prediction_act.csv"),as.is=T,row.names=1,header=T)
smy <- as.matrix(smy[, setdiff(colnames(smy),notGenomeWide_alt) ])

smy_test1 <- read.csv(paste0(predictionPath,"prediction_act_genomeWide_test.csv"),row.names=1,header=T)
pv <- signif(smy_test1[gene,"p.value"],2)

if(gene=="FCN1")
{
	labelText <- paste0("p = ",pv)
}else{
	labelText <- paste0(pv)
}


temp <- unlist(smy[gene,])
temp_sorted <- sort(temp)

fg.df <- data.frame(
	gene = names(temp_sorted),
	value = temp_sorted
)
fg.df[,1] <- factor(fg.df[,1],levels=names(temp_sorted))
print(range(fg.df[,2]))
library(ggplot2)
p1 <- ggplot(fg.df,aes(x = gene, y = value))+
	#geom_bar(stat="identity", width=0.01, color="grey66")+
	geom_hline(yintercept=0, color = "grey", linewidth=0.8, linetype="dashed")+
	geom_point(color="#619CFF",size=0.5)+
	annotate("text", x = 20, y=1.3, label = labelText)+
	xlab(gene)+
	ylab(paste0("Risk score"))+
	ylim(-3.38,1.56)+
	theme_classic()+ 
	theme(
	  panel.grid = element_blank(),
	  panel.background = element_blank(),
	  axis.text = element_text(colour = "black"),
	  #axis.text.x = element_text(size=8, angle = 90, hjust = 1, vjust = 0.5),
	  axis.text.x = element_blank(),
	  axis.ticks.x = element_blank(),
	  axis.title = element_text(colour = "black")
	)

ggsave(paste0(predictionPath,"prediction_",gene,".png"), p1, width = 4, height =5.4, dpi=400, units = "cm",)
}













library(ggplot2)

pro <- read.csv("/Users/rub2/Workspace/figure_SecAct/panel/wet/proliferation.csv")

p1 <- ggplot(pro,aes(x=Day, y=Value, group=Rep, color=Group))+
	geom_point()+
	geom_line(aes(x=Day, group=interaction(Rep, Group), color= Group), alpha=0.5)+
	scale_color_manual(values=c("red","blue"))+
	ylab("Metabolic activity (492/620)")+
	theme_classic()+
	theme(
		legend.position = c(0.3, 0.8)
	)

ggsave("/Users/rub2/Workspace/figure_SecAct/panel/wet/proliferation1.pdf", p1, width = 2.6, height = 2.5)
	
pro4 <- pro[pro[,1]==4,]

p1 <- ggplot(pro4, aes(x=Group, y=Value, group=Rep, color=Group))+
	geom_bar(stat="identity", position="dodge", fill="white", width=0.7)+
	scale_color_manual(values=c("red","blue"))+
	coord_cartesian(ylim = c(min(pro[,4]), NA)) +
	ylab("Metabolic activity (492/620)")+
	theme_classic()+
	theme(
		legend.position = "none"
	)

ggsave("/Users/rub2/Workspace/figure_SecAct/panel/wet/proliferation2.pdf", p1, width = 2.66, height = 2.5)


   

smy_test1 <- read.csv(paste0(predictionPath,"prediction_act_genomeWide_test.csv"),row.names=1,header=T)
smy_test2 <- read.csv(paste0(predictionPath,"prediction_exp_genomeWide_test.csv"),row.names=1,header=T)

genes1 <- rownames(smy_test1)[smy_test1[,"is.signif"]==1]
genes2 <- rownames(smy_test2)[smy_test2[,"is.signif"]==1]

genes_act_neg <- rownames(smy_test1)[smy_test1[,"is.signif"]==1&smy_test1[,"smy.median"]<0]
genes_act_pos <- rownames(smy_test1)[smy_test1[,"is.signif"]==1&smy_test1[,"smy.median"]>0]

shown_genes <- intersect(genes1,genes2)

CIDE <- rio::import_list(paste0(predictionPath,"CIDE.xlsx")) 
CIDE <- CIDE[[2]]
rownames(CIDE) <- CIDE[,1]

table(shown_genes%in%rownames(CIDE))

shown_genes[!shown_genes%in%rownames(CIDE)]

SecAct <- data.frame(
	gene=shown_genes,
	ICB=ifelse(shown_genes%in%genes_act_neg,"Anti-Tumor","Pro-Tumor"),
	Anno=CIDE[match(SecAct[,1],CIDE[,1]),2],
	Note=CIDE[match(SecAct[,1],CIDE[,1]),3]
)

write.table(SecAct,paste0(predictionPath,"SecAct_SP_function.tsv"),quote=F,sep="\t",row.names=F)

}

