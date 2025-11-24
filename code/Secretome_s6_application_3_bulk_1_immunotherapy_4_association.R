source("Secretome_s0_path.R")


inputPath <- paste0(dataAppPath,"Immunotherapy/")
outputPath <- paste0(applicationPath,"Immunotherapy/")



allFiles <- list.files(inputPath)
allFiles_OS <- allFiles[grepl("OS",allFiles)]
allFiles_PFS <- allFiles[grepl("PFS",allFiles)]
allFiles_RFS <- allFiles[grepl("RFS",allFiles)]
allFiles_RECIST <- allFiles[grepl("RECIST",allFiles)]
allFiles_response <- allFiles[grepl("response",allFiles)]
allFiles <- sort(c(allFiles_OS,allFiles_PFS,allFiles_RFS,allFiles_RECIST,allFiles_response))

x1 <- sapply(strsplit(allFiles,".",fixed=T),function(x) return(x[1]))
x2 <- sapply(strsplit(allFiles,".",fixed=T),function(x) return(x[2]))
mat <- data.frame(study=x1, meta=x2)
rownames(mat) <- allFiles

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
	}else if("RFS"%in%mat_sub[,2]){
		datasets <- c(datasets, rownames(mat_sub)[mat_sub[,2]%in%"RFS"])
		next
	}else if("RECIST"%in%mat_sub[,2]){
		datasets <- c(datasets, rownames(mat_sub)[mat_sub[,2]%in%"RECIST"])
		next
	}else{
		datasets <- c(datasets, rownames(mat_sub)[mat_sub[,2]%in%c("response","response_Basal","response_Luminal")])
	}
}

# not Solid Tumor; from blood
datasets <- datasets[!datasets%in%"Melanoma_ICB_Lozano2022.response"]

# already have subtype
datasets <- datasets[!datasets%in%"Breast_Pembrolizumab_Wolf2022.response"]

datasets <- c(datasets,"NSCLC_ICB_Ravi2023.Adeno.OS","NSCLC_ICB_Ravi2023.Squamous.OS")
datasets <- datasets[!datasets%in%"NSCLC_ICB_Ravi2023.OS"]

#duplicated with later cohort
datasets <- datasets[!datasets%in%"Urothelial_Atezolizumab_Snyder2017.OS"] # Urothelial_Atezolizumab_Mariathasan2018

#low exon coverage
datasets <- datasets[!datasets%in%"Colorectal_ICB_Thibaudin2023.RECIST"]
datasets <- datasets[!datasets%in%"Esophageal_Atezolizumab_VanDenEnde2021.response"]

#low sample number 
datasets <- datasets[!datasets%in%"NSCLC_PD1_Yan2025.RECIST"]
datasets <- datasets[!datasets%in%"Melanoma_CTLA4_Campbell2023.RECIST"]
datasets <- datasets[!datasets%in%"Melanoma_PD1-to-CTLA4_Campbell2023.RECIST"]
datasets <- datasets[!datasets%in%"Melanoma_CTLA4_Roh2017.response"] # targeted
datasets <- datasets[!datasets%in%"Melanoma_PD1_Roh2017.response"] # targeted

# Chemotherapy
datasets <- datasets[!datasets%in%"Urothelial_Atezolizumab_Hamidi2024-IMvigor010.OS"]
datasets <- datasets[!datasets%in%"Urothelial_Atezolizumab_Hamidi2024-IMvigor210.OS"]
datasets <- datasets[!datasets%in%"Urothelial_Chemotherapy_Hamidi2024-IMvigor130.OS"]
datasets <- datasets[!datasets%in%"Urothelial_Chemotherapy_Hamidi2024-IMvigor211.OS"]
datasets <- datasets[!datasets%in%"Urothelial_Surveillance_Hamidi2024-IMvigor010.OS"]
datasets <- datasets[!datasets%in%"Urothelial_Observation_Powles2024.OS"]
datasets <- datasets[!datasets%in%"SCLC_Placebo_Nabet2024.OS"]

datasets <- datasets[!datasets%in%"NSCLC_Docetaxel_Patil2022-OAK.OS"]
datasets <- datasets[!datasets%in%"NSCLC_Docetaxel_Patil2022-POPLAR.OS"]

datasets <- datasets[!datasets%in%"RCC_Sunitinib_Motzer2020-IMmotion151.PFS"]
datasets <- datasets[!datasets%in%"RCC_Sunitinib_Motzer2020.PFS"]
datasets <- datasets[!datasets%in%"mRCC_Sunitinib_McDermott2018.PFS"]

datasets <- datasets[!datasets%in%"Colorectal_Regorafenib_Eng2019.OS"]

datasets <- datasets[!datasets%in%"CCRCC_Everolimus_Braun2020.OS"]
datasets <- datasets[!datasets%in%"Hepatocellular_Sorafenib_Finn2020.OS"]
                  
# Company
datasets <- datasets[!datasets%in%"Pediatric_Atezolizumab_Nabbi2023.RECIST"]
datasets <- datasets[!datasets%in%"Colorectal_Atezolizumab_Eng2019.OS"]
datasets <- datasets[!datasets%in%"Colorectal_Atezo+Cobimetinib_Eng2019.OS"]

# treatment failed
datasets <- datasets[!datasets%in%"Pancreatic_Nivolumab_Padron2022.OS_Sotiga+Chemo"]


# do not have expression
# "HeadNeck_Pembrolizumab_Cristescu2018", # with CNV and mutation
# "Melanoma_Pembrolizumab_Cristescu2018", # with CNV and mutation
# "PanCancer_Pembrolizumab_Cristescu2018", # with CNV and mutation
# "Melanoma_CTLA4_Snyder2014", # with CNV and mutation
# "NSCLC_Pembrolizumab_Rizvi2015" # with CNV and mutation

datasets <- sort(datasets)

datasets


#SWARM -t 2 -g 20 --time 2:00:00
args = commandArgs(trailingOnly=TRUE)
n <- args[1]


for(i in 1:n_downsampling_rep)
{
outputPath <- paste0(applicationPath,"Immunotherapy/downsampling/",n,"/")

for(riskType in c("act"))
{
	smy <- data.frame()
	smy2 <- data.frame()
	ctype_vec <- c()
	for(dataset in datasets)
	{
		if(!dataset%in%c("NSCLC_ICB_Ravi2023.Adeno.OS","NSCLC_ICB_Ravi2023.Squamous.OS" ))
		{
			study <- strsplit(dataset,".",fixed=T) [[1]][1]
			clini <- strsplit(dataset,".",fixed=T) [[1]][2]
			ctype <- strsplit(clini,"_",fixed=T) [[1]][1]
			surv <- read.csv(gzfile(paste0(inputPath,dataset)),as.is=T,sep="\t",check.names=F)
		}else{
			study <- strsplit(dataset,".O",fixed=T) [[1]][1]
			ctype <- strsplit(dataset,".",fixed=T) [[1]][3]
			surv <- read.csv(gzfile(paste0(inputPath,"NSCLC_ICB_Ravi2023.OS")),as.is=T,sep="\t",check.names=F)
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
		if(ctype%in%c("RFS"))
		{
			colnames(surv)[1] <- "RFS"
			colnames(surv)[2] <- "RFS.Event"	
		}
		if(ctype%in%c("RECIST"))
		{
			colnames(surv)[1] <- "RECIST"
		}
		
		ctype_vec <- c(ctype_vec, ctype)
		
		
		if(riskType=="act")
		{
			if(!file.exists(paste0(outputPath,"rep_",i,"_",study,".RData"))) next
			load(paste0(outputPath,"rep_",i,"_",study,".RData"))
	
			Act <- as.matrix(res$zscore)
			Act <- t(Act)
			
			matx <- Act
		}else if(riskType=="exp"){
			load(paste0(outputPath,study,".RData"))
			Act <- as.matrix(res$zscore)
			Act <- t(Act)
						
			Exp <- t(as.matrix(read.csv(paste0(outputPath,study,".csv.gz"),check.names=F,header=TRUE,row.names=1)	))
			
			matx <- Exp[,colnames(Exp)%in%colnames(Act)]
		}else{
			load(paste0(outputPath,study,".RData"))
			Act <- as.matrix(res$zscore)
			Act <- t(Act)
			
			Exp <- t(as.matrix(read.csv(paste0(outputPath,study,".csv.gz"),check.names=F,header=TRUE,row.names=1)	))
			
			# download from https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/h.all.v2025.1.Hs.symbols.gmt
			NFKB_pathway <- c("ABCA1","ACKR3","AREG","ATF3","ATP2B1","B4GALT1","B4GALT5","BCL2A1","BCL3","BCL6","BHLHE40","BIRC2","BIRC3","BMP2","BTG1","BTG2","BTG3","CCL2","CCL20","CCL4","CCL5","CCN1","CCND1","CCNL1","CCRL2","CD44","CD69","CD80","CD83","CDKN1A","CEBPB","CEBPD","CFLAR","CLCF1","CSF1","CSF2","CXCL1","CXCL10","CXCL11","CXCL2","CXCL3","CXCL6","DENND5A","DNAJB4","DRAM1","DUSP1","DUSP2","DUSP4","DUSP5","EDN1","EFNA1","EGR1","EGR2","EGR3","EHD1","EIF1","ETS2","F2RL1","F3","FJX1","FOS","FOSB","FOSL1","FOSL2","FUT4","G0S2","GADD45A","GADD45B","GCH1","GEM","GFPT2","GPR183","HBEGF","HES1","ICAM1","ICOSLG","ID2","IER2","IER3","IER5","IFIH1","IFIT2","IFNGR2","IL12B","IL15RA","IL18","IL1A","IL1B","IL23A","IL6","IL6ST","IL7R","INHBA","IRF1","IRS2","JAG1","JUN","JUNB","KDM6B","KLF10","KLF2","KLF4","KLF6","KLF9","KYNU","LAMB3","LDLR","LIF","LITAF","MAFF","MAP2K3","MAP3K8","MARCKS","MCL1","MSC","MXD1","MYC","NAMPT","NFAT5","NFE2L2","NFIL3","NFKB1","NFKB2","NFKBIA","NFKBIE","NINJ1","NR4A1","NR4A2","NR4A3","OLR1","PANX1","PDE4B","PDLIM5","PER1","PFKFB3","PHLDA1","PHLDA2","PLAU","PLAUR","PLEK","PLK2","PLPP3","PMEPA1","PNRC1","PPP1R15A","PTGER4","PTGS2","PTPRE","PTX3","RCAN1","REL","RELA","RELB","RHOB","RIGI","RIPK2","RNF19B","SAT1","SDC4","SERPINB2","SERPINB8","SERPINE1","SGK1","SIK1","SLC16A6","SLC2A3","SLC2A6","SMAD3","SNN","SOCS3","SOD2","SPHK1","SPSB1","SQSTM1","STAT5A","TANK","TAP1","TGIF1","TIPARP","TLR2","TNC","TNF","TNFAIP2","TNFAIP3","TNFAIP6","TNFAIP8","TNFRSF9","TNFSF9","TNIP1","TNIP2","TRAF1","TRIB1","TRIP10","TSC22D1","TUBB2A","VEGFA","YRDC","ZBTB10","ZC3H12A","ZFP36")
			NFKB_pathway <- transferSymbol(NFKB_pathway)
			gmt <- list(NFKB_pathway=NFKB_pathway)
			gsva_scores <- GSVA::gsva(t(Exp), gmt, method = "gsva", kcdf = "Gaussian")
			
			# combine
			colnames(Exp) <- paste0(colnames(Exp),"_exp")
			Exp <- cbind(Exp,t(gsva_scores) )
			Exp <- Exp[,colnames(Exp)%in%c("LY86_exp","CD180_exp","NFKB_pathway"),drop=F]	
			
			matx <- cbind(LY86_act=Act[,"LY86"], Exp)
			
			smy2[colnames(Exp),dataset] <- cor(Act[,"LY86",drop=F], Exp)[1,]
		}
				
		olp <- intersect(rownames(surv),rownames(matx))
		X_olp <- surv[olp,,drop=F]
		Y_olp <- matx[olp,,drop=F]
		
		
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
				
				}else if(ctype=="RFS"){
					errflag <- F
					coxmodel_fit <- tryCatch(
						coxph(Surv(RFS, RFS.Event) ~ ., data = comb),
						
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
				write.table(Y_olp,paste0(outputPath,"X"),quote=F,sep="\t")
				write.table(X_olp[,"Response",drop=F],paste0(outputPath,"Y"),quote=F,sep="\t")
				write.table(X_olp[,2:ncol(X_olp),drop=F],paste0(outputPath,"B"),quote=F,sep="\t")
				system(paste0("/data/rub2/software/logis/bin/logis_batch -B ",outputPath,"B -X ",outputPath,"X -Y ",outputPath,"Y -cntthres 3 -out ",outputPath,"output"))
				
				zfile <- read.table(paste0(outputPath,"output.zscore"),sep="\t",check.names=F)	
				
				smy[rownames(zfile),dataset] <- - zfile[,1]
				
				system(paste0("rm ",outputPath,"X"))
				system(paste0("rm ",outputPath,"Y"))
				system(paste0("rm ",outputPath,"B"))
				system(paste0("rm ",outputPath,"output.coef"))
				system(paste0("rm ",outputPath,"output.zscore"))
				system(paste0("rm ",outputPath,"output.pvalue"))
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
	
	write.csv(smy,paste0(outputPath,"prediction_",riskType,"_",i,".csv"),quote=F)
	
	if(riskType=="downstream") write.csv(smy2,paste0(outputPath,"prediction_",riskType,"_cor.csv"),quote=F)
}	


} #i
