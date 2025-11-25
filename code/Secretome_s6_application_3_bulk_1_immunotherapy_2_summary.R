source("Secretome_s0_path.R")


inputPath <- paste0(dataAppPath,"Immunotherapy/")
outputPath <- paste0(applicationPath,"Immunotherapy/")

rename <- read.csv(paste0(dataAppPath,"cohort_rename.txt"),header=T,sep="\t")


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
	
	#Exp <- read.table(paste0(outputPath,study,".diff"),sep="\t",check.names=F)
	#print(nrow(Exp))
	#patient_count <- c(patient_count, dim(res$zscore)[2] )

	load(paste0(outputPath,study,".RData"))
	
	Act <- as.matrix(res$zscore)	
	Act <- expand_rows(Act)

	matx <- t(Act)
	olp <- intersect(rownames(surv),rownames(matx))
	
	patient_count <- c(patient_count, length(olp))
	
	coverage_vec <- c(coverage_vec, ifelse(dataset%in%notGenomeWide_alt, "Targeted", "Genome wide") )
}	


stat <- data.frame(
	ID=datasets, 
	Name=rename[match(datasets,rename[,1]),2],
	clinical=ctype_vec, 
	expr_count=patient_count, 
	Coverage=coverage_vec)

stat <- stat[order(stat[,"Coverage"]),]

sum(patient_count)
#5174

write.table(stat,paste0(outputPath,"patient_stat.txt"),quote=F,sep="\t")





###################
### risk score 
###################

for(riskType in c("exp","act","downstream"))
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
			load(paste0(outputPath,study,".RData"))
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
	
	write.csv(smy,paste0(outputPath,"prediction_",riskType,".csv"),quote=F)
	
	if(riskType=="downstream") write.csv(smy2,paste0(outputPath,"prediction_",riskType,"_cor.csv"),quote=F)
}	


for(riskType in c("act","exp","downstream"))
{
	smy <- read.csv(paste0(outputPath,"prediction_",riskType,".csv"),as.is=T,row.names=1,header=T,check.names=F)
	
	for(vers in c("genomeWide","all"))
	{
		if(vers == "genomeWide")
		{
			smy_v <- smy[, !colnames(smy)%in%notGenomeWide_alt ]
			smy_v <- smy_v[apply(smy_v, 1, function(x) sum(is.na(x)))<2,]
	
			ctype_v <- ctype_vec[!colnames(smy)%in%notGenomeWide_alt]
		}else{
			smy_v <- smy
			ctype_v <- ctype_vec
		}	
	
	
		smy.median <- sort(apply(smy_v,1,function(x) median(x,na.rm=T)))
		
		smy_v <- smy_v[names(smy.median),]
		write.csv(smy_v,paste0(outputPath,"prediction_",riskType,"_",vers,".csv"),quote=F)
		
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
		
		write.csv(smy_test,paste0(outputPath,"prediction_",riskType,"_",vers,"_test.csv"),quote=F)
		
		
		smy_top <- as.matrix(smy_v[smy_test[,"is.signif"]==TRUE,])
		
		
		geneCol <- ifelse(rownames(smy_top)%in%c("CR1L","AOAH","LY86","COLQ","ADAMTS7"),"green","black")
		dataCol <- ctype_v
		dataCol[dataCol=="OS"] <- "black"
		dataCol[dataCol=="PFS"] <- "red"
		dataCol[dataCol=="RFS"] <- "orange"
		dataCol[dataCol=="RECIST"] <- "cyan"
		dataCol[dataCol=="response"] <- "blue"
		
		png(paste0(outputPath,"prediction_",riskType,"_",vers,"_test_signif.png"), width = 28, height = 120, res=200, units = "cm")
		
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
}



###########
# fig s10 b
###########

riskType <- "downstream"
mycolors <- c("#C77CFF","#F8766D","#7CAD00","#00BFC4")
names(mycolors) <- c("LY86_activity","LY86_expression","CD180_expression","NFKB_pathway_activity")

smy <- read.csv(paste0(outputPath,"prediction_",riskType,".csv"),row.names=1)
rownames(smy)[rownames(smy)=="LY86_act"] <- "LY86_activity"
rownames(smy)[rownames(smy)=="LY86_exp"] <- "LY86_expression"
rownames(smy)[rownames(smy)=="CD180_exp"] <- "CD180_expression"
rownames(smy)[rownames(smy)=="NFKB_pathway"] <- "NFKB_pathway_activity"

smy <- as.matrix(smy[, setdiff(colnames(smy),notGenomeWide_alt) ])
smy <- smy[1:4,]

fg.df <- reshape::melt(smy)

sigOrder <- data.frame()
for(i in unique(fg.df[,"X1"]))
{
	sigOrder[i,"r_mean"] <- median(fg.df[fg.df[,"X1"]==i,"value"],na.rm=TRUE)
}
sigOrder <- sigOrder[order(sigOrder[,"r_mean"]),,drop=F]

fg.df[,1] <- factor(fg.df[,1], levels=rownames(sigOrder))


library(ggplot2)
p1 <- ggplot(fg.df,aes(x=X1,y=value)) + 
	geom_hline(yintercept=0, color = "grey", linewidth=0.6, linetype="dashed")+
	geom_jitter(aes(color=X1),alpha=0.5, size=1, width=0.1)+
	geom_boxplot(color="black", alpha=0, width=0.4, outlier.shape = NA)+
	scale_color_manual(values=mycolors[rownames(sigOrder)])+
	ylab("Risk score z")+
	xlab("")+
	theme_classic()+ 
	theme(
		panel.background = element_blank(),
		panel.grid = element_blank(),
		axis.title = element_text(colour = "black"),
		axis.text = element_text(colour = "black"),
		axis.title.y = element_text(size=10),
		axis.text.x = element_text(size=9, angle = 90, hjust = 1, vjust = 0.5),
		axis.text.y = element_text(size=8),
		#axis.ticks = element_line(colour = "darkgrey"),
		#axis.line = element_line(colour = "darkgrey"),
		axis.title.x = element_blank(),
		legend.position="none"
	)
	
ggsave(paste0(outputPath,"prediction_LY86_NFKB_z.png"), p1, width = 5.5, height = 10.5, dpi=200, units = "cm")
write.csv(fg.df, paste0(outputPath,"prediction_LY86_NFKB_z.csv"), quote=FALSE)



smy2 <- read.csv(paste0(outputPath,"prediction_",riskType,"_cor.csv"),row.names=1)
rownames(smy2)[rownames(smy2)=="LY86_exp"] <- "LY86_expression"
rownames(smy2)[rownames(smy2)=="CD180_exp"] <- "CD180_expression"
rownames(smy2)[rownames(smy2)=="NFKB_pathway"] <- "NFKB_pathway_activity"

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
p2 <- ggplot(fg.df,aes(x=X1,y=value)) + 
	geom_hline(yintercept=0, color = "grey", linewidth=0.6, linetype="dashed")+
	geom_jitter(aes(color=X1),alpha=0.5, size=1, width=0.1)+
	geom_boxplot(color="black", alpha=0, width=0.4, outlier.shape = NA)+
	scale_color_manual(values=mycolors[rownames(sigOrder)])+
	ylab("Correlation r ")+
	xlab("")+
	theme_classic()+ 
	theme(
		panel.background = element_blank(),
		panel.grid = element_blank(),
		axis.title = element_text(colour = "black"),
		axis.title.y = element_text(size=10),
		axis.text = element_text(colour = "black"),
		axis.text.x = element_text(size=9, angle = 90, hjust = 1, vjust = 0.5),
		axis.text.y = element_text(size=8),
		#axis.ticks = element_line(colour = "darkgrey"),
		#axis.line = element_line(colour = "darkgrey"),
		axis.title.x = element_blank(),
		legend.position="none"
	)
	
ggsave(paste0(outputPath,"prediction_LY86_NFKB_cor.png"), p2, width = 4.7, height = 10.5, dpi=200, units = "cm")
write.csv(fg.df, paste0(outputPath,"prediction_LY86_NFKB_cor.csv"), quote=FALSE)



##############################
# risk score water fall plot
##############################

gene <- "LY86"
smy <- read.csv(paste0(outputPath,"prediction_act.csv"),as.is=T,row.names=1,header=T)
smy <- as.matrix(smy[, setdiff(colnames(smy),notGenomeWide_alt) ])

smy_test1 <- read.csv(paste0(outputPath,"prediction_act_genomeWide_test.csv"),row.names=1,header=T)
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

ggsave(paste0(outputPath,"prediction_",gene,".png"), p1, width = 5.6, height =5.5, dpi=400, units = "cm",)
ggsave(paste0(outputPath,"prediction_",gene,".pdf"), p1, width = 6.8, height =6.8, dpi=400, units = "cm",)



######
# fig 6d  risk score box plot 
######

library(ggplot2)
p1 <- ggplot(fg.df,aes(x=gene,y=value)) + 
	geom_hline(yintercept=0, color = "black", linewidth=0.8, linetype="dashed")+
	geom_jitter(color="darkgrey",alpha=0.5, size=1, width=0.2)+
	geom_boxplot(color="black", alpha=0, width=0.5, outlier.shape = NA)+
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
		axis.ticks.x = element_blank(),
		legend.position="none"
	)
	
ggsave(paste0(outputPath,"prediction_",gene,"_box.pdf"), p1, width = 3.3, height = 5, dpi=400, units = "cm")
write.csv(fg.df, paste0(outputPath,"prediction_",gene,"_box.csv"), quote=FALSE)



######
# fig 6b  risk score box plot 
######

genes <- c("GZMA","IL11")

smy <- read.csv(paste0(outputPath,"prediction_act.csv"),as.is=T,row.names=1,header=T)
smy <- as.matrix(smy[, setdiff(colnames(smy),notGenomeWide_alt) ])

fg.df <- reshape::melt(smy[genes,])
fg.df <- cbind(fg.df, group=fg.df[,1]%in%genes[1:6])

smy_test1 <- read.csv(paste0(outputPath,"prediction_act_genomeWide_test.csv"),row.names=1,header=T)
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
	
ggsave(paste0(outputPath,"prediction_GZMA_IL11.pdf"), p1, width = 5.5, height = 4.8, dpi=500, units = "cm")
write.csv(t(smy[genes,]), paste0(outputPath,"prediction_GZMA_IL11.csv"), quote=FALSE)




######
# survival plot Fig 6a, e, and Fig s11a
######

gene_dataset <- data.frame(
	gene=c("GZMA","IL11","LY86","LY86","LY86"),
	dataset=c(
		"mRCC_Atezo+Bev_McDermott2018.PFS",
		"Hepatocellular_Atezo+Bev_Finn2020.OS",
		"Melanoma_Ipilimumab_VanAllen2015.OS",
		"NSCLC_Atezolizumab_Patil2022-OAK.OS",
		"Urothelial_Atezo+Chemo_Hamidi2024-IMvigor130.OS"
	)
)

library(survival)

for(n in 1:nrow(gene_dataset))
{
	gene <- gene_dataset[n,"gene"]
	dataset <- gene_dataset[n,"dataset"]
	
	study <- strsplit(dataset,".",fixed=T) [[1]][1]
	clini <- strsplit(dataset,".",fixed=T) [[1]][2]
	survival <- read.csv(gzfile(paste0(inputPath,dataset)),as.is=T,sep="\t",check.names=F)
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
	
	
	
	load(paste0(outputPath,study,".RData"))
	Act <- as.matrix(res$zscore)
	Act <- t(Act)
	
	
	data <- Act[,c("AOAH","VEGFA","GZMA","IFNG","IGF2","IL11","LY86")]
	
	
	out <- run_CoxPH_best_separation(data, survival, margin=5)
	data <- out[[1]]
	survival <- out[[2]]
	result <- out[[3]]
	cutoff <- round(result[gene,"thres.opt"],3)
	
	
	olp <- intersect(rownames(survival),rownames(data))
	X_olp <- survival[olp,,drop=F]
	Y_olp <- data[olp,,drop=F]
	
	comb <- cbind(X_olp, Act=Y_olp[,gene])
	comb <- cbind(comb, group=Y_olp[,gene]>cutoff)
	write.csv(comb,paste0(outputPath,gene,"_",dataset,".csv"),quote=FALSE)


	coxmodel_fit <- coxph(Surv(OS, OS.Event) ~ ., data = comb)
	coxmodel_obj <- summary(coxmodel_fit)
	zs <- coxmodel_obj$coefficients["Act","z"]
	pv <- coxmodel_obj$coefficients["Act","Pr(>|z|)"]
	
	
	surv_o <- Surv(survival[,1],survival[,2])
	groups <- as.character(data[,gene]>cutoff)
	surv_c  <- survfit(surv_o ~ groups)
	
	
	library(survminer)
	surv.df <- cbind(survival,groups)
	
	fit <- survfit(Surv(OS, OS.Event) ~ groups, data = surv.df)
	
	p2 <- ggsurvplot(fit, data=surv.df, palette = c("#0000FF","#FF0000"), legend.labs = c(paste0("Low (n=",sum(groups=="FALSE"),")"),paste0("High (n=",sum(groups=="TRUE"),")")) )$plot+
		xlab(xtext)+
		ylab("Fraction")+
		annotate("text", x = max(survival[,1])*0.7, y=0.88, label = paste0("z = ", round(zs,2), "\np = ", signif(pv,2)), size=4 )

	ggsave(paste0(outputPath,gene,"_",dataset,".pdf"), p2, width = 3.1, height = 3.5)
	
}




###################################################
# compare significant SP count from act and expr
###################################################

smy_test1 <- read.csv(paste0(outputPath,"prediction_act_genomeWide_test.csv"),row.names=1,header=T)
smy_test2 <- read.csv(paste0(outputPath,"prediction_exp_genomeWide_test.csv"),row.names=1,header=T)

genes_act <- rownames(smy_test1)[smy_test1[,"is.signif"]==1]
genes_exp <- rownames(smy_test2)[smy_test2[,"is.signif"]==1]

olp <- intersect(genes_act,genes_exp)

length(genes_act)
length(genes_exp)
length(olp)


write.csv(smy_test1[genes_act,], paste0(outputPath,"prediction_act_genomeWide_test_signif.csv"),quote=F)
write.csv(smy_test2[genes_exp,], paste0(outputPath,"prediction_exp_genomeWide_test_signif.csv"),quote=F)

write.csv(smy_test1[rownames(smy_test1)%in%olp,], paste0(outputPath,"prediction_act_genomeWide_test_signif_overlap.csv"),quote=F)
write.csv(smy_test2[rownames(smy_test2)%in%olp,], paste0(outputPath,"prediction_exp_genomeWide_test_signif_overlap.csv"),quote=F)



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


##########
# fig 6c #
##########

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

ggsave(paste0(outputPath,"prediction_act_exp_count_compare.pdf"), p1, width = 5.3, height = 4.8, units = "cm")
write.csv(fg.df, paste0(outputPath,"prediction_act_exp_count_compare.csv"),quote=F)






# same anti and pro number
genes1_top <- c(genes_act_neg[1:length(genes_exp_neg)],rev(genes_act_pos[1:length(genes_exp_pos)]))

length(unique(c(genes1_top,genes_exp)))

# same total number
#smy_test1_alt <- smy_test1
#smy_test1_alt[,1] <- abs(smy_test1_alt[,1])
#smy_test1_alt <- smy_test1_alt[order(smy_test1_alt[,1],decreasing=T),]
#genes1_top_alt <- rownames(smy_test1_alt)[1:length(genes_exp)]


to_be_annotated <- unique(c(genes1_top, genes_exp, genes_olp_pos, genes_olp_neg))




library(rio) 
xlsx <- import_list(paste0(outputPath,"soluble_factors.doubleblind.Lanqi.xlsx")) 

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



# sheet e
write.csv(comb[genes1_top, 2:3], paste0(outputPath,"prediction_act_genomeWide_test_signif_top_annotation.csv"))
# sheet f
write.csv(comb[genes_exp, 2:3], paste0(outputPath,"prediction_exp_genomeWide_test_signif_annotation.csv"))


smy_test1 <- smy_test1[rev(rownames(smy_test1)),]
smy_test2 <- smy_test2[rev(rownames(smy_test2)),]
genes_olp <- c( genes_olp_neg, rev(genes_olp_pos) )


lqlist <- read.csv(paste0(outputPath,"list_annotation.csv"),row.names=1,header=T)
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

# sheet c
write.csv(olp_annotation, paste0(outputPath,"prediction_olp_genomeWide_test_signif_annotation.csv"))





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
		
		
		png(paste0(outputPath,i,"_",j,"_ROC.png"), width = 8, height = 9.5, res=400, units = "cm")
		
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
		ggsave(paste0(outputPath,i,"_",j,"_AUC.png"), p2, width = 8, height = 8, dpi=200, units = "cm",limitsize = FALSE)

	}
}

smy <- smy[c(1,4,2,5,3,6),]

write.csv(smy, paste0(outputPath,"performance_compare_act_exp.csv"),quote=FALSE)


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
	
ggsave(paste0(outputPath,"prediction_act_exp_accuracy_compare.png"), p0, width = 5.3, height = 7, dpi=400, units = "cm")
write.csv(fg.df, paste0(outputPath,"prediction_act_exp_accuracy_compare.csv"),quote=FALSE)





################ 
# long heatmap
################


smy1 <- read.csv(paste0(outputPath,"prediction_act.csv"),row.names=1,header=T,check.names=F)
smy2 <- read.csv(paste0(outputPath,"prediction_exp.csv"),row.names=1,header=T,check.names=F)

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


write.csv(smy1[olp,], paste0(outputPath,"prediction_together_heatmap_comb_act.csv"),quote=FALSE)
write.csv(smy2[olp,], paste0(outputPath,"prediction_together_heatmap_comb_exp.csv"),quote=FALSE)



smy1 <- read.csv(paste0(outputPath,"prediction_act.csv"),row.names=1,header=T,check.names=F)
smy2 <- read.csv(paste0(outputPath,"prediction_exp.csv"),row.names=1,header=T,check.names=F)

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

ggsave(paste0(outputPath,"prediction_together_heatmap_l.pdf"), p1, width = 16, height = 42, units = "cm")




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

ggsave(paste0(outputPath,"prediction_together_heatmap_s.pdf"), p2, width = 4.4, height = 42, units = "cm")



################ 
# rank
################

smy <- read.csv(paste0(outputPath,"prediction_act.csv"),row.names=1,header=T,check.names=F)
smy.median <- apply(smy,1,function(x) median(x,na.rm=T))

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

ggsave(paste0(outputPath,"prediction_together_heatmap_rank_l.pdf"), p11, width = 4.5, height =32, units = "cm",)


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

ggsave(paste0(outputPath,"prediction_together_heatmap_rank_s.pdf"), p22, width = 4.5, height =32, units = "cm",)





################ 
# annotation
################


lqlist <- read.csv(paste0(outputPath,"list_annotation.csv"),row.names=1,header=T)

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
ggsave(paste0(outputPath,"prediction_together_anno.pdf"), p3, width = 8.6, height = 40, units = "cm")



olp_annotation <- read.csv(paste0(outputPath,"prediction_olp_genomeWide_test_signif_annotation.csv"),row.names=1,as.is=T)
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

ggsave(paste0(outputPath,"prediction_together_heatmap_comb.pdf"), p_comb, width = 27, height =45, units = "cm")
ggsave(paste0(outputPath,"prediction_together_heatmap_comb.png"), p_comb, width = 27, height =45, units = "cm")






###################################################
# compare SP activity from genome-wide and targeted
###################################################


smy_test1 <- read.csv(paste0(outputPath,"prediction_act_genomeWide_test.csv"),row.names=1,header=T)

genes_act_neg <- rownames(smy_test1)[smy_test1[,"is.signif"]==1&smy_test1[,"smy.median"]<0]
genes_act_pos <- rev(rownames(smy_test1)[smy_test1[,"is.signif"]==1&smy_test1[,"smy.median"]>0])


smy <- read.csv(paste0(outputPath,"prediction_act.csv"),as.is=T,row.names=1,header=T)


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
ggsave(paste0(outputPath,"prediction_genome_targeted_compare.pdf"), p1, width = 6, height = 7, dpi=400, units = "cm")
write.csv(fg.df, paste0(outputPath,"prediction_genome_targeted_compare.csv"),quote=FALSE)




smy <- data.frame()
for(cancer in notGenomeWide)
{
	expr <- as.matrix(read.csv(gzfile(paste0(inputPath,cancer,".expression.gz")),as.is=T,sep="\t",check.names=F))
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
ggsave(paste0(outputPath,"Targeted_gene_count.pdf"), p1, width = 7, height = 10, dpi=500, units = "cm", limitsize =FALSE)
write.csv(fg.df, paste0(outputPath,"Targeted_gene_count.csv"),quote=FALSE)

