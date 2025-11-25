source("Secretome_s0_path.R")

version <- "vst"

##########
# fig 1b #
##########

sts <- "SKCM_2022_Sudmeier"
#for(st in sts)
#{
	visiumPath.st <- paste0(visiumPath,st,"/")
	dir.create(visiumPath.st)
	
	signaturePath.st <- paste0(signaturePath,st,"/")
	dir.create(signaturePath.st)
	
	sampleNames <- meta[meta[,"Study"]==st,"Sample_Name"]
	#for(sampleName in sampleNames)
	#{	
		visiumPath.st.sample <- paste0(visiumPath.st,sampleName,"/")
		dir.create(visiumPath.st.sample)
		
		signaturePath.st.sample <- paste0(signaturePath.st,sampleName,"/")
		dir.create(signaturePath.st.sample)
		
		st.matrix.data <- as.matrix(read.table(gzfile(paste0(visiumPath.st.sample,"counts.tsv.gz")),check.names=FALSE))
		rownames(st.matrix.data) <- transferSymbol(rownames(st.matrix.data))
		st.matrix.data <- rm_duplicates(st.matrix.data)
		st.matrix.data <- st.matrix.data[,colSums(st.matrix.data)>100]
	
	#}
#}
fg.df <- t(st.matrix.data[c("TGFB1","SERPINE1","LTBP2"),])

write.csv(fg.df,paste0(fig1Path,"TGFB1_SERPINE1_LTBP2.csv"), quote=FALSE)








##########
# fig 1d #
##########

# TGFB1 to targets or partners
gene <- "TGFB1"
targets <- c("SERPINE1","CCN2","TGFBI")
mediators <- c("LTBP2","TGFBR2","ITGAV")

# all TGFB1 signature from different visium
comb <- data.frame()
sts <- unique(meta[,"Study"])
for(st in sts)
{
	signaturePath.st <- paste0(signaturePath,st,"/")
	
	sampleNames <- meta[meta[,"Study"]==st,"Sample_Name"]
	for(sampleName in sampleNames)
	{	
		signaturePath.st.sample <- paste0(signaturePath.st,sampleName,"/")
		dir.create(signaturePath.st.sample)
		
		signaturePath.st.sample.singleSig <- paste0(signaturePath.st.sample,"singleSig_",version,"/")
		dir.create(signaturePath.st.sample.singleSig)
		
		if(file.exists(paste0(signaturePath.st.sample.singleSig,gene,".tsv.gz")))
		{
			m_sub <- read.table(paste0(signaturePath.st.sample.singleSig,gene,".tsv.gz"), sep="\t")
			comb[rownames(m_sub),paste0(st,"@",sampleName)] <- m_sub[,1]
		}
	} #sampleName
} #st


fg.df <- data.frame()
pvs <- c()
for(pt in c(targets,mediators))
{
	temp <- unlist(comb[pt,])
	temp <- temp[!is.na(temp)]
	if(length(temp)==0) next
	
	wilcox_res <- wilcox.test(temp, mu=0)
	pv <- signif(wilcox_res$p.value,2)
	pvs <- c(pvs, pv)
	
	fg.df[paste0(gene,pt,names(temp)),"x"] <- pt
	fg.df[paste0(gene,pt,names(temp)),"y"] <- temp
	fg.df[paste0(gene,pt,names(temp)),"z"] <- gene
}
names(pvs) <- c(targets,mediators)

fg.df[,1] <- factor(fg.df[,1],levels= c(targets, setdiff(names(sort(pvs)),targets) ) )
fg.df <- cbind(fg.df,group="none")
fg.df[fg.df[,1]%in%targets,"group"] <- "target"
fg.df[fg.df[,1]%in%mediators,"group"] <- "mediator"


library(ggplot2)
p1 <- ggplot(fg.df,aes(x=x,y=y,colour=group)) + 
	geom_hline(yintercept=0, color = "grey", linewidth=0.8, linetype="dashed")+
	geom_jitter(alpha=0.3, size=0.2, width=0.25)+
	geom_boxplot(color="black",fill="white", alpha=0.5, width=0.5, outlier.shape = NA)+
	scale_colour_manual(values=c("#41b9c1","#8080FF"))+
	annotate("text", x = 0.75, y=0.2, label = paste0("p = "), size=2.8)+
	annotate("text", x = 1, y=0.16, label = paste0("",pvs[1]), size=2.8)+
	annotate("text", x = 2, y=0.16, label = paste0("",pvs[2]), size=2.8)+
	annotate("text", x = 3, y=0.16, label = paste0("",pvs[3]), size=2.8)+
	annotate("text", x = 4, y=0.16, label = paste0("",pvs[4]), size=2.8)+
	annotate("text", x = 5, y=0.16, label = paste0("",pvs[5]), size=2.8)+
	annotate("text", x = 6, y=0.16, label = paste0("",pvs[6]), size=2.8)+
	coord_cartesian(ylim = c(-0.1, 0.2))+
	ylab("Moran's I")+
	theme_classic()+ 
	theme(
		axis.text.x = element_text(angle = 38, hjust = 1),
		axis.title.x = element_blank(),
		legend.position="none"
	)
ggsave(paste0(fig1Path,"moranI_example_single_",gene,".png"), p1, width = 10.8, height = 5.3, dpi=300, units = "cm",limitsize = FALSE)
write.csv(fg.df,paste0(fig1Path,"moranI_example_single_",gene,".csv"), quote=FALSE)



###########
# fig s1d #
###########


# TGFB1-targets vs TGFB1-random
comb <- comb[rowSums(is.na(comb))<ncol(comb)*0.3,,drop=F]	
comb_vec <- apply(comb,1,function(x) median(x,na.rm=T))

#https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_TGF_BETA_SIGNALING.html
TGFB1_targets <- c("ACVR1","APC","ARID4B","BCAR3","BMP2","BMPR1A","BMPR2","CDH1","CDK9","CDKN1C","CTNNB1","ENG","FKBP1A","FNTA","FURIN","HDAC1","HIPK2","ID1","ID2","ID3","IFNGR2","JUNB","KLF10","LEFTY2","LTBP2","MAP3K7","NCOR2","NOG","PMEPA1","PPM1A","PPP1CA","PPP1R15A","RAB31","RHOA","SERPINE1","SKI","SKIL","SLC20A1","SMAD1","SMAD3","SMAD6","SMAD7","SMURF1","SMURF2","SPTBN1","TGFB1","TGFBR1","TGIF1","THBS1","TJP1","TRIM33","UBE2D3","WWTR1","XIAP")
TGFB1_targets <- unique(transferSymbol(TGFB1_targets))
TGFB1_targets <- intersect(names(comb_vec),TGFB1_targets)

set.seed(123)
TGFB1_nonTargets <- sample(setdiff(names(comb_vec),TGFB1_targets),length(TGFB1_targets))

mean(comb_vec[TGFB1_targets])
mean(comb_vec[TGFB1_nonTargets])
wilcox_res <- wilcox.test(comb_vec[TGFB1_targets], comb_vec[TGFB1_nonTargets])
pv <- signif(wilcox_res$p.value,3)

fg.df <- data.frame(x=c(rep("Targets",length(TGFB1_targets)),rep("Random",length(TGFB1_nonTargets))),y=c(comb_vec[TGFB1_targets], comb_vec[TGFB1_nonTargets]))

library(ggplot2)
p2 <- ggplot(fg.df,aes(x=x, y=y, fill=x)) + 
	geom_hline(yintercept=0, color = "grey", linewidth=0.8, linetype="dashed")+
	geom_jitter(alpha=0.2, size=1, width=0.08)+
	geom_boxplot(colour="grey1", alpha=0.6, outlier.shape = NA)+
	scale_fill_manual(values=c("#8dd3c7","#8080FF"))+
	annotate("text", x = 1.5, y=0.025, label = paste0("p = ",pv), size=4.2)+
	ylab("Median Moran's I")+
	xlab(" ")+
	theme_bw()+ 
	theme(
		plot.background = element_blank(),
		panel.grid = element_blank(),
		axis.title = element_text(size=13,colour = "black"),
		axis.text = element_text(size=12,colour = "black"),
		legend.position="none"
	)

ggsave(paste0(fig1Path,"moranI_example_TGFB1.png"), p2, width = 10, height = 9.5, dpi=400, units = "cm")
write.csv(fg.df,paste0(fig1Path,"moranI_example_TGFB1.csv"),quote=F)




###########
# fig s1e #
###########

fg.df <- data.frame()
sts <- unique(meta[,"Study"])
for(st in sts)
{	
	signaturePath.st <- paste0(signaturePath,st,"/")
	dir.create(signaturePath.st)
		
	sampleNames <- meta[meta[,"Study"]==st,"Sample_Name"]
	for(sampleName in sampleNames)
	{	
		signaturePath.st.sample <- paste0(signaturePath.st,sampleName,"/")
		dir.create(signaturePath.st.sample)
		
		smyMat.median <- readLines(paste0(signaturePath.st.sample,"cor_of_two_vst_version.txt"))

		fg.df[paste0(st,"@",sampleName),"x"] <- "Visium"
		fg.df[paste0(st,"@",sampleName),"y"] <- as.numeric(smyMat.median)

	}
}

st <- "SKCM_2022_Sudmeier"
sampleName <- "pt16"
signaturePath.st.sample <- paste0(signaturePath,st,"/",sampleName,"/")
		
fg.df[["z"]] <- rownames(fg.df)==paste0(st,"@",sampleName)
fg.df <- fg.df[order(fg.df[["z"]]),]		

library(ggplot2)
p1 <- ggplot(fg.df,aes(x=x, y=y)) + 
	geom_hline(yintercept=0, color = "grey", linewidth=0.8, linetype="dashed")+
	geom_jitter(aes(colour=z),size=1)+
	scale_colour_manual(values=c("#fb8072","green"))+
	ylab("Median correlation r")+
	xlab(" ")+
	theme_bw()+ 
	theme(
		plot.background = element_blank(),
		panel.grid = element_blank(),
		axis.title = element_text(colour = "black"),
		axis.text = element_text(colour = "black"),
		axis.text.x = element_blank(), 
		axis.ticks.x = element_blank(), 
		legend.position="none"
	)

ggsave(paste0(fig1Path,"cor_of_two_vst_version_all_visiums.png"), p1, width = 3.5, height = 7, dpi=300, units = "cm")
write.csv(fg.df,paste0(fig1Path,"cor_of_two_vst_version_all_visiums.csv"),quote=F)
	


gene <- "TGFB1"

st.matrix.data.vst.1 <- as.matrix(read.table(gzfile(paste0(signaturePath.st.sample,"vst.tsv.gz")), sep="\t"))
st.matrix.data.vst.2 <- as.matrix(read.table(gzfile(paste0(signaturePath.st.sample,"vst_condition_logUMI_cellType.tsv.gz")), sep="\t"))
smyMat <- cor.two.matrix.same.rows(st.matrix.data.vst.1, st.matrix.data.vst.2)

fg.df <- data.frame(x=st.matrix.data.vst.1[gene,], y=st.matrix.data.vst.2[gene,])
	
p2 <- ggplot(fg.df,aes(x=x, y=y)) + 
	geom_point(color="#bebada")+
	ggtitle(paste0(gene," Expr (transformed)"))+
	ylab("Seq depth")+
	xlab("Seq depth and cell composition")+
	theme_bw()+ 
	theme(
		plot.background = element_blank(),
		panel.grid = element_blank(),
		plot.title = element_text(hjust = 0.5),
		axis.title = element_text(colour = "black"),
		axis.text = element_text(colour = "black"),
		legend.position="none"
	)

ggsave(paste0(fig1Path, "cor_of_two_vst_version_",gene,".png"), p2, width = 8, height = 8, dpi=300, units = "cm")
write.csv(fg.df,paste0(fig1Path, "cor_of_two_vst_version_",gene,".csv"),quote=F,row.names=F)


fg.df <- data.frame(x=sampleName, y=smyMat, z="Yes")
fg.df[["z"]] <- names(smyMat)==gene
fg.df <- fg.df[order(fg.df[["z"]]),]
			
p3 <- ggplot(fg.df,aes(x=x, y=y)) + 
	geom_hline(yintercept=0, color = "grey", linewidth=0.8, linetype="dashed")+
	geom_jitter(aes(colour=z),size=1)+
	scale_colour_manual(values=c("#8dd3c7","blue"))+
	ylab("Correlation r")+
	xlab(" ")+
	theme_bw()+ 
	theme(
		plot.background = element_blank(),
		panel.grid = element_blank(),
		axis.title = element_text(colour = "black"),
		axis.text = element_text(colour = "black"),
		axis.text.x = element_blank(), 
		axis.ticks.x = element_blank(), 
		legend.position="none"
	)
	
ggsave(paste0(fig1Path, "cor_of_two_vst_version_all_genes.png"), p3, width = 3.5, height = 7, dpi=300, units = "cm")
write.csv(fg.df,paste0(fig1Path, "cor_of_two_vst_version_all_genes.csv"),quote=F)




###########
# fig 1e #
###########
gene <- "TGFB1"
st <- "SKCM_2022_Sudmeier"
sampleName <- "pt16"
QCPath.st.sample <- paste0(QCPath,st,"/",sampleName,"/")
QCFilterPath.st.sample <- paste0(QCFilterPath,st,"/",sampleName,"/")
					
fg.df <- read.csv(paste0(QCPath.st.sample,"TCGA_SKCM_",gene,"_expr_vs_signatureSoce.csv"),row.names=1)

dim(fg.df)
## [1] 469   2

cor_res <- cor.test(fg.df[,1],fg.df[,2])
rv <- round(cor_res$estimate,2)
pv <- signif(cor_res$p.value,2)

library(ggplot2)
p1 <- ggplot(fg.df,aes(x=x, y=y)) + 
	geom_point(color="#ff8080", alpha=0.5, size=0.6)+
	annotate("text", x=-2.1, y=0.48, label=paste0("r = ",rv,"\np = ",pv) ,size=3)+
	xlab("TGFB1 Expression")+
	ylab("Signature Score")+
	theme_classic()+ 
	theme(
		panel.grid = element_blank(),
		panel.background = element_blank(),
		axis.title = element_text(colour = "black"),
		axis.text = element_text(colour = "black")
	)#+geom_smooth(method = 'lm',color="brown")		
ggsave(paste0(fig1Path, "QC0_scatter_",gene,"_expr_sigScore.png"), p1, width = 5.6, height = 5.2, dpi=200, units = "cm")
write.csv(fg.df,paste0(fig1Path, "QC0_scatter_",gene,"_expr_sigScore.csv"),quote=F)



###########
# fig 1f #
###########

# real vs random in SKCM
cancer <- "TCGA_SKCM"

comb <- data.frame()
for(group in c("Random","Real"))
{
	if(group=="Random")
	{
		smyMat <- read.csv(paste0(QCPath.st.sample,cancer,"_random.csv"),row.names=1)
	}else{
		smyMat <- read.csv(paste0(QCPath.st.sample,cancer,"_vst.csv"),row.names=1)
	}
	
	smyMat <- smyMat[!is.na(smyMat[,1]),,drop=F]
	
	comb[rownames(smyMat),group] <- smyMat[,1]
}
comb <- comb[rownames(comb)%in%SPs,]

dim(comb)
## [1] 929   2

comb.m <- reshape2::melt(as.matrix(comb))
comb.m <- comb.m[!is.na(comb.m[,3]),]

comb.m[comb.m[,1]=="TGFB1",]
## 	      Var1   Var2       value
## 808  TGFB1 Random -0.03791815
## 1736 TGFB1   Real  0.59266157

x1 <- comb.m[comb.m[,"Var2"]=="Random",3]
x2 <- comb.m[comb.m[,"Var2"]=="Real",3]
pv <- signif(wilcox.test(x1, x2)$p.value,2)
	
library(ggplot2)
p2 <- ggplot(comb.m,aes(x = Var2, y = value, colour=Var2))+
	geom_hline(yintercept=0, color = "grey", linewidth=0.6, linetype="dashed")+
	geom_jitter(alpha=0.3, size=0.2, width=0.25,)+
	geom_boxplot(color="black",fill="white", alpha=0.5, width=0.5, outlier.shape = NA)+
	scale_color_manual(values = c("grey","#ff8080"))+
	annotate("text", x=1.2, y=0.85, label=paste0("p = ",pv) ,size=3)+
	xlab("Signature")+
	ylab("Correlation r")+
	theme_classic()+ 
	theme(
	  panel.grid = element_blank(),
	  panel.background = element_blank(),
	  axis.text = element_text(colour = "black"),
	  axis.title = element_text(colour = "black"),
	  legend.position="none"
	) 
	
ggsave(paste0(fig1Path, "QC0_pt16_SKCM_SPs.png"), p2, width = 5, height = 5, dpi=200, units = "cm",limitsize = FALSE)
write.csv(comb.m,paste0(fig1Path, "QC0_pt16_SKCM_SPs.csv"),quote=F)


	
###########
# fig s2a #
###########

# real vs random in all cancer cohorts
cancers <- list.files(QCPath.st.sample)
cancers <- cancers[grepl("_random",cancers)]
cancers <- sapply(strsplit(cancers,"_random",fixed=T),function(x) return(x[1]))

comb <- data.frame()
for(cancer in cancers)
{
	for(group in c("Random","Real"))
	{
		if(group=="Random")
		{
			smyMat <- read.csv(paste0(QCPath.st.sample,cancer,"_random.csv"),row.names=1)
		}else{
			smyMat <- read.csv(paste0(QCPath.st.sample,cancer,"_vst.csv"),row.names=1)
		}
		
		smyMat <- smyMat[!is.na(smyMat[,1]),,drop=F]
		
		comb[paste0(cancer,group,rownames(smyMat)),"cancer"] <- cancer
		comb[paste0(cancer,group,rownames(smyMat)),"Group"] <- group
		comb[paste0(cancer,group,rownames(smyMat)),"gene"] <- rownames(smyMat)
		comb[paste0(cancer,group,rownames(smyMat)),"value"] <- smyMat[,1]
	}
}

comb.m <- comb[comb[,"gene"]%in%SPs,]


comb.m[,"cancer"] <- factor(
	comb.m[,"cancer"],
	levels=c(cancers[grepl("TCGA",cancers)],cancers[grepl("ICGC",cancers)])
)

library(ggplot2)
library(ggsignif)
p2 <- ggplot(comb.m,aes(x = Group, y = value, colour=Group))+
	geom_hline(yintercept=0, color = "grey", linewidth=0.6, linetype="dashed")+
	geom_jitter(alpha=0.3, size=0.2, width=0.25,)+
	geom_boxplot(color="black",fill="white", alpha=0.5, width=0.5, outlier.shape = NA)+
	geom_signif( comparisons = list(c("Random","Real")) ,y_position = c(0.72), test = "wilcox.test", color="black")+
	scale_color_manual(values = c("grey","#ff8080"))+
	xlab("Signature")+
	ylab("Pearson r value")+
	guides(colour = guide_legend(override.aes = list(size=5)))+
	theme_bw()+ 
	theme(
	  panel.grid = element_blank(),
	  panel.background = element_blank(),
	  axis.text = element_text(size=15,colour = "black"),
	  axis.text.x = element_blank(),
	  axis.ticks.x = element_blank(),
	  axis.title.x = element_blank(),
	  axis.title = element_text(size=15,colour = "black"),
	  legend.position="bottom",
	  legend.text=element_text(size=15),
	  legend.title=element_text(size=15),
	  strip.text = element_text(size = 13)
	)+facet_wrap(~ cancer, ncol =11) #scales = "free"

ggsave(paste0(fig1Path, "QC0_pt16_allCancerType_SPs.png"), p2, width = 55, height = 23, dpi=200, units = "cm", limitsize = FALSE)

write.csv(comb.m,paste0(fig1Path, "QC0_pt16_allCancerType_SPs.csv"),quote=F)

comb.m[,4] <- round(comb.m[,4],3)


mat = reshape2::dcast( comb.m[comb.m[,"Group"]=="Random", ], cancer~gene )
rownames(mat) <- mat[,1]
mat <- t(mat[,-1])

write.csv(mat,paste0(fig1Path, "QC0_pt16_allCancerType_SPs_Random.csv"),quote=F)

mat = reshape2::dcast( comb.m[comb.m[,"Group"]=="Real", ], cancer~gene )
rownames(mat) <- mat[,1]
mat <- t(mat[,-1])

write.csv(mat,paste0(fig1Path, "QC0_pt16_allCancerType_SPs_Real.csv"),quote=F)


###########
# fig 1g #
###########

# Secreted vs Intracellular
cancer <- "TCGA_SKCM"

combComp_r <- read.csv(paste0(QCFilterPath.st.sample,"Secreted_filterBar_vst_comp_r.csv"), row.names=1)

smyMat <- read.csv(paste0(QCPath.st.sample,cancer,"_vst.csv"),row.names=1)
smyMat <- smyMat[!is.na(smyMat[,1]),,drop=F]

smyMat <- cbind(smyMat,Group="NA")
smyMat[rownames(smyMat)%in%IPs,"Group"] <- "Intracellular"
smyMat[rownames(smyMat)%in%MPs,"Group"] <- "Membrane"
smyMat[rownames(smyMat)%in%SPs,"Group"] <- "Secreted"
smyMat <- smyMat[smyMat[,"Group"]!="NA",]

smyMat <- smyMat[smyMat[,"Group"]%in%c("Secreted","Intracellular"),]

table(smyMat[,"Group"])
## Intracellular      Secreted 
##          9340           929
         
x1 <- smyMat[smyMat[,"Group"]=="Secreted",1]
x2 <- smyMat[smyMat[,"Group"]=="Intracellular",1]
pv <- signif(wilcox.test(x1, x2, alternative="greater")$p.value,2)

library(ggplot2)
library(ggsignif)
p2 <- ggplot(smyMat,aes(x = Group, y = x))+
	geom_violin(aes(group=Group, fill=Group),trim=FALSE)+
	scale_fill_manual(values = c("skyblue","#ff8080"))+
	geom_boxplot(width=0.15,outlier.shape=NA)+
	geom_hline(yintercept = combComp_r[cancer,"new_0.25"], colour="black", linewidth=0.6, linetype="dashed")+
	annotate("text", x=1.5, y=1.1, label=paste0("p = ",pv) ,size=3)+
	xlab("Group")+
	ylab("Correlation r")+
	theme_classic()+ 
	theme(
	  panel.grid = element_blank(),
	  panel.background = element_blank(),
	  axis.text = element_text(colour = "black"),
	  axis.title = element_text(colour = "black"),
	  axis.title.x = element_blank(),
	  legend.position="none"
	) 
ggsave(paste0(fig1Path,"QC2_pt16_SKCM_SPs_vs_IPs.png"), p2, width = 6, height = 5.2, dpi=200, units = "cm",limitsize = FALSE)
write.csv(smyMat,paste0(fig1Path, "QC2_pt16_SKCM_SPs_vs_IPs.csv"),quote=F)




###########
# fig 1h #
###########
						
fg.df <- data.frame()

x2_sorted <- sort(x2,decreasing=T)
tempMin <- 1
for(i in length(x2_sorted):1 )
{
	X <- sum( x2 >= x2_sorted[[i]] )/length(x2)
	Y <- sum( x1 >= x2_sorted[[i]] )/length(x1)
	
	fg.df[i,"r"] <- x2_sorted[[i]]
	fg.df[i,"fraction"] <- (length(x2_sorted)-i+1)/length(x2_sorted)
	fg.df[i,"Intracellular"] <- X
	fg.df[i,"Secreted"] <- Y
	fg.df[i,"FDR"] <- min(X/Y,tempMin)
	
	tempMin <- min(X/Y,tempMin)
}

fg.df[,"FDR"] <- 1 - fg.df[,"FDR"]

library(ggplot2)
p2 <- ggplot(fg.df, aes(x = r))+
	geom_point(aes(y = Secreted), color="#ff8080", size=0.1)+
	geom_point(aes(y = Intracellular), color="skyblue", size=0.1)+
	geom_point(aes(y = FDR), color="orange", size=0.1)+
	#scale_colour_manual(values = c("orange","skyblue","#ff8080"))+
	geom_vline(xintercept = combComp_r["TCGA_SKCM","new_0.25"], colour="black", linewidth=0.6, linetype="dashed")+
	xlab("Correlation r")+
	ylab("Fraction")+
	scale_y_continuous(
		  # Features of the first axis
		  name = "Fraction",
		  # Add a second axis and specify its features
		  sec.axis = sec_axis(~1-., name="False Discovery Rate")
		  #trans = "reverse"
		)+
	theme_bw()+ 
	theme(
	  panel.grid = element_blank(),
	  panel.background = element_blank(),
	  panel.border = element_rect(colour = "black"),
	  axis.text = element_text(colour = "black"),
	  axis.title = element_text(colour = "black")
	) 
ggsave(paste0(fig1Path,"QC2_pt16_SKCM_SPs_vs_IPs_FDR_new.png"), p2, width = 7.0, height = 5.2, dpi=200, units = "cm",limitsize = FALSE)
write.csv(fg.df,paste0(fig1Path, "QC2_pt16_SKCM_SPs_vs_IPs_FDR_new.csv"),quote=F)




###########
# fig s2b #
###########

# Secreted vs Intracellular in all cancer cohorts
cancers <- list.files(QCPath.st.sample)
cancers <- cancers[grepl("_random",cancers)]
cancers <- sapply(strsplit(cancers,"_random",fixed=T),function(x) return(x[1]))

average_samples_in_one_cancer <- function(cancer)
{
	combCancer <- data.frame()
	
	sts <- unique(meta[,"Study"])
	for(st in sts)
	{
		QCPath.st <- paste0(QCPath,"/",st,"/")
		dir.create(QCPath.st)
		
		sampleNames <- meta[meta[,"Study"]==st,"Sample_Name"]
		for(sampleName in sampleNames)
		{	
			QCPath.st.sample <- paste0(QCPath.st,sampleName,"/")
			dir.create(QCPath.st.sample)
			
			smyMat <- read.csv(paste0(QCPath.st.sample,cancer,"_vst.csv"),row.names=1)
			smyMat <- smyMat[!is.na(smyMat[,1]),,drop=F]
				
			combCancer[rownames(smyMat),paste0(st,"@",sampleName)] <- smyMat[,1]
		}
	}
	
	combCancer <- apply(combCancer,1,function(x) mean(x,na.rm=T) )
}	

combCancer_list <- parallel::mclapply(cancers, average_samples_in_one_cancer, mc.cores=15) 
names(combCancer_list) <- cancers

comb <- data.frame()
for(cancer in cancers)
{
	comb[names(combCancer_list[[cancer]]),cancer] <- combCancer_list[[cancer]]
}


comb.m <- reshape2::melt(as.matrix(comb))
comb.m <- comb.m[!is.na(comb.m[,3]),]
comb.m <- cbind(comb.m,Group="NA")

comb.m[comb.m[,1]%in%IPs,"Group"] <- "Intracellular"
comb.m[comb.m[,1]%in%MPs,"Group"] <- "Membrane"
comb.m[comb.m[,1]%in%SPs,"Group"] <- "Secreted"

comb.m <- comb.m[!comb.m[,"Group"]%in%c("NA"),]
comb.m <- comb.m[!comb.m[,"Group"]%in%c("Membrane","NA"),]


comb.m[,2] <- factor(
	comb.m[,2],
	levels=c(cancers[grepl("TCGA",cancers)],cancers[grepl("ICGC",cancers)])
)

library(ggplot2)
library(ggsignif)
p2 <- ggplot(comb.m,aes(x = Group, y = value))+
	geom_hline(yintercept=0, color = "grey", linewidth=0.6, linetype="dashed")+
	geom_violin(aes(group=Group, fill=Group), trim=FALSE)+
	geom_boxplot(width=0.1, outlier.shape=NA)+
	#geom_signif( comparisons = list(c("Secreted", "Intracellular"),c("Membrane", "Intracellular")) ,y_position = c(0.72,-0.5), test = "wilcox.test", test.args ="greater")+
	#scale_fill_manual(values = c("skyblue","green","#ff8080"))+
	geom_signif( comparisons = list(c("Secreted", "Intracellular")) ,y_position = c(0.72), test = "wilcox.test", test.args ="greater")+
	scale_fill_manual(values = c("skyblue","#ff8080"))+
	ylab("Pearson r value")+
	theme_bw()+ 
	theme(
	  panel.grid = element_blank(),
	  panel.background = element_blank(),
	  axis.text = element_text(size=15,colour = "black"),
	  axis.text.x = element_blank(),
	  axis.ticks.x = element_blank(),
	  axis.title.x = element_blank(),
	  axis.title = element_text(size=15,colour = "black"),
	  legend.position="bottom",
	  legend.text=element_text(size=15),
	  legend.title=element_text(size=15),
	  strip.text = element_text(size = 13)
	) + facet_wrap(~ Var2, ncol =11) #scales = "free"
	
ggsave(paste0(fig1Path,"QC2_all_SPs_vs_IPs.png"), p2, width = 55, height = 23, dpi=200, units = "cm", limitsize = FALSE)


comb.m[,3] <- round(comb.m[,3],3)

mat = reshape2::dcast( comb.m , Var2~Var1 )
rownames(mat) <- mat[,1]
mat <- t(mat[,-1])

write.csv(mat,paste0(fig1Path, "QC2_all_SPs_vs_IPs.csv"),quote=F)


comb.m[,"flag"] <- 1
comb.m <- comb.m[,c("Var1","Group","flag")]
mat = reshape2::dcast( comb.m , Group~Var1)
rownames(mat) <- mat[,1]
mat <- t(mat[,-1])

mat[mat>0] <- 1
write.csv(mat[,2,drop=FALSE],paste0(fig1Path, "QC2_all_SPs_vs_IPs_anno.csv"),quote=F)

