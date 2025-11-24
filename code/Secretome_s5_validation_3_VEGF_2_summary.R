source("Secretome_s0_path.R")

inputPath <- paste0(dataValPath,"VEGF/")
outputPath <- paste0(validationPath,"VEGF/")

compSigs <- paste0(sigNames,"LigandExp","ReceptorExp","LRsumExp")
items <- c("VEGFA_GSE72951","VEGFA_E-MTAB-3267")


clinicalComb <- data.frame()
compSigs <- c(sigNames,"LigandExp","ReceptorExp","LRsumExp")

for(compSig in compSigs)
{	
	for(item in items)
	{
		if(file.exists(paste0(outputPath,item,"/",compSig,".RData")))
		{
			load(paste0(outputPath,item,"/",compSig,".RData"))
			
			clinicalComb[paste0(item,"_",compSig),"item"] <- item
			clinicalComb[paste0(item,"_",compSig),"compSig"] <- compSig
			clinicalComb[paste0(item,"_",compSig),"value"] <- z
		}
	}
}

mat = reshape2::dcast( clinicalComb , compSig~item )
rownames(mat) <- mat[,1]
mat <- t(mat[,-1])


mat <- mat[order(apply(mat,1,function(x) median(x,na.rm=T))),]
mat <- mat[,order(apply(mat,2,function(x) median(x,na.rm=T)))]

write.csv(t(mat),paste0(outputPath,"validation_VEGF_summary_heatmap.csv"),quote=F)

png(paste0(outputPath,"validation_VEGF_summary_heatmap.png"), width = 30, height = 10, res=200, units = "cm")

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
ggsave(paste0(outputPath,"validation_VEGF_summary_bar.png"), p1, width = 12.5, height = 8, dpi=500, units = "cm", limitsize =FALSE)


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
ggsave(paste0(outputPath,"validation_VEGF_summary_bar2.png"), p1, width = 7.2, height = 7.8, dpi=500, units = "cm", limitsize =FALSE)






sigName <- "SecAct"

for(item in items)
{
	load(paste0(outputPath,"/",item,"/",sigName,"_res.RData"))
	Act <- as.matrix(res$zscore)
	data <- t(as.matrix(Act)) [,c("LY86","VEGFA")]
	
	if(item=="VEGFA_GSE72951")
	{
		survival <- read.csv(paste0(inputPath,"GSE72951.OS.Bevacizumab"),as.is=T,sep="\t")
		xtext <- "Overall (Months)"
		widthValue <- 9.7
		addtext <- 20
	}else{
		survival <- read.csv(paste0(inputPath,"E-MTAB-3267.PFS"),as.is=T,sep="\t")
		xtext <- "Progression-Free (Months)"
		widthValue <- 9.5
		addtext <- 50
	}
	
	
	library(survival)

	out <- run_CoxPH_best_separation(data, survival, margin=5)
	data <- out[[1]]
	survival <- out[[2]]
	result <- out[[3]]
	cutoff <- round(result["VEGFA","thres.opt"],3)

		
	olp <- intersect(rownames(survival),rownames(data))
	X_olp <- survival[olp,,drop=F]
	Y_olp <- data[olp,,drop=F]
	
	comb <- cbind(X_olp, Act=Y_olp[,"VEGFA"])
	comb <- cbind(comb, group=Y_olp[,"VEGFA"]>cutoff)
	write.csv(comb,paste0(outputPath,item,"_",sigName,".csv"),quote=F)

	
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
	
	
	jpeg(file=paste0(outputPath,item,"_",sigName,".jpg"),width=widthValue, height=11, units="cm",res=300)
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
	
	
	library(survminer)
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

	ggsave(paste0(outputPath,item,"_",sigName,".pdf"), p2, width = 3.2, height = 3.5)
}


