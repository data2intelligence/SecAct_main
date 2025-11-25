source("Secretome_s0_path.R")	
library(SecAct)
library(SpaCET)


inputPath <- paste0(dataAppPath,"ST_Visium_LIHC/")
outputPath <- paste0(applicationPath,"ST_LIHC/")
dir.create(outputPath)

sampleNames <- c("HCC-R1","HCC-R2")

runParallel <- function(sampleName)
{


# load ST data to create an SpaCET object
visiumPath <- paste0(inputPath,sampleName,"/")
SpaCET_obj <- create.SpaCET.object.10X(visiumPath = visiumPath)

# filter out spots with less than 1000 expressed genes
SpaCET_obj <- SpaCET.quality.control(SpaCET_obj, min.genes = 1000)

# deconvolution
SpaCET_obj <- SpaCET.deconvolution(
  SpaCET_obj, 
  cancerType = "LIHC", 
  coreNo = 10
)



# infer activity
SpaCET_obj <- SecAct.activity.inference.ST(
    inputProfile=SpaCET_obj,
    scale.factor = 1e+05,
)

SpaCET_obj <- SecAct.signaling.pattern(
	SpaCET_obj, 
	scale.factor = 1e+05, 
	radius=200, 
	k=3
)

SpaCET_obj <- SpaCET.GeneSetScore(SpaCET_obj, GeneSets="Hallmark")

saveRDS(SpaCET_obj, file = paste0(outputPath,"/Visium_",sampleName,"_SpaCET_obj.rds"))



# viz cell fraction
SpaCET_obj <- readRDS(paste0(outputPath,"/Visium_",sampleName,"_SpaCET_obj.rds"))

p <- SpaCET.visualize.spatialFeature(
  SpaCET_obj, 
  spatialType = "CellFraction", 
  spatialFeatures= "All", 
  sameScaleForFraction = FALSE,
  imageBg = FALSE,
  legend.position = "none",
  pointSize = 0.5,
  nrow = 5
)

ggsave(paste0(outputPath,"/Visium_",sampleName,"_cellFraction.png"), p, width = 50, height = 30, dpi=500, units = "cm")



p <- SpaCET.visualize.spatialFeature(
  SpaCET_obj, 
  spatialType = "CellFraction", 
  spatialFeatures=list(Malignant="Malignant",Fibroblast="CAF",Endothelial="Endothelial",Macrophage="Macrophage",Hepatocyte="Hepatocyte",'B cell'="B cell", 'T cell'=c("T CD4","T CD8")), 
  sameScaleForFraction = TRUE,
  imageBg = FALSE,
  legend.position = "none",
  pointSize = 0.5, 
  nrow = 1
)
#p <- p & theme(plot.margin = margin(0, 0, 0, 0))
ggsave(paste0(outputPath,"/Visium_",sampleName,"_cellFraction_sub.png"), p, width = 35, height = 5.8, dpi=200, units = "cm")
write.csv(
	t(SpaCET_obj@results$deconvolution$propMat[c("Malignant","CAF","Endothelial","Macrophage","Hepatocyte","B cell","T CD4","T CD8"),]),
	paste0(outputPath,"/Visium_",sampleName,"_cellFraction.csv"),quote=F)



# viz Signaling Pattern
p <- SpaCET.visualize.spatialFeature(
	SpaCET_obj, 
	spatialType = "SignalingPattern", 
	spatialFeatures="All", 
	imageBg=FALSE,
	legend.position = "none",
	colors=c("#03c383","#fbbf45","#ef6a32"),
	pointSize=0.8
)
ggsave(paste0(outputPath,"/Visium_",sampleName,"_pattern_1row.png"), p, width = 25, height = 9, dpi=100, units = "cm")


genes <- c("HDGF","MYDGF","COL1A1","TGFB1","HGFAC","LCAT")
for(gene in genes)
{
	p <- SpaCET.visualize.spatialFeature(
		SpaCET_obj, 
		spatialType = "SecretedProteinActivity", 
		spatialFeatures=gene,
		imageBg=FALSE,
		colors=c("#1A9850","#1A9850","#1A9850","#1A9850","#FDAE61","#FC8D59","#D7191C"),
		legend.position = "none",
		pointSize=0.8
	)
	ggsave(paste0(outputPath,"/Visium_",sampleName,"_activity_",gene,".png"), p, width = 7.3, height = 7.9, dpi=100, units = "cm")
}
write.csv(
	t(SpaCET_obj@results$SecAct_output$SecretedProteinActivity$zscore[genes,]),
	paste0(outputPath,"/Visium_",sampleName,"_activity_6genes.csv"),quote=F)



#library(SecAct)
#SpaCET_obj <- readRDS(paste0("/Users/rub2/Downloads/Visium_HCC-R1_SpaCET_obj.rds"))
p <- SecAct.signaling.velocity.spotST(SpaCET_obj, gene = "SPARC", signalMode="receiving", contourMap=TRUE, coutourBins=9) 
ggsave(paste0(outputPath,"/Visium_",sampleName,"_activity_SPARC.png"), p, width = 22, height = 20, dpi=300, units = "cm")
write.csv(
	t(SpaCET_obj@results$SecAct_output$SecretedProteinActivity$zscore["SPARC",,drop=FALSE]),
	paste0(outputPath,"/Visium_",sampleName,"_activity_SPARC.csv"),quote=F)



p <- SpaCET.visualize.spatialFeature(
  SpaCET_obj, 
  spatialType = "GeneSetScore", 
  spatialFeatures = rownames(SpaCET_obj@results$GeneSetScore),
  legend.position = "none",
  nrow=5
)

ggsave(paste0(outputPath,"/Visium_",sampleName,"_hallmark.png"), p, width = 100, height = 50, dpi=300, units = "cm")



p <- SpaCET.visualize.spatialFeature(
  SpaCET_obj, 
  spatialType = "GeneSetScore", 
  spatialFeatures = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  legend.position = "none",
  imageBg=FALSE,
  nrow=1
)

ggsave(paste0(outputPath,"/Visium_",sampleName,"_EMT.png"), p, width = 8.88, height = 9, dpi=300, units = "cm")




if(sampleName=="HCC-R1")
{
	newOrder <- c(3,2,1)
}else{
	newOrder <- c(1,3,2)
}

weight.W <- SpaCET_obj @results $SecAct_output $pattern $ weight.W
weight.W <- weight.W[,newOrder]
colnames(weight.W) <- as.character(1:3)
weight.W -> SpaCET_obj @results $SecAct_output $pattern $ weight.W


signal.H <- SpaCET_obj @results $SecAct_output $pattern $ signal.H
signal.H <- signal.H[newOrder,]
rownames(signal.H) <- as.character(1:3)
signal.H -> SpaCET_obj @results $SecAct_output $pattern $ signal.H



p <- SpaCET.visualize.spatialFeature(
	SpaCET_obj, 
	spatialType = "SignalingPattern", 
	spatialFeatures="All", 
	imageBg=FALSE,
	legend.position = "none",
	colors=c("#03c383","#fbbf45","#ef6a32"),
	pointSize=0.8
)
ggsave(paste0(outputPath,"/Visium_",sampleName,"_pattern_1row.png"), p, width = 25, height = 9, dpi=100, units = "cm")

write.csv(
	t(SpaCET_obj@results$SecAct_output$pattern$signal.H),
	paste0(outputPath,"/Visium_",sampleName,"_pattern_1row.csv"),quote=F)


# scatter emt and pattern 2
fg.df <- data.frame(
	x=SpaCET_obj@results$GeneSetScore["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",], 
	y=SpaCET_obj@results$SecAct_output$pattern$signal.H[2,]
)

cor_res <- cor.test(fg.df[,1],fg.df[,2])
rv <- round(cor_res$estimate,2)
pv <- signif(cor_res$p.value,2)

if(pv==0)
{
	labelText <- paste0("r = ",rv,"\np < 2.2e-16")
}else{
	labelText <- paste0("r = ",rv,"\np = ",pv)
}

library(ggplot2)
p1 <- ggplot(fg.df,aes(x=x, y=y)) + 
	geom_point(color="#ff8080", alpha=0.5, size=0.6)+
	annotate("text", x=0.17, y=0.28, label=labelText, size=4)+
	xlab("EMT gene expression")+
	ylab("Pattern 2 score")+
	theme_classic()+ 
	theme(
		panel.background = element_rect(fill='transparent'),
		plot.background = element_rect(fill='transparent', color=NA),
		axis.title = element_text(colour = "black"),
		axis.text = element_text(colour = "black")
	)
ggsave(paste0(outputPath,"/Visium_",sampleName,"_EMT_pattern2.png"), p1, width = 6.1, height = 6, dpi=300, units = "cm")
write.csv(fg.df,paste0(outputPath,"/Visium_",sampleName,"_EMT_pattern2.csv"),quote=F)



for(n in 1:3)
{
	pattern.gene <- SecAct.signaling.pattern.gene(SpaCET_obj, n)
	
	#sortrownames
	write.csv(pattern.gene, paste0(outputPath,"/Visium_",sampleName,"_pattern_",n,"_gene.csv"), quote=F)
}


corRes <- cor(t(SpaCET_obj@results$GeneSetScore), t(SpaCET_obj@results$SecAct_output$pattern$signal.H))
write.csv(corRes, paste0(outputPath,"/Visium_",sampleName,"_pattern_hallmark.csv"), quote=F)



corRes[order(corRes[,1]),]

n <- 10
SPs.1 <- rownames(corRes[order(corRes[,1],decreasing=T),])[1:n]
SPs.2 <- rownames(corRes[order(corRes[,2],decreasing=T),])[1:n]
SPs.3 <- rownames(corRes[order(corRes[,3],decreasing=T),])[1:n]

SPs <- unique(c(SPs.1, SPs.2, SPs.3))

fg.mat <- corRes[SPs,]

p <- SecAct.heatmap.plot(fg.mat)

ggsave(paste0(outputPath,"/Visium_",sampleName,"_pattern_hallmark_top.png"), p, width = 12, height = 17, dpi=200, units = "cm")



}
	
parallel::mclapply(sampleNames, runParallel, mc.cores=2) 







hallmark <- read.csv(paste0(data_MSigDB_path,"hallmark_rename.csv"),header=F,row.names=1)
hallmark[hallmark[,1]=="Epithelial mesenchymal transition",1] <- "EMT"
hallmark[hallmark[,1]=="Reactive oxygen species pathway",1] <- "Reactive oxygen species"


sampleName <- "HCC-R1"
corRes <- read.csv(paste0(outputPath,"/Visium_",sampleName,"_pattern_hallmark.csv"),row.names=1,check.names=F)
rownames(corRes) <- hallmark[rownames(corRes),1]
SPs.1 <- rownames(corRes[order(corRes[,1],decreasing=T),])[1:8]
SPs.2 <- rownames(corRes[order(corRes[,2],decreasing=T),])[1:8]
SPs.3 <- rownames(corRes[order(corRes[,3],decreasing=T),])[1:5]
SPs <- unique(c(SPs.1, SPs.2, SPs.3))


for(sampleName in sampleNames)
{
	corRes <- read.csv(paste0(outputPath,"/Visium_",sampleName,"_pattern_hallmark.csv"),row.names=1,check.names=F)
	rownames(corRes) <- hallmark[rownames(corRes),1]
	
	mat <- t(corRes[SPs,])
	
	png(paste0(outputPath,"/Visium_",sampleName,"_pattern_hallmark_show.png"), width = 14, height = 6.5, res=200, units = "cm")
	
	library(ComplexHeatmap)
	ht <- Heatmap(as.matrix(mat),
		name = "Cor",
		rect_gp = gpar(col = "white", lwd = 2),
		col = circlize::colorRamp2(c(-0.8, 0, 0.8), c("#03c383","white","#ef6a32")),
		row_names_max_width = max_text_width(
	        rownames(mat), 
	        gp = gpar(fontsize = 10)
	        ),
	    column_names_max_height = max_text_width(
	        colnames(mat), 
	        gp = gpar(fontsize = 12)
	        ),
	    column_names_gp = gpar(fontsize = 10),
	    show_row_names = FALSE,
	    show_column_names = TRUE,
		cluster_rows = FALSE,
		cluster_columns = FALSE,
		show_heatmap_legend = FALSE
		)
		
	draw(ht)
	
	dev.off()
}
