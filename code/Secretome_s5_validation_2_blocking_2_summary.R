source("Secretome_s0_path.R")

inputPath <- paste0(dataValPath,"blocking/")
outputPath <- paste0(validationPath,"blocking/")


items <- c(
"TGFB_GSE174686","TGFB3_GSE174686",
"IL1B_GSE80060","IL1B_GSE57253","IL1B_GSE70019",
"IFNG_GSE100093","IFNG_GSE78193",
"TNFSF12_GSE42049","TNFSF12_GSE42048", # tumor
"IL6_GSE62941","IL6_GSE45867","IL6_GSE61201", # tumor x x 
"TNF_GSE48498","TNF_GSE11903",
"IL1A_GSE57253","IL1A_GSE70019",
"NTN1_GSE225691", # tumor
"IL11_Widjaja2024"
#"IL17A_GSE31652","IL17A_GSE11903","IL17A_GSE55201", # SecAct does not have such signatures.
#"IL22_GSE99802",
#"IL4_GSE130588","IL4_GSE59294",
#"IL13_GSE130588","IL13_GSE59294"
)


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


# TNFSF12_GSE42048 miss ligand expression
TNFSF12_GSE42048 <- read.csv(paste0(inputPath,"GSE42048.HG-U133_Plus_2.rma"),sep="\t",row.names=1)
TNFSF12_GSE42049 <- read.csv(paste0(inputPath,"GSE42049.HG-U133_Plus_2.rma"),sep="\t",row.names=1)

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


meta <- read.csv(paste0(inputPath,"cytokine_disease_GEOID.csv"),as.is=T,row.names=1)
for(i in 1:nrow(meta))
{
	if(meta[i,"Organism"]=="Mouse")
	(
		meta[i,"Dataset_alt"] <- paste0(meta[i,"Dataset_alt"], "^")
	)
}
meta[meta[,"Dataset_alt"]=="TGFB1&2&3_Breast cancer^","Dataset_alt"] <- "TGFB1&2&3_Breast cancer^#"
meta_sorted <- meta[rownames(mat),c(3,2,1,4)]
write.csv(meta_sorted,paste0(inputPath,"cytokine_disease_GEOID_sorted.csv"))

rownames(mat) <- meta[rownames(mat),"Dataset_alt"]
rownames(mat) <- paste0("Anti-",rownames(mat))


write.csv(mat,paste0(outputPath,"validation_blocking.csv"),quote=F)

pdf(paste0(outputPath,"validation_blocking.pdf"), width = 8.1, height = 7.8)

library(ComplexHeatmap)
disease_vec <- meta_sorted[,"Disease"]
dataset_vec <- rownames(mat)
	
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
