source("Secretome_s0_path.R")

inputPath <- paste0(dataAppPath, "LY86_RNAseq/")
outputPath <- paste0(applicationPath,"LY86_RNAseq/")

TPM <- read.csv(gzfile(paste0(inputPath,"GSE291028_RSEM.TPM.processed.gz")),sep="\t")
TPM <- TPM[,sort(colnames(TPM))]


expr <- data.frame()
expr[rownames(TPM),"Ly86_1"] <- rowMeans(TPM[,c("Ly86_1_1","Ly86_1_2")])
expr[rownames(TPM),"Ly86_2"] <- TPM[,c("Ly86_2_1")]
expr[rownames(TPM),"Ly86_3"] <- TPM[,c("Ly86_3_1")]
expr[rownames(TPM),"Ly86_4"] <- TPM[,c("Ly86_4")]

expr[rownames(TPM),"Vec_1"] <- rowMeans(TPM[,c("Vector_1_1","Vector_1_2")])
expr[rownames(TPM),"Vec_2"] <- rowMeans(TPM[,c("Vector_2_1","Vector_2_2")])
expr[rownames(TPM),"Vec_3"] <- TPM[,c("Vector_3_1")]

expr <- expr[rowSums(expr)>1,]



############
# QC
############

logCPMs <- expr

# Calculate rowwise variance
rv <- apply(logCPMs, 1, var)

# Sort decreasingly and take top 1000
o <- order(rv, decreasing=TRUE)
top1000 <- head(o, 1000)

# From the logCPMs subset for the top-1000
logCPM_top1000 <- logCPMs[top1000,]

# Run PCA
pca <- prcomp(t(logCPM_top1000))

# Combine PCA coordinates with the metadata from the DGEList
to_plot <- data.frame(sample=rownames(pca$x),batch=c(rep("Ly86",4),rep("Vec",3)),pca$x, Ly86_expr=unlist(expr["Ly86",]))
to_plot[["batch"]] <- factor(to_plot[["batch"]], levels=c("Vec","Ly86"))

# Calculate how many % of total variance is explained by each principal component
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100

# We focus here on PC1 and PC2
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
labs <- paste0(paste0("PC", use.pcs, " - "), paste0(round(percentVar[use.pcs], 2), "%"))

library(ggplot2)
p <- ggplot(to_plot, aes(x=PC1, y=PC2, color=batch)) + 
	geom_point(size=2) +
	#geom_text(aes(label=sample)) +
	scale_color_manual(values=c("blue","red"))+
	xlab(labs[1]) + 
	ylab(labs[2]) + 
	theme_classic()+ 
	theme(
	  panel.grid = element_blank(),
	  panel.background = element_blank(),
	  axis.text = element_text(colour = "black"),
	  axis.title = element_text(colour = "black")
	)	
ggsave(paste0(outputPath,"RNAseq_pca.pdf"), p, width =8.5, height =6, units = "cm")


pv <- wilcox.test(unlist(expr["Ly86",1:4]),unlist(expr["Ly86",5:7]), "greater")$p.value

p1 <- ggplot(to_plot,aes(x=batch, y=Ly86_expr, group=batch, color=batch))+
	geom_boxplot(color="grey",width=0.6)+
	geom_point()+
	scale_color_manual(values=c("blue","red"))+
	annotate("text", x = 2, y=5, label = paste0("p = ",signif(pv,2)))+
	ylab("Ly86 Expression")+
	theme_classic()+
	theme(
		legend.position = "right"
	)

ggsave(paste0(outputPath,"RNAseq_Ly86.pdf"), p1, width =7.4, height =5.5, units = "cm")
write.csv(to_plot,paste0(outputPath,"RNAseq_Ly86.csv"), quote=FALSE)



logCPMs <- transferMouseToHuman(as.matrix(logCPMs))
logCPMs <- rm_duplicates(logCPMs)


expr.treatment <- logCPMs[,1:4]
expr.control <- logCPMs[,5:7]

expr.diff <- rowMeans(expr.treatment) - rowMeans(expr.control)

# generate input matrix
expr.diff <- as.matrix(expr.diff, ncol=1)
colnames(expr.diff) <- "Diff"

write.table(expr.diff, paste0(outputPath,"Ly86.overexpression_vs_Vector.diff"),  quote=F)

library(SecAct)
# Run against to differential profile
res <- SecAct.activity.inference(
  inputProfile = expr.diff, 
  is.differential = TRUE
)


 
 
############
# GSEA
############

library(limma)
TT <- c(1,1,1,1,0,0,0)
WT <- c(0,0,0,0,1,1,1)
design <- cbind(TT,WT)
fit <- lmFit(expr,design)
cont.matrix <- makeContrasts(TTvsWT=TT-WT,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
res <- topTable(fit2,coef=1,number=nrow(expr))
res <- res[order(res[,"t"],decreasing=T),]

write.csv(res,paste0(outputPath,"limma_res.csv"), quote=FALSE)




rnk <- res[,3]
names(rnk) <- rownames(res)

names(rnk) <- transferSymbolFromMouseToHuman(names(rnk))

rnk <- rnk[!names(rnk)%in%c("humanGeneNotExist","humanGeneMultiple")]



library(fgsea)
library(SpaCET)

gmtNames <- list.files(data_MSigDB_path)	

for(gmtName in gmtNames)
{
	gmt <- read.gmt(paste0(data_MSigDB_path,gmtName))
	gmt_length <- lapply(gmt, length)
	gmt <- gmt[gmt_length>=15&gmt_length<=500]
	
	fgseaRes <- fgsea(pathways = gmt, stats = rnk, nperm=1000)
	fgseaRes <- fgseaRes[order(fgseaRes[,5],decreasing=T),]
	
	write.csv(as.matrix(fgseaRes[,1:7]),paste0(outputPath,gmtName,".csv"), quote=F)
}







gsNames <- c("HALLMARK_INTERFERON_GAMMA_RESPONSE","GOBP_DNA_REPLICATION")
gmtNames <- c("h.all.v2023.2.Hs.symbols.gmt","c5.go.bp.v2023.2.Hs.symbols.gmt")
leadingEdges <- list()

for(i in 1:2)
{
	set.seed(123)

	gmt <- read.gmt(paste0(data_MSigDB_path,gmtNames[i]))
	gmt_length <- lapply(gmt, length)
	gmt <- gmt[gmt_length>=15&gmt_length<=500]
	
	fgseaRes <- fgsea(pathways = gmt, stats = rnk, nperm=1000)
	fgseaRes <- fgseaRes[order(fgseaRes[,5],decreasing=T),]	
	
	NES_value <- round(fgseaRes[c(fgseaRes[,"pathway"]==gsNames[i]),"NES"],2)
	Padj_value <- signif(fgseaRes[c(fgseaRes[,"pathway"]==gsNames[i]),"padj"],2)
	leadingEdges[[gsNames[i]]] <- unlist(fgseaRes[c(fgseaRes[,"pathway"]==gsNames[i]),"leadingEdge"])
	
	p1 <- plotEnrichment(gmt[[gsNames[i]]],rnk) + labs(title="")+
	    #geom_point(color="skyblue", size=0.1)+
	    #geom_line(color="skyblue")+
	    ggtitle(paste0("NES = ",NES_value,", Padj = ",Padj_value))+
	    xlab("Ranked gene list")+
	    ylab("Enrichment score")+
	    theme(
			panel.grid.minor = element_blank(),
			panel.grid.major = element_blank(),
			axis.text = element_text(size=11,colour = "black"),
			axis.title = element_text(size=13,colour = "black"),
			axis.line.y.left = element_line(color = 'black')
	)
	
	ggsave(paste0(outputPath,"GSEA_",gsNames[i],".pdf"), p1, width = 8.2, height = 6, units = "cm")
	
}
	




p1 <- c(
	"HALLMARK_INFLAMMATORY_RESPONSE",
	"HALLMARK_INTERFERON_ALPHA_RESPONSE",
	"HALLMARK_INTERFERON_GAMMA_RESPONSE",
	"HALLMARK_TNFA_SIGNALING_VIA_NFKB",
	"HALLMARK_MTORC1_SIGNALING",
	"HALLMARK_MYC_TARGETS_V1",
	"HALLMARK_E2F_TARGETS",
	"HALLMARK_G2M_CHECKPOINT")
p2 <- c(
	"GOBP_MYELOID_LEUKOCYTE_ACTIVATION",
	"GOBP_REGULATION_OF_CELL_KILLING",
	"GOBP_TUMOR_NECROSIS_FACTOR_SUPERFAMILY_CYTOKINE_PRODUCTION",
	"GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION",
	"GOBP_RIBOSOME_BIOGENESIS",
	"GOBP_DNA_REPLICATION",
	"GOBP_RNA_SPLICING",
	"GOBP_REGULATION_OF_TRANSLATIONAL_FIDELITY")
	
p3 <- c(
	"KEGG_COMPLEMENT_AND_COAGULATION_CASCADES",
	"KEGG_ARACHIDONIC_ACID_METABOLISM",
	"KEGG_CELL_CYCLE",
	"KEGG_MISMATCH_REPAIR"
)

fg.df <- data.frame()

fgseaRes1 <- read.csv(paste0(outputPath,"h.all.v2023.2.Hs.symbols.gmt.csv"),row.names=2)
fgseaRes2 <- read.csv(paste0(outputPath,"c5.go.bp.v2023.2.Hs.symbols.gmt.csv"),row.names=2)
fgseaRes3 <- read.csv(paste0(outputPath,"c2.cp.kegg_legacy.v2023.2.Hs.symbols.gmt.csv"),row.names=2)

fg.df[p1,"NES"] <- fgseaRes1[p1,"NES"]
fg.df[p2,"NES"] <- fgseaRes2[p2,"NES"]
fg.df[p3,"NES"] <- fgseaRes3[p3,"NES"]


rownames(fg.df) <- gsub("HALLMARK_","",rownames(fg.df)) 
rownames(fg.df) <- gsub("GOBP_","",rownames(fg.df)) 
rownames(fg.df) <- gsub("KEGG_","",rownames(fg.df)) 


rownames(fg.df) <- gsub("_"," ",rownames(fg.df)) 
rownames(fg.df) <- tolower(rownames(fg.df)) 

rownames(fg.df) <- stringr::str_to_title(rownames(fg.df)) 
rownames(fg.df) <- gsub("Tnfa","TNFA",rownames(fg.df)) 
rownames(fg.df) <- gsub("Nfkb","NF-kB",rownames(fg.df)) 
rownames(fg.df) <- gsub("E2f","E2F",rownames(fg.df)) 
rownames(fg.df) <- gsub("Myc","MYC",rownames(fg.df)) 
rownames(fg.df) <- gsub("G2m","G2M",rownames(fg.df)) 
rownames(fg.df) <- gsub("Rna","RNA",rownames(fg.df)) 
rownames(fg.df) <- gsub("Dna","DNA",rownames(fg.df)) 
rownames(fg.df) <- gsub("Of","of",rownames(fg.df)) 
rownames(fg.df) <- gsub("And","&",rownames(fg.df)) 
rownames(fg.df) <- gsub("Tumor Necrosis Factor","TNF",rownames(fg.df)) 
rownames(fg.df) <- gsub("Mtorc1","mTORC1",rownames(fg.df)) 



fg.df <- fg.df[order(fg.df[,"NES"]),,drop=F]
fg.vec <- fg.df[,"NES"]

fg.df <- cbind(pw=rownames(fg.df), fg.df)

  fg.df <- cbind(fg.df, dir=ifelse(fg.vec<0,"down","up"))
  fg.df <- cbind(fg.df, y=ifelse(fg.vec<0,0.1,-0.1))
  fg.df <- cbind(fg.df, hjust=ifelse(fg.vec<0,0,1))
  fg.df[["pw"]] <- factor(fg.df[["pw"]], levels=rownames(fg.df))

library(ggplot2)
p <- ggplot(fg.df, aes(pw, NES, label=pw)) +
    geom_col(aes(fill=dir), width = .88, color = "white", alpha=0.6) +
    geom_text(aes(y = y, hjust=hjust), angle = 0, size = 2.5) +
    scale_fill_manual(values=c("#66bd63","#f46d43"))+
    geom_hline(yintercept=0)+
    ggtitle("NES (Normalized Enrichment Score)")+
    theme_classic()+
    theme(
      plot.background = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5),
      axis.title = element_blank(),
      axis.text.x = element_text(color="black", vjust=0.5),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.line.y = element_blank(),
      legend.position = "none"
    )+
    coord_flip()+
    scale_y_continuous(position = "right",limits=c(-3.5,3.5))


ggsave(paste0(outputPath,"GSEA_summary.pdf"), p, width = 9.2, height = 12, units = "cm")
write.csv(fg.df,paste0(outputPath,"GSEA_summary.csv"), quote=F)





fg.df <- res[,c("logFC","P.Value")]
fg.df <- fg.df[order(fg.df[,1]),]
fg.df[,"neg_log10_pval"] <- -log10(fg.df[,"P.Value"])


# Define significance threshold
pval_cutoff <- 0.05
log2fc_cutoff <- 0.5

# Label significant genes
fg.df$significance <- ifelse(fg.df$P.Value < pval_cutoff & abs(fg.df$logFC) > log2fc_cutoff, 
                            ifelse(fg.df$logFC > 0, "Upregulated", "Downregulated"), 
                            "Not Significant")

fg.df$gene <- rownames(fg.df)
fg.df[abs(fg.df$logFC) * fg.df$neg_log10_pval <5.5,"gene"] <- ""


fg.df$human <- transferSymbolFromMouseToHuman(rownames(fg.df))

get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x, y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  density <- dens$z[cbind(ix, iy)]
  return(density)
}

# Add density to the data
fg.df$density <- get_density(fg.df$logFC, fg.df$neg_log10_pval)

library(ggplot2)
library(ggrepel)

p <- ggplot(fg.df, aes(x = logFC, y = neg_log10_pval, color=density, label=gene)) +
  geom_point(alpha = 0.3, size=0.1) +  # Scatter points
  geom_text_repel(size=2,max.overlaps=20,color="black")+
  #scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  #geom_vline(xintercept = 0, linetype = "dashed") +  # Add vertical lines for cutoff
  #geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed") +  # Add vertical lines for cutoff
  #geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed") +  # Add horizontal line for p-value threshold
  labs(title = " ", 
       x = "Log2 Fold Change",
       y = "-log10(P-value)")+
  theme_classic()+
  theme(
  	plot.title = element_text(hjust = 0.5),
  	legend.position = "None"
  )

ggsave(paste0(outputPath,"volcano.png"), p, width = 8.2, height = 6, dpi=300, units = "cm")




















if(FALSE)
{
#################
# deconvolution #
#################


library(Seurat)
library(ggplot2)


count <- read.csv(gzfile(paste0(dataPath,"GSE228014_table_processed.txt.gz")),sep="\t",row.names=1)
TPM <- exp(count)-1
TPM <- TPM*100
TPM <- round(TPM)

sc <- CreateSeuratObject(counts = TPM)
sc <- NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = 10000)
sc <- FindVariableFeatures(sc, nfeatures = 4000)
sc <- ScaleData(sc)

sc <- RunPCA(sc, npcs = 30, verbose = FALSE)
sc <- RunUMAP(sc, dims = 1:30)
sc <- FindNeighbors(sc, dims = 1:30)
sc <- FindClusters(sc, resolution = 2.5) #V1

	g <- DimPlot(sc, reduction = reduc, group.by = "seurat_clusters",label = TRUE)
	#ggsave(paste0(CombPath,"doubletsRemove_",doubletsRemove,"_","All_",reduc,"_group.by_clustering.jpg"), g, width = 18, height = 16, dpi=300, units = "cm")
	ggsave(paste0(dataPath,"All_",reduc,"_group.by_clustering.pdf"), g, width = 18, height = 16, dpi=300, units = "cm")
	


Idents(sc) <-  "seurat_clusters"
new.cluster.ids <- c("c1", "c1", "c1", "c1", "c1", "c1", "c1", "c711", "c8", "c9", "c10", "c711", "c12", "c13")
names(new.cluster.ids) <- levels(sc)
sc <- RenameIdents(sc, new.cluster.ids)

	g <- DimPlot(sc, reduction = reduc, label = TRUE) + NoLegend()
	#ggsave(paste0(CombPath,"doubletsRemove_",doubletsRemove,"_","All_",reduc,"_group.by_clustering.jpg"), g, width = 18, height = 16, dpi=300, units = "cm")
	ggsave(paste0(dataPath,"All_",reduc,"_group.by_clustering_new.pdf"), g, width = 18, height = 16, dpi=300, units = "cm")
	


# https://www.sciencedirect.com/science/article/pii/S1074761319301268?via%3Dihub
# https://www.cellsignal.com/pathways/immune-cell-markers-mouse
genes <- c(
"Jchain", # plasma
"Cd19", # B
"Nkx1-1","Ncr1","Klrk1", # NK
"Cd3e",
"Mki67", # proliferating
"Ccr7","Tcf7","Lef1","Sell", # naive
"Nkg7","Ifng","Gzma","Gzmb","Prf1","Ccl4", # cytotoxic
"Pdcd1","Tigit","Lag3","Havcr2", # (Tim3 - Havcr2)
"Cd4","Foxp3","Cd8a","Cd8b1",
"Itgam",# Myeloid(CD11b)
"Csf1r", # MonoMacroDC
"Cd68", # MonoMacro
"Cd14", # Monocyte
"Adgre1", # Macrophage(F4/80)
"Cd86","Cd80", # M1
"Cd163","Mrc1", # M2 (Cd206 - Mrc1)
"Itgax",# cDC(CD11c)
"Xcr1","Clec9a",# cDC1
"Sirpa",# cDC2
"Siglech","Bst2", # pDC (Cd317 - Bst2)
"Arg1","Arg2", # MDSC
"Ly6c1", # M-MDSC
"Ly6g", # PMN-MDSC
"Csf3r", # Neutrophil
"Ccr3","Siglecf", # Eosinophil
"Fcer1a", # Basophil
"Kit","Fcer2a", # Mast cell (Cd117 - Kit, Cd23 - Fcer2a)
"Nos2","Tgfb1","Tgfb2","Tgfb3",
"Tyrp1","Pecam1","Col1a1"
)

reduc <- "umap"
	
	g <- FeaturePlot(sc, reduction = reduc, features = genes, ncol =8)
	#ggsave(paste0(CombPath,"doubletsRemove_",doubletsRemove,"_","All_",reduc,"_marker.jpg"), g, width = 18*8, height = 16*6, dpi=300, units = "cm", limitsize = FALSE)
	ggsave(paste0(dataPath,"All_",reduc,"_marker.pdf"), g, width = 18*8, height = 16*6, dpi=300, units = "cm", limitsize = FALSE)
	



sc_counts <- as.matrix(sc@assays$ RNA@ layers$ counts  )
rownames(sc_counts) <- rownames(TPM)
colnames(sc_counts) <- colnames(TPM)

sc_annotation <- data.frame(cellID=colnames(sc_counts), bio_celltype=as.character(sc@ active.ident) )

sc_lineageTree <- as.list(unique(new.cluster.ids))
names(sc_lineageTree) <- unique(new.cluster.ids)


TPM <- read.csv(gzfile("/data/rub2/project/Secretome/results/LY86_RNAseq/RSEM.TPM.processed.gz"),sep="\t")
TPM <- TPM[,sort(colnames(TPM))]
TPM <- 2^TPM-1

expr <- data.frame()
expr[rownames(TPM),"Ly86_1"] <- rowMeans(TPM[,c("Ly86_1_1","Ly86_1_2")])
expr[rownames(TPM),"Ly86_2"] <- TPM[,c("Ly86_2_1")]
expr[rownames(TPM),"Ly86_3"] <- TPM[,c("Ly86_3_1")]
expr[rownames(TPM),"Ly86_4"] <- TPM[,c("Ly86_4")]

expr[rownames(TPM),"Vec_1"] <- rowMeans(TPM[,c("Vector_1_1","Vector_1_2")])
expr[rownames(TPM),"Vec_2"] <- rowMeans(TPM[,c("Vector_2_1","Vector_2_2")])
expr[rownames(TPM),"Vec_3"] <- TPM[,c("Vector_3_1")]

st_counts <- as.matrix(expr)
spotCoordinates <- as.matrix(cbind(X=1:7,Y=1:7))
rownames(spotCoordinates) <- colnames(st_counts)

library(SpaCET)
SpaCET_obj <- create.SpaCET.object(
  counts=st_counts,
  spotCoordinates=spotCoordinates,
  imagePath=NA,
  platform = "bulk"
)

SpaCET_obj <- SpaCET.deconvolution.matched.scRNAseq(
  SpaCET_obj, 
  sc_counts=sc_counts, 
  sc_annotation=sc_annotation, 
  sc_lineageTree=sc_lineageTree, 
  coreNo=1
)


}