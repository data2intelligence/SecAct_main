source("Secretome_s0_path.R")
library(Seurat)
library(ggplot2)
library(patchwork)
library(SpaCET)
library(SecAct)

inputPath <- paste0(dataAppPath,"CosMx/")
outputPath <- paste0(applicationPath,"LIHC/")
sampleName <- "CancerousLiver"

load(paste0(inputPath,"/LIHC_CosMx_data.rda"))

SpaCET_obj <- create.SpaCET.object(
  counts=counts,
  spotCoordinates=spotCoordinates,
  metaData=metaData,
  imagePath=NA,
  platform = "CosMx"
)

SpaCET_obj <- SpaCET.quality.control(SpaCET_obj, min.genes=50)


my_cols <- c(
	'B'='#C88888',
	'Erythrocyte'='#fe666d',
	'T.alpha.beta'='#B95FBB',
	'T.gamma.delta'='#3288bd',
	'NK'='#bb8761',
	
	'Hepatocyte'='#63636d',
	'Cholangiocyte'='#de77ae',
	
	'Endothelial'='#D4D915',
	'Fibroblast'='#66c2a5',
	'Macrophage'='#ff9a36',
	
	'Tumor_core'='#A4DFF2',
	'Tumor_boundary'='blue'
)


p <- SpaCET.visualize.spatialFeature(
  SpaCET_obj, 
  spatialType = "metaData", 
  spatialFeatures= "cellType",
  colors = my_cols,
  legend.position="none",
  pointSize = 0.6
)
ggsave(paste0(outputPath,"/CosMx_",sampleName,"_cellType.png"), p, width = 13, height = 13, dpi=500, units = "cm")



Sys.time()

SpaCET_obj <- SecAct.activity.inference.ST(
    inputProfile=SpaCET_obj,
    scale.factor = 1000,
    sigFilter=TRUE
)

Sys.time()



SpaCET_obj <- SecAct.CCC.scST(
    SpaCET_obj,
    cellType_meta = "cellType",
    scale.factor = 1000,
    radius = 20,
    ratio_cutoff = 0.2,
    padj_cutoff = 0.01
)
saveRDS(SpaCET_obj, file = paste0(outputPath,"/CosMx_",sampleName,"_SpaCET_obj.rds"))









SpaCET_obj <- readRDS(paste0(outputPath,"/CosMx_",sampleName,"_SpaCET_obj.rds"))


png(paste0(outputPath,"/CosMx_",sampleName,"_heatmap.png"), width = 15, height = 15, res=500, units = "cm")

SecAct.CCC.heatmap(SpaCET_obj, row.sorted=TRUE, column.sorted=TRUE, colors_cellType=my_cols)

dev.off()


png(paste0(outputPath,"/CosMx_",sampleName,"_circle.png"), width = 13, height = 13, res=300, units = "cm")

SecAct.CCC.circle(SpaCET_obj, colors_cellType=my_cols)

dev.off()








ccc <- SpaCET_obj @results $SecAct_output $SecretedProteinCCC

write.csv(ccc, paste0(outputPath,"/CosMx_",sampleName,"_SecretedProteinCCC.csv"), quote=F)




# overlap between visium and cosmx
cellTypes <- c("Tumor_boundary","Fibroblast","Macrophage","Endothelial")
ccc_boundary <- ccc[ccc[,1]%in%cellTypes&ccc[,3]%in%cellTypes,]

cosmx <- sort(unique(ccc_boundary[,2]))

visium1 <- read.csv(paste0(outputPath,"/Visium_HCC-R1_pattern_2_gene.csv"),row.names=1)
visium2 <- read.csv(paste0(outputPath,"/Visium_HCC-R2_pattern_2_gene.csv"),row.names=1)

visium <- intersect(rownames(visium1), rownames(visium2))
visium <- sort(unique(c(rownames(visium1), rownames(visium2))))

length(visium)
length(cosmx)

visium_specific <- setdiff(visium, cosmx)
cosmx_specific <- setdiff(cosmx, visium)
olp <- intersect(visium,cosmx)


writeLines(visium_specific, paste0(outputPath,"/visium_specific.csv"))
writeLines(cosmx_specific, paste0(outputPath,"/cosmx_specific.csv"))
writeLines(olp, paste0(outputPath,"/olp.csv"))



Xfile<- file.path(system.file(package = "SecAct"), "extdata/SecAct.tsv.gz")
X <- read.table(Xfile,sep="\t",check.names=F)
no0 <- ncol(X)

no1 <- length(visium_specific)
no2 <- length(olp)
no3 <- length(cosmx_specific)

p <- phyper(length(olp),length(visium),no0-length(visium),length(cosmx),FALSE)

hyperTest <- c(
	paste0("visium_specific ", no1), 
	paste0("olp ", no2), 
	paste0("cosmx_specific ", no3),  
	paste0("p.value ", p)  
)
hyperTest

writeLines(hyperTest, paste0(outputPath,"/hyperTest_visium_cosmx.txt"))



# download from https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/h.all.v2025.1.Hs.symbols.gmt
EMT <- c("ABI3BP","ACTA2","ADAM12","ANPEP","APLP1","AREG","BASP1","BDNF","BGN","BMP1","CADM1","CALD1","CALU","CAP2","CAPG","CCN1","CCN2","CD44","CD59","CDH11","CDH2","CDH6","COL11A1","COL12A1","COL16A1","COL1A1","COL1A2","COL3A1","COL4A1","COL4A2","COL5A1","COL5A2","COL5A3","COL6A2","COL6A3","COL7A1","COL8A2","COLGALT1","COMP","COPA","CRLF1","CTHRC1","CXCL1","CXCL12","CXCL6","CXCL8","DAB2","DCN","DKK1","DPYSL3","DST","ECM1","ECM2","EDIL3","EFEMP2","ELN","EMP3","ENO2","FAP","FAS","FBLN1","FBLN2","FBLN5","FBN1","FBN2","FERMT2","FGF2","FLNA","FMOD","FN1","FOXC2","FSTL1","FSTL3","FUCA1","FZD8","GADD45A","GADD45B","GAS1","GEM","GJA1","GLIPR1","GPC1","GPX7","GREM1","HTRA1","ID2","IGFBP2","IGFBP3","IGFBP4","IL15","IL32","IL6","INHBA","ITGA2","ITGA5","ITGAV","ITGB1","ITGB3","ITGB5","JUN","LAMA1","LAMA2","LAMA3","LAMC1","LAMC2","LGALS1","LOX","LOXL1","LOXL2","LRP1","LRRC15","LUM","MAGEE1","MATN2","MATN3","MCM7","MEST","MFAP5","MGP","MMP1","MMP14","MMP2","MMP3","MSX1","MXRA5","MYL9","MYLK","NID2","NNMT","NOTCH2","NT5E","NTM","OXTR","P3H1","PCOLCE","PCOLCE2","PDGFRB","PDLIM4","PFN2","PLAUR","PLOD1","PLOD2","PLOD3","PMEPA1","PMP22","POSTN","PPIB","PRRX1","PRSS2","PTHLH","PTX3","PVR","QSOX1","RGS4","RHOB","SAT1","SCG2","SDC1","SDC4","SERPINE1","SERPINE2","SERPINH1","SFRP1","SFRP4","SGCB","SGCD","SGCG","SLC6A8","SLIT2","SLIT3","SNAI2","SNTB1","SPARC","SPOCK1","SPP1","TAGLN","TFPI2","TGFB1","TGFBI","TGFBR3","TGM2","THBS1","THBS2","THY1","TIMP1","TIMP3","TNC","TNFAIP3","TNFRSF11B","TNFRSF12A","TPM1","TPM2","TPM4","VCAM1","VCAN","VEGFA","VEGFC","VIM","WIPF1","WNT5A")
EMT <- transferSymbol(EMT)

sort(unique(ccc_boundary[,2])) -> bb
intersect(bb,EMT)




secretedProtein <- c("BGN","COL1A1","COL1A2","DCN","IGFBP5","LGALS1","LGALS9","LYZ","LUM","MGP","SPP1","THBS1","THBS2")

p <- SecAct.CCC.dot(SpaCET_obj, sender=cellTypes, secretedProtein=secretedProtein, receiver=cellTypes)
ggsave(paste0(outputPath,"/CosMx_",sampleName,"_dot.png"), p, width = 14, height = 14, dpi=300, units = "cm")

write.csv(
	ccc[ccc[,2]%in%secretedProtein,],
	paste0(outputPath,"/CosMx_",sampleName,"_dot.csv"), quote=FALSE)



p <- SecAct.signaling.velocity.scST(
  SpaCET_obj, 
  sender="Fibroblast", 
  secretedProtein="THBS2", 
  receiver="Tumor_boundary", 
  cellType_meta="cellType"
)
ggsave(paste0(outputPath,"/CosMx_",sampleName,"_velocity.png"), p, width = 16, height = 16, dpi=300, units = "cm")



# background
# n_cellType <- length(unique(SpaCET_obj@input$metaData[,1]))
# n_secretedProtein <- nrow(SpaCET_obj@results$SecAct_output$ SecretedProteinActivity$ zscore)
# no0 <- n_cellType * (n_cellType-1) * n_secretedProtein

spacia <- read.csv("/data/Jiang_Lab/Data/Seongyong/spacia/aggregated_results_long.csv")
rownames(spacia) <- paste0(spacia[,1],"_",spacia[,3],"_",spacia[,2])
dim(spacia)
spacia[,"protein"] <- transferSymbol(spacia[,"protein"])

no0 <- nrow(spacia)

spacia <- spacia[spacia[,"fdr"]<0.01&spacia[,"beta"]>0,]
dim(spacia)

secact_specific <- setdiff(rownames(ccc), rownames(spacia))
spacia_specific <- setdiff(rownames(spacia), rownames(ccc))
olp <- intersect(rownames(ccc), rownames(spacia))
length(olp)

no1 <- length(secact_specific)
no2 <- length(olp)
no3 <- length(spacia_specific)

p <- phyper(length(olp),nrow(ccc),no0-nrow(ccc),nrow(spacia),FALSE)

hyperTest <- c(
	paste0("secact_specific ", no1), 
	paste0("olp ", no2), 
	paste0("spacia_specific ", no3),  
	paste0("p.value ", p)  
)
hyperTest

writeLines(hyperTest, paste0(outputPath,"/hyperTest_secact_spacia.txt"))


Fibroblast_LGALS1_Macrophage
[104] "Tumor_core_CCL5_Macrophage"           
