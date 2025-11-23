source("Secretome_s0_path.R")
source("../../SpaCE/code/ifun.R")
library(Seurat)
library(ggplot2)
library(patchwork)
library(SpaCET)
library(SecAct)

objPath <- "../data/"


cancer <- "LIHC"
sampleName <- "CancerousLiver"

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


seurat_obj_LIHC <- readRDS(paste0(objPath,"CosMx/LiverDataReleaseSeurat_newUMAP_",sampleName,".rds"))

counts <- seurat_obj_LIHC@assays$RNA@counts
niche_vec <- seurat_obj_LIHC@meta.data[,c("niche")]
cellType_vec <- seurat_obj_LIHC@meta.data[,c("cellType")]
coordinate_mat <- seurat_obj_LIHC@meta.data[,c("x_slide_mm","y_slide_mm")]


cellType_vec[grepl("Antibody.secreting.B.cells",cellType_vec)] <- "B"
cellType_vec[grepl("Mature.B.cells",cellType_vec)] <- "B"

cellType_vec[grepl("CD3+.alpha.beta.T.cells",cellType_vec,fixed=T)] <- "T.alpha.beta"
cellType_vec[grepl("gamma.delta.T.cells.1",cellType_vec)] <- "T.gamma.delta"
cellType_vec[grepl("NK.like.cells",cellType_vec)] <- "NK"

cellType_vec[grepl("Hep",cellType_vec)] <- "Hepatocyte"
cellType_vec[grepl("Cholangiocytes",cellType_vec)] <- "Cholangiocyte"
cellType_vec[grepl("Erthyroid.cells",cellType_vec)] <- "Erythrocyte"

cellType_vec[grepl("Inflammatory.macrophages",cellType_vec)] <- "Macrophage"
cellType_vec[grepl("Non.inflammatory.macrophages",cellType_vec)] <- "Macrophage"

cellType_vec[grepl("Central.venous.LSECs",cellType_vec)] <- "Endothelial"
cellType_vec[grepl("Periportal.LSECs",cellType_vec)] <- "Endothelial"
cellType_vec[grepl("Portal.endothelial.cells",cellType_vec)] <- "Endothelial"
cellType_vec[grepl("Stellate.cells",cellType_vec)] <- "Fibroblast"

cellType_vec[niche_vec=="interface"] <- "Tumor_boundary"	
cellType_vec[grepl("tumor_1",cellType_vec)] <- "Tumor_core"
cellType_vec[grepl("tumor_2",cellType_vec)] <- "Tumor_core"




cellType_mat <- data.frame(cellType=cellType_vec)
rownames(cellType_mat) <- colnames(counts)

spotCoordinates <- coordinate_mat
metaData <- cbind(cellType_mat,niche=niche_vec)
save(counts, spotCoordinates, metaData, file = paste0(applicationPath,cancer,"/LIHC_CosMx_data.rda"))


SpaCET_obj <- create.SpaCET.object(
  counts=counts,
  spotCoordinates=coordinate_mat,
  metaData=cellType_mat,
  imagePath=NA,
  platform = "CosMx"
)


SpaCET_obj <- SpaCET.quality.control(SpaCET_obj, min.genes=50)


#
#
#p <- SpaCET.visualize.spatialFeature(
#  SpaCET_obj, 
#  spatialType = "metaData", 
#  spatialFeatures= "cellType",
#  colors = my_cols,
#  legend.position="none",
#  pointSize = 0.6
#)
#
#
#ggsave(paste0(applicationPath,cancer,"/CosMx_",sampleName,"_cellType.png"), p, width = 13, height = 13, dpi=500, units = "cm")
#
#
#
ref <- paste0("/data/rub2/project/Secretome/results/signatureComb/AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3_1k_vst_condition_logUMI_cellType_0.9.tsv")

Sys.time()

SpaCET_obj <- SecAct.activity.inference.ST(
    inputProfile=SpaCET_obj,
    inputProfile_control = NULL,
    scale.factor = 1000,
    sigMatrix=ref,
    lambda=5e+5,
    nrand=1000,
    sigFilter=TRUE
)

Sys.time()



SpaCET_obj <- readRDS(paste0(applicationPath,cancer,"/CosMx_",sampleName,"_SpaCET_obj_618.rds"))

SpaCET_obj <- SecAct.CCC.scST(
    SpaCET_obj,
    cellType_meta = "cellType",
    scale.factor = 1000,
    radius = 0.02,
    ratio_cutoff = 0.2,
    padj_cutoff = 0.01
)

saveRDS(SpaCET_obj, file = paste0(applicationPath,cancer,"/CosMx_",sampleName,"_SpaCET_obj.rds"))




png(paste0(applicationPath,cancer,"/CosMx_",sampleName,"_heatmap.png"), width = 15, height = 15, res=500, units = "cm")

SecAct.CCC.heatmap(SpaCET_obj, row.sorted=TRUE, column.sorted=TRUE, colors_cellType=my_cols)

dev.off()


png(paste0(applicationPath,cancer,"/CosMx_",sampleName,"_circle.png"), width = 13, height = 13, res=300, units = "cm")

SecAct.CCC.circle(SpaCET_obj, colors_cellType=my_cols)

dev.off()




ccc <- SpaCET_obj @results $SecAct_output $SecretedProteinCCC


write.csv(ccc, paste0(applicationPath,cancer,"/SecretedProteinCCC.csv"), quote=F)






cellTypes <- c("Tumor_boundary","Fibroblast","Macrophage","Endothelial")
ccc_boundary <- ccc[ccc[,1]%in%cellTypes&ccc[,3]%in%cellTypes,]

cosmx <- sort(unique(ccc_boundary[,2]))

visium1 <- read.csv(paste0(applicationPath,cancer,"/Visium_HCC-R1_pattern_2_gene.csv"),row.names=1)
visium2 <- read.csv(paste0(applicationPath,cancer,"/Visium_HCC-R2_pattern_2_gene.csv"),row.names=1)

visium <- intersect(rownames(visium1), rownames(visium2))
visium <- sort(unique(c(rownames(visium1), rownames(visium2))))

length(visium)
length(cosmx)

visium_specific <- setdiff(visium, cosmx)
cosmx_specific <- setdiff(cosmx, visium)
olp <- intersect(visium,cosmx)

length(visium_specific)
length(cosmx_specific)
length(olp)

phyper(length(olp),length(visium),1155-length(visium),length(cosmx),FALSE)

writeLines(visium_specific, paste0(applicationPath,cancer,"/visium_specific.csv"))
writeLines(cosmx_specific, paste0(applicationPath,cancer,"/cosmx_specific.csv"))
writeLines(olp, paste0(applicationPath,cancer,"/olp.csv"))




gmt <- read.gmt(paste0("/data/rub2/data/MSigDB/h.all.v2025.1.Hs.symbols.gmt"))
EMT <- gmt[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]]
EMT <- transferSymbol(EMT)






sort(unique(ccc_boundary[,2])) -> bb
intersect(bb,EMT)

secretedProtein <- c("BGN","COL1A1","COL1A2","DCN","IGFBP5","LGALS1","LGALS9","LYZ","LUM","MGP","SPP1","THBS1","THBS2")


p <- SecAct.CCC.dot(SpaCET_obj, colors_cellType, sender=cellTypes, secretedProtein=secretedProtein, receiver=cellTypes)
ggsave(paste0(applicationPath,cancer,"/CosMx_dot.png"), p, width = 14, height = 14, dpi=300, units = "cm")



sp1 <- ccc[ccc[,1]=="Fibroblast"&ccc[,3]=="Tumor_interface",2]
sp2 <- ccc[ccc[,1]=="Fibroblast"&ccc[,3]=="Tumor_core",2]


