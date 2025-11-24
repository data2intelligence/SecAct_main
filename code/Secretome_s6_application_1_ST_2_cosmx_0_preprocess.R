source("Secretome_s0_path.R")
library(Seurat)
library(ggplot2)
library(patchwork)
library(SpaCET)
library(SecAct)

inputPath <- paste0(dataAppPath,"CosMx/")
outputPath <- paste0(applicationPath,"LIHC/")

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


# read ST data
seurat_obj_LIHC <- readRDS(paste0("../_raw/CosMx/LiverDataReleaseSeurat_newUMAP_",sampleName,".rds"))

# extract count and cell type data
counts <- seurat_obj_LIHC@assays$RNA@counts
niche_vec <- seurat_obj_LIHC@meta.data[,c("niche")]
cellType_vec <- seurat_obj_LIHC@meta.data[,c("cellType")]
coordinate_mat <- seurat_obj_LIHC@meta.data[,c("x_slide_mm","y_slide_mm")]

coordinate_mat <- coordinate_mat*1000
colnames(coordinate_mat) <- c("coordinate_x_um","coordinate_y_um")


# rename cell type
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


# create SpaCET object
cellType_mat <- data.frame(cellType=cellType_vec)
rownames(cellType_mat) <- colnames(counts)

spotCoordinates <- coordinate_mat
metaData <- cbind(cellType_mat,niche=niche_vec)
save(counts, spotCoordinates, metaData, file = paste0(inputPath,"/LIHC_CosMx_data.rda"))

