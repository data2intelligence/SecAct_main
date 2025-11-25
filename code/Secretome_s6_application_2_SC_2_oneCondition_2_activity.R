source("Secretome_s0_path.R")
library(anndata)
library(SecAct)
library(Seurat)
library(ggplot2)
library(patchwork)

outputPath <- paste0(applicationPath,"scRNAseq_PanCancer/")


#SWARM -t 2 -g 200 --time 15:00:00
args = commandArgs(trailingOnly=TRUE)
cancer <- args[1]

adata <- read_h5ad(paste0(outputPath,cancer,".h5ad"))
	
meta <- adata$obs
meta$newCluster <- "Others"

meta[meta$majorCluster=="CD4T","newCluster"] <- "CD4T"
meta[meta$majorCluster=="CD8T","newCluster"] <- "CD8T"
meta[meta$majorCluster=="Endothelial","newCluster"] <- "Endothelial"

meta[meta$subCluster=="CD4T08_Treg_FOXP3","newCluster"] <- "Treg"
meta[meta$subCluster=="CD8T07_Tex_CXCL13","newCluster"] <- "Tex"

meta[sapply(strsplit(as.character(meta$subCluster),"_",fixed=T),function(x) return(x[1]))%in%c("B01","B02","B03","B04","B05","B06","B07","B013","B014"),"newCluster"] <- "B"
meta[sapply(strsplit(as.character(meta$subCluster),"_",fixed=T),function(x) return(x[1]))%in%c("Epi"),"newCluster"] <- "Tumor"

meta[sapply(strsplit(as.character(meta$subCluster),"_",fixed=T),function(x) return(x[2]))%in%c("CD16hiNK","CD16loNK"),"newCluster"] <- "NK"
meta[sapply(strsplit(as.character(meta$subCluster),"_",fixed=T),function(x) return(x[2]))%in%c("cDC","cDC1","cDC2"),"newCluster"] <- "cDC"
meta[sapply(strsplit(as.character(meta$subCluster),"_",fixed=T),function(x) return(x[2]))%in%c("pDC"),"newCluster"] <- "pDC"
meta[sapply(strsplit(as.character(meta$subCluster),"_",fixed=T),function(x) return(x[2]))%in%c("Mph"),"newCluster"] <- "Macrophage"
meta[sapply(strsplit(as.character(meta$subCluster),"_",fixed=T),function(x) return(x[2]))%in%c("Fb","iCAF","myCAF","apCAF"),"newCluster"] <- "Fibroblast"

meta$newCluster <- factor(meta$newCluster)

Seurat_obj <- CreateSeuratObject(counts = t(as.matrix(adata$X)), meta.data = meta)
Seurat_obj <- subset(Seurat_obj, subset = newCluster != "Others")
Seurat_obj <- subset(Seurat_obj, subset = tissue == "Tumor")

Seurat_obj <- SecAct.CCC.scRNAseq(
  Seurat_obj, 
  cellType_meta="newCluster",
  condition_meta=NULL
)  

saveRDS(Seurat_obj, file = paste0(outputPath,cancer,"_single_condition_CCC_Seurat_obj.rds"))

