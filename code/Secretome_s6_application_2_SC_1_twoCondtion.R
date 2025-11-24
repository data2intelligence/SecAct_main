source("Secretome_s0_path.R")

library(Seurat)
library(ggplot2)
library(patchwork)
library(SecAct)


inputPath <- paste0(dataAppPath,"OV/")
outputPath <- paste0(applicationPath,"OV/")

Seurat_obj <- readRDS(paste0(inputPath, "OV_scRNAseq_Seurat.rds"))

my_cols <- c(
	'B'='#C88888',
	'Tex'='#fe666d',
	'CD8T'='#B95FBB',
	'CD4T'='#3288bd',
	'Treg'='#E6C122',
	'NK'='#bb8761',
	'DC'='#63636d',
	'Endothelial'='#D4D915',
	'Fibroblast'='#66c2a5',
	'Macrophage'='#ff9a36',
	'Tumor'='#A4DFF2'
)

# draw CD4 first, otherwise covered by CD8
cellID_cellType <- data.frame(
	cellID = Cells(Seurat_obj),
	cellType = Seurat_obj@meta.data$MyCellType
)
cellID_cellType <- cellID_cellType[order(match(cellID_cellType[,"cellType"],names(my_cols))),]

write.csv(cellID_cellType,paste0(outputPath,"/cellType.csv"),quote=F)

p1 <- FeaturePlot(Seurat_obj, reduction = "umap", features = c("CLEC9A", "CD1C", "LAMP3", "LGALS2", "LILRA4"), ncol=5)
ggsave(paste0(outputPath,"OV_marker.png"), p1, width = 60, height = 12, dpi=500, units = "cm")


p1 <- DimPlot(Seurat_obj, reduction = "umap", cells = cellID_cellType[,1], cols = my_cols, group.by = "MyCellType")+NoLegend()+ggtitle(' \n \n ')+NoAxes()
p2 <- DimPlot(Seurat_obj, reduction = "umap", cells = sample(Cells(Seurat_obj)), group.by = "Groups")+ NoLegend()+ggtitle(' \n \n ')+NoAxes()
ggsave(paste0(outputPath,"OV_UMAP.png"), p1/p2, width = 12, height = 25, dpi=500, units = "cm")


p1 <- DimPlot(Seurat_obj, reduction = "umap", group.by = "MyCellType", split.by="Groups", cols = my_cols)+NoLegend()+ggtitle(NULL)+NoAxes()
ggsave(paste0(outputPath,"OV_UMAP2.png"), p1, width = 24, height = 12, dpi=500, units = "cm")


p1 <- DimPlot(Seurat_obj, reduction = "umap", group.by = "Patients")
ggsave(paste0(outputPath,"OV_UMAP3.png"), p1, width = 14, height = 12, dpi=500, units = "cm")



Patients <- c("HGSOC1","HGSOC3","HGSOC4","HGSOC6")

cell_stat <- data.frame()
for(Patient in Patients)
{
	Seurat_obj_patient <- subset(Seurat_obj, Patients %in% Patient )
	
	Seurat_obj_patient_Primary <- subset(Seurat_obj_patient, Groups %in% c("Primary") )
	Seurat_obj_patient_Metastatic <- subset(Seurat_obj_patient, Groups %in% c("Metastatic") )
	
	stat1 <- table(Seurat_obj_patient_Primary@meta.data[,"MyCellType"])
	stat2 <- table(Seurat_obj_patient_Metastatic@meta.data[,"MyCellType"])
	
	cell_stat[names(stat1),paste0(Patient,"_Primary")] <- stat1
	cell_stat[names(stat2),paste0(Patient,"_Metastatic")] <- stat2
}
write.csv(cell_stat, paste0(outputPath,"OV_cell_stat.csv"), quote=FALSE)


Seurat_obj <- SecAct.CCC.scRNAseq(
  Seurat_obj, 
  cellType_meta="MyCellType",
  condition_meta="Groups", 
  conditionCase="Metastatic", 
  conditionControl="Primary",
  act_diff_cutoff=2,
  exp_logFC_cutoff=0.2,
  exp_mean_all_cutoff=2,
  exp_fraction_case_cutoff=0.1,
  padj_cutoff=0.01
)  

saveRDS(Seurat_obj, file = paste0(outputPath,"/OV_Metastatic_Seurat_obj.rds"))


png(paste0(outputPath,"/OV_heatmap_Metastatic.png"), width = 12, height = 12, res=500, units = "cm")

SecAct.CCC.heatmap(Seurat_obj, row.sorted=TRUE, column.sorted=TRUE, colors_cellType=my_cols)

dev.off()


png(paste0(outputPath,"/OV_circlize_Metastatic.png"), width = 9, height = 9, res=500, units = "cm")

SecAct.CCC.circle(Seurat_obj, colors_cellType=my_cols)

dev.off()


png(paste0(outputPath,"/OV_circlize_tumor_Metastatic.png"), width = 9, height = 9, res=500, units = "cm")

SecAct.CCC.circle(Seurat_obj, colors_cellType=my_cols, receiver="Tumor")

dev.off()




ccc <- Seurat_obj @misc $SecAct_output $SecretedProteinCCC
write.csv(ccc, paste0(outputPath,"/OV_SecretedProteinCCC.csv"), quote=F)
ccc <- read.csv(paste0(outputPath,"/OV_SecretedProteinCCC.csv"))


# download from https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/h.all.v2025.1.Hs.symbols.gmt
EMT <- c("ABI3BP","ACTA2","ADAM12","ANPEP","APLP1","AREG","BASP1","BDNF","BGN","BMP1","CADM1","CALD1","CALU","CAP2","CAPG","CCN1","CCN2","CD44","CD59","CDH11","CDH2","CDH6","COL11A1","COL12A1","COL16A1","COL1A1","COL1A2","COL3A1","COL4A1","COL4A2","COL5A1","COL5A2","COL5A3","COL6A2","COL6A3","COL7A1","COL8A2","COLGALT1","COMP","COPA","CRLF1","CTHRC1","CXCL1","CXCL12","CXCL6","CXCL8","DAB2","DCN","DKK1","DPYSL3","DST","ECM1","ECM2","EDIL3","EFEMP2","ELN","EMP3","ENO2","FAP","FAS","FBLN1","FBLN2","FBLN5","FBN1","FBN2","FERMT2","FGF2","FLNA","FMOD","FN1","FOXC2","FSTL1","FSTL3","FUCA1","FZD8","GADD45A","GADD45B","GAS1","GEM","GJA1","GLIPR1","GPC1","GPX7","GREM1","HTRA1","ID2","IGFBP2","IGFBP3","IGFBP4","IL15","IL32","IL6","INHBA","ITGA2","ITGA5","ITGAV","ITGB1","ITGB3","ITGB5","JUN","LAMA1","LAMA2","LAMA3","LAMC1","LAMC2","LGALS1","LOX","LOXL1","LOXL2","LRP1","LRRC15","LUM","MAGEE1","MATN2","MATN3","MCM7","MEST","MFAP5","MGP","MMP1","MMP14","MMP2","MMP3","MSX1","MXRA5","MYL9","MYLK","NID2","NNMT","NOTCH2","NT5E","NTM","OXTR","P3H1","PCOLCE","PCOLCE2","PDGFRB","PDLIM4","PFN2","PLAUR","PLOD1","PLOD2","PLOD3","PMEPA1","PMP22","POSTN","PPIB","PRRX1","PRSS2","PTHLH","PTX3","PVR","QSOX1","RGS4","RHOB","SAT1","SCG2","SDC1","SDC4","SERPINE1","SERPINE2","SERPINH1","SFRP1","SFRP4","SGCB","SGCD","SGCG","SLC6A8","SLIT2","SLIT3","SNAI2","SNTB1","SPARC","SPOCK1","SPP1","TAGLN","TFPI2","TGFB1","TGFBI","TGFBR3","TGM2","THBS1","THBS2","THY1","TIMP1","TIMP3","TNC","TNFAIP3","TNFRSF11B","TNFRSF12A","TPM1","TPM2","TPM4","VCAM1","VCAN","VEGFA","VEGFC","VIM","WIPF1","WNT5A")
EMT <- transferSymbol(EMT)


ccc_tumor_receiver <- ccc[ccc[,"receiver"]=="Tumor",]

ccc_tumor_receiver_stat <- sort(table(ccc_tumor_receiver[,2]),decreasing=T)
topSPs <- names(ccc_tumor_receiver_stat)[ccc_tumor_receiver_stat>1]

otherSPs <- ccc_tumor_receiver[!ccc_tumor_receiver[,"sender"]%in%c("Fibroblast","Macrophage"),"secretedProtein"]

fibroSPs <- ccc_tumor_receiver[ccc_tumor_receiver[,"sender"]%in%c("Fibroblast"),"secretedProtein"]
fibroSPs <- intersect(fibroSPs, EMT)

macroSPs <- ccc_tumor_receiver[ccc_tumor_receiver[,"sender"]%in%c("Macrophage"),"secretedProtein"]
macroSPs <- intersect(macroSPs, EMT)

# from fibroSPs and macroSPs
selectedSPs <- c("TIMP1","BGN","FAP","THBS2","SPARC","COL1A1","COL1A2","TGFBI","VCAN","POSTN","FBN1","MXRA5","HTRA1","BMP1","FSTL1")
selectedSPs <- intersect(selectedSPs, EMT)

secretedProtein <- sort(unique(c(topSPs, otherSPs, selectedSPs)))



sender <- unique(ccc_tumor_receiver[,1])
secretedProtein <- secretedProtein
receiver <- c("Tumor")

p <- SecAct.CCC.sankey(Seurat_obj, my_cols, sender=sender, secretedProtein=secretedProtein, receiver=receiver)
ggsave(paste0(outputPath,"/OV_sankey_Tumor.png"), p, width = 25, height = 20, dpi=300, units = "cm")
ggsave(paste0(outputPath,"/OV_sankey_Tumor.pdf"), p, width = 25, height = 20, dpi=300, units = "cm")


