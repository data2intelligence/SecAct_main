source("Secretome_s0_path.R")
source("../../SpaCE/code/ifun.R")
library(Seurat)
library(ggplot2)
library(patchwork)
library(SecAct)

cancer <- "OV"

Seurat_obj <- readRDS(paste0(applicationPath, "OV_618/OV_scRNAseq.rds"))

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


p1 <- DimPlot(Seurat_obj, reduction = "umap", cells = cellID_cellType[,1], cols = my_cols, group.by = "MyCellType")+NoLegend()+ggtitle(' \n \n ')+NoAxes()
p2 <- DimPlot(Seurat_obj, reduction = "umap", cells = sample(Cells(Seurat_obj)), group.by = "Groups")+ NoLegend()+ggtitle(' \n \n ')+NoAxes()

ggsave(paste0(applicationPath,cancer,"/OV_UMAP.png"), p1/p2, width = 12, height = 25, dpi=500, units = "cm")



p1 <- DimPlot(Seurat_obj, reduction = "umap", group.by = "MyCellType", split.by="Groups", cols = my_cols)+NoLegend()+ggtitle(NULL)+NoAxes()
ggsave(paste0(applicationPath,cancer,"/OV_UMAP2.png"), p1, width = 24, height = 12, dpi=500, units = "cm")



p1 <- DimPlot(Seurat_obj, reduction = "umap", group.by = "Patients")

ggsave(paste0(applicationPath,cancer,"/OV_UMAP3.png"), p1, width = 14, height = 12, dpi=500, units = "cm")




cell_stat <- data.frame()
Patients <- c("HGSOC1","HGSOC3","HGSOC4","HGSOC6")
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


ref <- paste0("/data/rub2/project/Secretome/results/signatureComb/AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3_1k_vst_condition_logUMI_cellType_0.9.tsv")


Patients <- c("Metastatic","Primary","HGSOC1","HGSOC3","HGSOC4","HGSOC6")
Patients <- c("HGSOC1","HGSOC3","HGSOC4","HGSOC6")
Patients <- c("Metastatic")

for(Patient in Patients)
{
	if(Patient%in%c("Metastatic","Primary"))
	{
		Seurat_obj_patient <- Seurat_obj
	}else{
		Seurat_obj_patient <- subset(Seurat_obj, Patients %in% Patient )
	}
	
	if(Patient%in%c("Primary"))
	{
		conditionCase="Primary"
		conditionControl="Metastatic"
	}else{
		conditionCase="Metastatic"
		conditionControl="Primary"
	}	
	
	Seurat_obj_patient <- SecAct.CCC.scRNAseq(
	  Seurat_obj_patient, 
	  cellType_meta="MyCellType",
	  condition_meta="Groups", 
	  conditionCase=conditionCase, 
	  conditionControl=conditionControl,
	  sigMatrix=ref,
	  act_diff_cutoff=2,
	  exp_logFC_cutoff=0.2,
	  exp_mean_all_cutoff=2,
	  exp_fraction_case_cutoff=0.1,
	  padj_cutoff=0.01
	)  
	
	
	png(paste0(applicationPath,cancer,"/OV_heatmap_",Patient,".png"), width = 12, height = 12, res=500, units = "cm")
	
	SecAct.CCC.heatmap(Seurat_obj_patient, row.sorted=TRUE, column.sorted=TRUE, colors_cellType=my_cols)
	
	dev.off()
	
	
	png(paste0(applicationPath,cancer,"/OV_circlize_",Patient,".png"), width = 9, height = 9, res=500, units = "cm")
	
	SecAct.CCC.circle(Seurat_obj_patient, colors_cellType=my_cols)
	
	dev.off()
	
	
	png(paste0(applicationPath,cancer,"/OV_circlize_tumor_",Patient,".png"), width = 9, height = 9, res=500, units = "cm")
	
	SecAct.CCC.circle(Seurat_obj_patient, colors_cellType=my_cols, receiver="Tumor")
	
	dev.off()
	
	saveRDS(Seurat_obj_patient, file = paste0(applicationPath,cancer,"/OV_",Patient,"_Seurat_obj.rds"))

}




if(FALSE)
{

Seurat_obj_patient <- readRDS(file = paste0(applicationPath,cancer,"/OV_",Patient,"_Seurat_obj.rds"))

ccc <- Seurat_obj_patient @misc $SecAct_output $SecretedProteinCCC
write.csv(ccc, paste0(applicationPath,cancer,"/SecretedProteinCCC.csv"), quote=F)


gmt <- read.gmt(paste0("/data/rub2/data/MSigDB/h.all.v2025.1.Hs.symbols.gmt"))
EMT <- gmt[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]]
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

p <- SecAct.CCC.sankey(Seurat_obj_patient, my_cols, sender=sender, secretedProtein=secretedProtein, receiver=receiver)
ggsave(paste0(applicationPath,cancer,"/OV_sankey_Tumor.png"), p, width = 25, height = 20, dpi=300, units = "cm")
ggsave(paste0(applicationPath,cancer,"/OV_sankey_Tumor.pdf"), p, width = 25, height = 20, dpi=300, units = "cm")

}





# compare patient-specific CCC
if(FALSE)
{


Patients <- c("HGSOC1","HGSOC3","HGSOC4","HGSOC6")
fg.df <- data.frame()
CCC_list <- list()
for(Patient in Patients)
{
	Seurat_obj_patient <- readRDS(paste0(applicationPath,cancer,"/OV_",Patient,"_Seurat_obj.rds"))
	
	CCC <- Seurat_obj_patient@ misc$ SecAct_output$ SecretedProteinCCC
	
	CCC_list[[Patient]] <- rownames(CCC)
	
	fg.df[Patient,"Patient"] <- Patient
	fg.df[Patient,"CCC"] <- nrow(CCC)
}


library(ggplot2)
p0 <- ggplot(fg.df, aes(x=Patient, y=CCC, fill=Patient)) +
  geom_bar(stat="identity", position="dodge", color="white", width=0.66, alpha=0.7) +
  ylab("# Cell-cell communication")+
  theme_classic()+
  theme(
	  panel.grid = element_blank(),
	  panel.background = element_blank(),
	  axis.text = element_text(colour = "black"),
	  axis.title = element_text(colour = "black"),
	  axis.title.x = element_blank(),
	  strip.background = element_rect(colour="white", fill="white"),
	  legend.position = "none"
	)
	
ggsave(paste0(applicationPath,cancer,"/OV_patient_CCC_count.png"), p0, width = 8, height = 6, dpi=400, units = "cm")


library(ggVennDiagram)

p0 <- ggVennDiagram(CCC_list)
ggsave(paste0(applicationPath,cancer,"/OV_patient_CCC_olp.png"), p0, width = 16, height = 12, dpi=400, units = "cm")



}