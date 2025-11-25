rawPath <- "../_raw/"
dataPath <- "../data/"
resultsPath <- "../results/"

dataCrePath <- paste0(dataPath, "1_Creation/")
dataValPath <- paste0(dataPath, "2_Validation/")
dataAppPath <- paste0(dataPath, "3_Application/")

TSTAPath <- paste0(dataCrePath, "TSTA/")
visiumPath <- paste0(dataCrePath, "TSTA/data/")
meta <- read.csv(paste0(dataCrePath,"TSTA/meta.csv"),row.names=1)

HPAPath <- paste0(dataCrePath, "HPA/")
load(paste0(HPAPath,"HPA_20230612_gene_location.RData"))

NCBIPath <- paste0(dataCrePath, "NCBI/")
TCGAPath <- paste0(dataCrePath, "TCGA/")
ICGCPath <- paste0(dataCrePath, "ICGC/")
GEOPath <- paste0(dataCrePath, "GEO/")
finalSignaturesPath <- paste0(dataCrePath, "Signatures/")

data_MSigDB_path <- paste0(dataAppPath, "MSigDB/")

fig1Path <- paste0(resultsPath,"fig1/")
deconvResPath <- paste0(resultsPath,"preprocessDataDeconv/")

preprocessDataStatPath <- paste0(resultsPath,"preprocessDataStat/")
preprocessDataSummaryPath <- paste0(resultsPath,"preprocessDataSummary/")

signaturePath <- paste0(resultsPath,"signature/")
signatureCombPath <- paste0(resultsPath,"signatureComb/")
signatureComparePath <- paste0(resultsPath,"signatureCompare/")

QCPath <- paste0(resultsPath,"QC/")
QCFilterPath <- paste0(resultsPath,"QCFilter/")
QCSummaryPath <- paste0(resultsPath,"QCSummary/")

lambdaPath <- paste0(resultsPath,"lambda/")
lambdaSummaryPath <- paste0(resultsPath,"lambdaSummary/")

validationPath <- paste0(resultsPath,"validation/")
applicationPath <- paste0(resultsPath,"application/")

rawDataPath <- "../_raw/ST/"



lambdas <- c(10000,50000,100000,500000,1000000,5000000,10000000)

colors_surv <- c("#4575b4","#d73027")


sigNames <- c("SecAct","CytoSig","NicheNet.v1","NicheNet.v2","ImmuneDic")

sigColors <- c("#bdb5e1","#f9d580","#7ac7e2","#b0d992","#fa9fb5","#e3716e","#eca680","#54beaa")
names(sigColors) <- c(sigNames,"LigandExp","ReceptorExp","LRsumExp")

n_downsampling <- c(5,10,15,20,40,60,80,100,200,400,600,800,1000,2000,4000,6000,8000,10000)
n_downsampling_rep <- 10

notGenomeWide <- c(
"HNSCC_ICB_Foy2022",
"NSCLC_PD1_Foy2022",
"PanCancer_PD1_Prat2017"
)
notGenomeWide_alt <- c(
"HNSCC_ICB_Foy2022.OS",
"NSCLC_PD1_Foy2022.OS",
"PanCancer_PD1_Prat2017.PFS"
)
notGenomeWide_alt2 <- c(
"HeadNeck_ICB_Foy 2022",
"Lung (Non-Small Cell)_PD1_Foy 2022",
"PanCancer_PD1_Prat 2017"
)

CPTAC_rename <- data.frame(
	x=c("BRCA","COAD","GBM","HNSC","KIRC","LUAD","LUSC","OV","PDAC","UCEC"),
	y=c("Breast","Colon","Glioblastoma","HeadNeck","Kidney","Lung-Adeno","Lung-Squamous","Ovarian","Pancreatic","Endometrial")	
)

source("Secretome_s0_ifun.R")
