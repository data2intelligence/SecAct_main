dataPath <- "/data/rub2/project/Secretome/data/"
resultsPath <- "/data/rub2/project/Secretome/results/"

dataCrePath <- paste0(dataPath, "1_Creation/")
dataValPath <- paste0(dataPath, "2_Validation/")
dataAppPath <- paste0(dataPath, "3_Application/")

finalSignaturesPath <- paste0(dataCrePath, "Signatures/")

TSTAPath <- paste0(dataCrePath, "TSTA/")
NCBIPath <- paste0(dataCrePath, "NCBI/")
TCGAPath <- paste0(dataCrePath, "TCGA/")
ICGCPath <- paste0(dataCrePath, "ICGC/")
HPAPath <- paste0(dataCrePath, "HPA/")
GEOPath <- paste0(dataCrePath, "GEO/")

visiumPath <- paste0(dataCrePath, "TSTA/data/")

meta <- read.csv(paste0(dataCrePath,"TSTA/meta.csv"),row.names=1)
load(paste0(HPAPath,"HPA_20230612_gene_location.RData"))

rawDataPath <- "../_raw/ST/"

fig1Path <- "../results/fig1/"
deconvResPath <- "../results/preprocessDataDeconv/"

preprocessDataStatPath <- "../results/preprocessDataStat/"
preprocessDataSummaryPath <- "../results/preprocessDataSummary/"

signaturePath <- "../results/signature/"
signatureCombPath <- "../results/signatureComb/"
signatureComparePath <- "../results/signatureCompare/"

QCPath <- "../results/QC/"
QCFilterPath <- "../results/QCFilter/"
QCSummaryPath <- "../results/QCSummary/"

lambdaPath <- "../results/lambda/"
lambdaSummaryPath <- "../results/lambdaSummary/"

validationPath <- "../results/validation/"
applicationPath <- "../results/application/"
predictionPath <- "../results/application/prediction/"

MSigDBPath <- "../results/MSigDB/"
MSigDB_genePath <- "../results/MSigDB_gene/"

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
