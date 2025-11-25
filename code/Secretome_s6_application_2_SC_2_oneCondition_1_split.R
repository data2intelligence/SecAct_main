source("Secretome_s0_path.R")	

inputPath <- paste0(dataAppPath,"scRNAseq_PanCancer/")
outputPath <- paste0(applicationPath,"scRNAseq_PanCancer/")
dir.create(outputPath)

#download data
system(paste0("curl -L https://zenodo.org/records/15554080/files/PanCancer_igt_s9_fine_counts.h5ad -o ", inputPath, "PanCancer_igt_s9_fine_counts.h5ad"))


library(anndata)
adata <- read_h5ad(paste0(inputPath,"PanCancer_igt_s9_fine_counts.h5ad"))
aa <- adata$obs

table(aa[,"tissue"])
#     Adjacent         Blood    Metastasis PleuralFluids     PreLesion 
#       390119        459866        379126         17749        169382 
#        Tumor 
#      2730733 

table(aa[,"preLesion"])
##        CAG          CG         CIS        HSIL         IAC         IMS 
##      14155       15943        2840       10240       37305        8163 
##        IMW Leukoplakia         NAG         NLH          No         OLK 
##       2351        5962        4242       14600     3977593       19898 
##      Polyp         iCD 
##      10239       23444

table(aa[,"metastasis"])
##     Abdomen        Brain    ChestWall        Liver         Lung    LymphNode 
##        3077        26832        18033       166370         2025       142016 
##          No      Omentum         PVTT   Peritoneal Subcutaneous 
##     3767849           84         5505         6286         8898

table(aa[,"tumorPhase"])
## Carcinoma In Situ                 I                IA               IA2 
##              2840            689707             42274             20856 
##               IA3                IB               IB1               IC2 
##             24220             63904             46931              7777 
##                II               IIA               IIB               III 
##            334271             29925             50920            454954 
##              IIIA              IIIB              IIIC                IV 
##             41975             47675             61019            302402 
##               IVA               IVB           Unknown 
##            129238             84506           1711581

table(aa[,"recurrence"])
##     No Unknown     Yes 
## 618348 3470702   57925

table(aa[,"treatment"])
##                                                      2 times (GT, GT) 
##                                                                  14266 
##                                                   3 times (GT, GT, GT) 
##                                                                  13612 
##                                                 3 times (MTX, IFO, AP) 
##                                                                  30959 
##                                               4 times (AP, IE, AP, IE) 
##                                                                   8157 
##                                            4 times (IFO, AP, MTX, MTX) 
##                                                                   8419 
##                                             4 times (MTX, AP, IFO, AP) 
##                                                                   4665 
##                                            4 times (MTX, AP, IFO, MTX) 
##                                                                   9345 
##                                   6 times (MTX, AP, MTX, AP, MTX, MTX) 
##                                                                   3807 
##                                                      Anti-PD-L1+ Chemo 
##                                                                  43562 
##                                                      Anti-PD-L2+ Chemo 
##                                                                  18549 
##                                                      Anti-PD-L3+ Chemo 
##                                                                  14151 
##                                                      Anti-PD-L4+ Chemo 
##                                                                  23487 
##                                                      Anti-PD-L5+ Chemo 
##                                                                  12311 
##                                                      Anti-PD-L6+ Chemo 
##                                                                   8775 
##                                                      Anti-PD-L7+ Chemo 
##                                                                  25036 
##                                                      Anti-PD-L8+ Chemo 
##                                                                  11602 
##                                                      Anti-PD-L9+ Chemo 
##                                                                  19437 
##                                                     Anti-PD-L10+ Chemo 
##                                                                  16340 
##                                                     Anti-PD-L11+ Chemo 
##                                                                  18508 
##                                   Bilateral orchidectomy, bicalutamide 
##                                                                   9175 
##                                                                  Chemo 
##                                                                 188762 
##                                                              Docetaxel 
##                                                                   7603 
##                                                                FFX(7M) 
##                                                                  18728 
##                                                        FFX,G/A(17M,3M) 
##                                                                   2636 
##                                                        FFX,G/A(20M,2M) 
##                                                                   3216 
##                                                           FFX/SBRT(5M) 
##                                                                   3932 
##                                                               FFX（3M) 
##                                                                   4552 
##                                                               FFX（4M) 
##                                                                  14003 
##                                                                G/A(4M) 
##                                                                   6879 
##                                                Goserelin, bicalutamide 
##                                                                   6353 
##                                 Hemi-thyroidectomy and TSH suppression 
##                                                                   7511 
##                                    Iodine ablation and TSH suppression 
##                                                                   2516 
##                                                                     No 
##                                                                2492461 
##                                Total thyroidectomy and TSH suppression 
##                                                                   4499 
##Total thyroidectomy, three times of iodine ablation and TSH suppression 
##                                                                  12512 
##                                                                Unknown 
##                                                                1024129 
##                                                                    Yes 
##                                                                   1807 
##                                                        aPD1 + VEGF TKI 
##                                                                   2190 
##                                                                  aPD-1 
##                                                                   4095 
##                                                        aPD-1 + aCTLA-4 
##                                                                   5493 
##                                                               anti-PD1 
##                                                                  18935 

table(aa[,"treatmentResponse"])
##     NE      No      PD      PR      SD Unknown 
##    165 2397409   18583  209063  179854 1341901

table(aa[,"treatmentPhase"])
##      After    Baseline          No Progression     Unknown 
##     355270      181437     2397409       66640     1146219 
     
table(aa[,"sampleType"])
##   Fresh 
## 4146975

table(aa[,"cellSort"])
##    Mix   Total 
## 195823 3951152

cohort_vec <- as.character(adata$obs[,"cohortName"])
cancers <- unique(cohort_vec)

cancers
##  [1] "HGSOC_GSE184880"       "HCC_GSE149614"         "ICC_GSE138709"        
##  [4] "PTC_GSE184362"         "ccRCC_GSE207493"       "Melanoma_GSE215120"   
##  [7] "NPC_GSE150430"         "OS_GSE152048"          "CC_E-MTAB-11948"      
## [10] "LUAD_GSE189357"        "PCa_GSE141445"         "PDAC_GSE205013"       
## [13] "Pre_CRC_GSE134809"     "GC_GSE183904"          "TNBC_GSE169246"       
## [16] "TGCT_GSE197778"        "ESCC_GSE160269"        "GCTB_GSE168664"       
## [19] "OSCC_GSE172577"        "UVM_GSE138433"         "NB_GSE137804"         
## [22] "NET_GSE140312"         "PDAC_GSE154778"        "LUAD_GSE131907"       
## [25] "UVM_GSE139829"         "RCC_GSE152938"         "PDAC_GSE155698"       
## [28] "NPC_GSE150825"         "Pre_BLCA_GSE225190"    "NPC_GSE162025"        
## [31] "PanCancer_E-MTAB-8107" "RCC_SCP1288"           "PCa_GSE137829"        
## [34] "BC_GSE148673"          "cSCC_GSE144236"        "PDAC_PRJCA001063"     
## [37] "LUAD_HRA000154"        "RCC_GSE159115"         "CRC_SYN26844071"      
## [40] "PLC_PRJCA007744"       "Pre_GC_GSE134520"      "Pre_OSCC_HRA001006"   
## [43] "Pre_CRC_GSE161277"     "Pre_CC_GSE208653"      "Pre_HNSCC_GSE181919"  

length(cancers)
# 45

cellNumber <- data.frame()
tissueNumber <- data.frame()
for(cancer in cancers)
{
	# Subset and copy to a new AnnDataR6 object
	adata_subset <- adata[cohort_vec == cancer, ]
	
	print(cancer)
	if("Metastasis"%in%unique(adata_subset$obs[,"tissue"]))
	{
		print(table(adata_subset$obs[,"tissue"]))
	}
	
	cellNumber[cancer,"count"] <- nrow(adata_subset)
	tissueNumber[cancer,"count"] <- length(table(adata_subset$obs[,"tissue"]))
	
	# Save to h5ad
	write_h5ad(adata_subset, paste0(outputPath,cancer,".h5ad"))
}
write.csv(cellNumber, paste0(outputPath, "cellNumber.csv"), quote=FALSE)
write.csv(tissueNumber, paste0(outputPath, "tissueNumber.csv"), quote=FALSE)

