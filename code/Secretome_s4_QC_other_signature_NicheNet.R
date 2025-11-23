source("importHPA.R")


# v1
NicheNet <- readRDS("/data/rub2/project/Secretome/data/NicheNet_ligand_target_matrix.rds")

ligand_0 <- colnames(NicheNet)[colSums(NicheNet)==0] # all value in this column are 0\
ligand_0

NicheNet <- NicheNet[,!colnames(NicheNet)%in%ligand_0] 
NicheNet <- NicheNet[,colnames(NicheNet)%in%SPs]

rownames(NicheNet) <- transferSymbol(rownames(NicheNet))
colnames(NicheNet) <- transferSymbol(colnames(NicheNet))
NicheNet <- rm_duplicates(NicheNet)

write.table(NicheNet,"/data/rub2/project/Secretome/code/CytoSig-master/CytoSig/signature.centroid.11",quote=F,sep="\t")

varVect <- sort(apply(NicheNet,1,var),decreasing=T)
NicheNet <- NicheNet[names(varVect)[1:10000],]

write.table(NicheNet,"/data/rub2/project/Secretome/code/CytoSig-master/CytoSig/signature.centroid.111",quote=F,sep="\t")


# v2
NicheNet <- readRDS("/data/rub2/project/Secretome/data/NicheNet_ligand_target_matrix_nsga2r_final.rds") 

ligand_0 <- colnames(NicheNet)[colSums(NicheNet)==0] # all value in this column are 0\
ligand_0

NicheNet <- NicheNet[,!colnames(NicheNet)%in%ligand_0] 
NicheNet <- NicheNet[,colnames(NicheNet)%in%SPs]

rownames(NicheNet) <- transferSymbol(rownames(NicheNet))
colnames(NicheNet) <- transferSymbol(colnames(NicheNet))
NicheNet <- rm_duplicates(NicheNet)

NicheNet <- NicheNet[rowSums(NicheNet)>0,]

write.table(NicheNet,"/data/rub2/project/Secretome/code/CytoSig-master/CytoSig/signature.centroid.112",quote=F,sep="\t")

