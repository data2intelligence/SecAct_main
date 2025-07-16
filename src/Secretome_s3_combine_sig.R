source("Secretome_s0_path.R")
source("../../SpaCE/code/ifun.R")
source("importHPA.R")

groupGenes <- MPs
groupName <- "MP"

groupGenes <- SPs
groupName <- "SP"

#SWARM -t 2 -g 5 --time 02:00:00
args = commandArgs(trailingOnly=TRUE)
id <- args[1]

id <- as.numeric(id)
id_end <- id+9

if(id_end>length(groupGenes))
{
	id_end <- length(groupGenes)
}

for(gene in groupGenes[id:id_end])
{
	for(version in c("vst","weighted","weighted2","vst_condition_logUMI_cellType","vst_condition_cellType"))
	{
		comb <- data.frame()
		
		sts <- unique(meta[meta[,"Method"]=="Visium",1])
		for(st in sts)
		{
			signaturePath.st <- paste0(signaturePath,st,"/")
			deconvResPath.st <- paste0(deconvResPath,st,"/")
			
			sampleNames <- meta[meta[,"Study"]==st&meta[,"Method"]=="Visium"&meta[,"Organism"]=="Human","Sample_Name"]
			for(sampleName in sampleNames)
			{	
				signaturePath.st.sample <- paste0(signaturePath.st,sampleName,"/")
				deconvResPath.st.sample <- paste0(deconvResPath.st,sampleName,"/")
				
				if(version%in%c("vst","vst_condition_logUMI_cellType","vst_condition_cellType"))
				{
					signaturePath.st.sample.singleSig <- paste0(signaturePath.st,sampleName,"/singleSig_mean_",groupName,"_s0_r3_",version,"/")
				
					fileName <- paste0(signaturePath.st.sample.singleSig,gene,".tsv")
					if(file.exists(fileName))
					{
						m_sub <- read.table(fileName, sep="\t")
						comb[rownames(m_sub),paste0(st,"@",sampleName)] <- m_sub[,1]
					}
				}
				
				if(version%in%c("weighted","weighted2"))
				{
					signaturePath.st.sample.singleSig <- paste0(signaturePath.st,sampleName,"/singleSig_mean_",groupName,"_s0_r3_vst/")
				
					fileName <- paste0(signaturePath.st.sample.singleSig,gene,".tsv")
					if(file.exists(fileName))
					{
						m_sub <- read.table(fileName, sep="\t")
						
						weight_vec <- read.table(paste0(deconvResPath.st.sample,version,".tsv"), sep="\t")
						weight_vec <- weight_vec[,1]
				
						comb[rownames(m_sub),paste0(st,"@",sampleName)] <- m_sub[,1]*weight_vec
					}
				}
				
			} #sampleName
		} #st
		
		dir.create(paste0(signatureCombPath,"combSig_mean_",groupName,"_s0_r3"))
		write.table(comb, paste0(signatureCombPath,"combSig_mean_",groupName,"_s0_r3/",gene,"_",version,".tsv"), quote=F, sep="\t")
	
	} # version
} #gene
