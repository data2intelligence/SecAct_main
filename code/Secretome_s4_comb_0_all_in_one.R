source("Secretome_s0_path.R")

groupGenes <- SPs
groupName <- "SP"

#SWARM -t 2 -g 5 --time 03:00:00
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
	for(version in c("vst_free","vst","vst_condition_logUMI_cellType"))
	{
		comb <- data.frame()
		
		sts <- unique(meta[,"Study"])
		for(st in sts)
		{
			signaturePath.st <- paste0(signaturePath,st,"/")
			
			sampleNames <- meta[meta[,"Study"]==st,"Sample_Name"]
			for(sampleName in sampleNames)
			{	
				signaturePath.st.sample <- paste0(signaturePath.st,sampleName,"/")
				dir.create(signaturePath.st.sample)
				
				signaturePath.st.sample.singleSig <- paste0(signaturePath.st.sample,"singleSig_",version,"/")
				dir.create(signaturePath.st.sample.singleSig)

				if(file.exists(paste0(signaturePath.st.sample.singleSig,gene,".tsv.gz")))
				{
					m_sub <- read.table(paste0(signaturePath.st.sample.singleSig,gene,".tsv.gz"), sep="\t")
					comb[rownames(m_sub),paste0(st,"@",sampleName)] <- m_sub[,1]
				}
			} #sampleName
		} #st
	
		dir.create(paste0(signatureCombPath,"combSig_all_in_one_",groupName))
		saveRDS(comb, paste0(signatureCombPath,"combSig_all_in_one_",groupName,"/",gene,"_",version,".rds") )

	} # version
} #gene
