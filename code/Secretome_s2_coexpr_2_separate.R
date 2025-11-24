source("Secretome_s0_path.R")

#SWARM -t 10 -g 50 --time 00:30:00
args = commandArgs(trailingOnly=TRUE)
st <- args[1]
sampleName <- args[2]

sts <- unique(meta[,"Study"])
#for(st in sts)
#{	
	signaturePath.st <- paste0(signaturePath,st,"/")
	dir.create(signaturePath.st)
	
	sampleNames <- meta[meta[,"Study"]==st,"Sample_Name"]
	#for(sampleName in sampleNames)
	#{	
		signaturePath.st.sample <- paste0(signaturePath.st,sampleName,"/")
		dir.create(signaturePath.st.sample)
		
		for(version in c("vst_free","vst","vst_condition_logUMI_cellType"))
		{
			signaturePath.st.sample.singleSig <- paste0(signaturePath.st.sample,"singleSig_",version,"/")
			dir.create(signaturePath.st.sample.singleSig)

			st.matrix.data.vst <- as.matrix(read.table(gzfile(paste0(signaturePath.st.sample, version, ".tsv.gz")), sep="\t",check.names=F))
			
			W <- calWeights(colnames(st.matrix.data.vst), radius=200, sigma=100, diagAsZero=TRUE)
					
			sigmat <- spatialCrossCorrelation(st.matrix.data.vst, W)
			
			sigmat <- sigmat[,colnames(sigmat)%in%SPs]		
			sigmat <- sigmat - rowMeans(sigmat, na.rm=TRUE)
			
			#for(gene in colnames(sigmat))
			runSeparate <- function(gene)
			{
				write.table(sigmat[,gene,drop=F], gzfile(paste0(signaturePath.st.sample.singleSig, gene, ".tsv.gz")), quote=F, sep="\t")
			}
			
			parallel::mclapply(colnames(sigmat), runSeparate, mc.cores=10)
		}
#	}
#}
