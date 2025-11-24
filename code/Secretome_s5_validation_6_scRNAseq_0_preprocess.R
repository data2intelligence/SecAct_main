source("Secretome_s0_path.R")

#SWARM -t 2 -g 200 --time 10:00:00
args = commandArgs(trailingOnly=TRUE)
cancer <- args[1]

inputPath <- paste0("../data/2_Validation/scRNAseq/")
outputPath <- paste0(validationPath,"scRNAseq/")
dir.create(outputPath)

transformFormat <- paste0(
	"python Secretome_s5_validation_6_scRNAseq_0_preprocess.py ",
	inputPath,cancer,".pickle.gz ",
	outputPath,cancer,"/expr.csv.gz"
	)
	
system(transformFormat)
