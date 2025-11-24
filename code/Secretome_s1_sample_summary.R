source("Secretome_s0_path.R")
	
### visium sample stat	
	
sts <- unique(meta[,1])
cancerTypes <- meta[,3]
cancerTypes_symbol <- sapply(strsplit(meta[,1],"_",fixed=T),function(x) return(x[1]))
cancerTypes_symbol_both <- paste0(cancerTypes, " (", cancerTypes_symbol,")")


meta_simple <- data.frame()
for(st in sts)
{
	meta_simple[st,"Cancer"] <- unique(meta[meta[,1]==st,3])
	meta_simple[st,"Sample"] <- nrow(meta[meta[,1]==st,,drop=F])
	
	if(grepl("10x",st)){
		meta_simple[st,"PMID"] <- "https://www.10xgenomics.com/resources/datasets"
	}else if(st%in%c("GBM_2025_DeJong","MM_2024_Sudupe")){
		meta_simple[st,"PMID"] <- paste0("https://www.biorxiv.org/content/10.1101/", unique(meta[meta[,1]==st,2]) )
	}else{
		meta_simple[st,"PMID"] <- paste0("https://pubmed.ncbi.nlm.nih.gov/", unique(meta[meta[,1]==st,2]) )
	}
}
write.csv(meta_simple,paste0(preprocessDataSummaryPath,"visium_stat_Table1.csv"),quote=F)
write.csv(meta_simple,paste0(visiumPath,"stat.csv"),quote=F)




fg.df <- as.data.frame(table(cancerTypes_symbol))

fg.df_alt <- fg.df
fg.df_alt[,2] <- fg.df_alt[,2]+0.5

library(ggplot2)
p <- ggplot(fg.df_alt, aes(cancerTypes_symbol, Freq)) +
    geom_bar(stat="identity",color="white", fill="#9abf88")+
    xlab(" ")+
    ylab("# Sample")+
    scale_y_continuous(trans='log2',breaks=c(2,5,10,20,50)) + 
    theme_bw()+ 
	theme(
		plot.background = element_blank(),
		panel.grid = element_blank(),
    	axis.title = element_text(color="black"),
    	axis.text.x = element_text(color="black", angle = 90, hjust = 1, vjust=0.5),
        axis.text.y = element_text(color="black", angle = 90, hjust = 0.5),
        legend.position = "none"
    )
ggsave(paste0(preprocessDataSummaryPath,"visium_stat_symbol.png"), p, width = 12.5, height = 7.5, dpi=400, units = "cm")




fg.df_comb <- fg.df
rownames(fg.df_comb) <- fg.df_comb[,1]
fg.df_comb[,1] <- as.character(fg.df_comb[,1])

tempList <- list(
	Brain=c("PCNSL","EPN","GBM","MB","PNST"),
	Nasopharyngeal=c("NPC"),
	Oral=c("THCA"),
	'Head and Neck'=c("HNAS","HNSC","OPC","OSCC"),
	Breast=c("BRCA"),
	Lung=c("LUAD","LUSC"),
	Gastric=c("STAD","GIST"),
	Liver=c("LIHC","CHOL","HB","GBC"),
	Kidney=c("KIRC","sRCC","WT"),
	Pancreatic=c("IPMN","PDAC","PanIN"),
	Ovarian=c("OV"),
	Colorectal=c("CRC"),
	Uterine=c("UCEC"),
	Bladder=c("BLCA"),
	Cervical=c("CESC"),
	Prostate=c("PRAD"),
	Skin=c("cSCC","SKCM"),
	'Soft tissue'=c("SARC","OS","MM")
)


fg.df_comb_new <- data.frame()
for(i in 1:length(tempList))
{
	temp <- tempList[[i]]
	tempP <- paste0(temp,collapse="/")
	fg.df_comb_new[names(tempList)[i],"Anatomical"] <- names(tempList)[i]
	fg.df_comb_new[names(tempList)[i],"CancerType"] <- tempP
	fg.df_comb_new[names(tempList)[i],"Count"] <- sum(fg.df_comb[temp,2])
}

write.csv(fg.df_comb_new,paste0(preprocessDataSummaryPath,"visium_biorender.csv"),quote=FALSE,row.names=FALSE)




fg.df <- as.data.frame(table(cancerTypes_symbol_both))
fg.df <- fg.df[order(fg.df[,2],decreasing=T),]
fg.df[,1] <- factor(fg.df[,1],levels=fg.df[,1])

library(ggplot2)
p <- ggplot(fg.df, aes(cancerTypes_symbol_both, Freq, label=cancerTypes_symbol_both)) +
    geom_col(width = 1, color = "white" ,fill="#C43E96", alpha=0.6) +
    geom_text(aes(y = 1), angle = 90, hjust = 0, size = 3) +
    #scale_fill_manual(values=mypalette)+
    ggtitle(" ")+
    xlab("Cancer type")+
    ylab(" ")+
    scale_y_continuous(limits = c(0, max(fg.df[,2])+5), expand = c(0, 0), position = "right")+
    theme_bw()+ 
	theme(
		plot.background = element_blank(),
		panel.grid = element_blank(),
    	axis.title.x = element_text(angle = 180), 
    	axis.title.y = element_text(angle = 270), 
    	axis.ticks.x = element_blank(), 
    	axis.text.x = element_blank(), 
        axis.text.y = element_text(color="black",angle = 90,vjust=0.5),
        legend.position = "none"
    )
    
ggsave(paste0(preprocessDataSummaryPath,"visium_stat_fullname.png"), p, width = 15.5, height = 10.5, dpi=400, units = "cm")
write.csv(fg.df,paste0(preprocessDataSummaryPath,"visium_stat_fullname.csv"),quote=F,row.names=F)


if(FALSE)
{
	# count FF and FFPE
	Preserv <- meta[,"Preservation"]
	fg.df <- as.data.frame(table(Preserv))
	fg.df[,1] <- as.character(fg.df[,1]) 
	
	FFPE_flag <- which(fg.df[,"Preserv"]=="Formalin-fixed paraffin-embedded (FFPE)")
	fg.df[FFPE_flag,"Preserv"] <- paste0("Formalin-fixed \n paraffin-embedded \n (FFPE, ", fg.df[FFPE_flag,"Freq"], ", ", round(fg.df[FFPE_flag,"Freq"]*100/sum(fg.df[,"Freq"]),1), "%)")
	FF_flag <- which(fg.df[,"Preserv"]=="Fresh frozen (FF)")
	fg.df[FF_flag,"Preserv"] <- paste0("Fresh frozen \n (FF, ", fg.df[FF_flag,"Freq"], ", ", round(fg.df[FF_flag,"Freq"]*100/sum(fg.df[,"Freq"]),1), "%)")
	
	
	library(ggplot2)
	library(scales)
	
	p1 <- ggplot(fg.df, aes(x="", y=Freq, fill=Preserv)) +
	  geom_bar(stat="identity", width=1, alpha=0.66) +
	  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
	  coord_polar("y", start=0)+
	  geom_text(aes(y = Freq/3 + c(0, cumsum(Freq)[-length(Freq)]), label = Preserv), color = "black", size=6)+
	  theme_void()+
	  theme(
			legend.position="none"
		)	 
	
	ggsave(paste0(preprocessDataSummaryPath,"visium_FF_FFPE.png"), p1, width = 12, height = 12, dpi=400, units = "cm")


	# comb gene and spot count
	smry <- data.frame()
	
	sts <- unique(meta[,"Study"])
	for(st in sts)
	{	
		preprocessDataStatPath.st <- paste0(preprocessDataStatPath,st,"/")
		dir.create(preprocessDataStatPath.st)
		
		sampleNames <- meta[meta[,"Study"]==st,"Sample_Name"]
		for(sampleName in sampleNames)
		{
			preprocessDataStatPath.st.sample <- paste0(preprocessDataStatPath.st,sampleName,"/")
			dir.create(preprocessDataStatPath.st.sample)
					
			stat <- read.csv(paste0(preprocessDataStatPath.st.sample,"stat_gene_spot.csv"), as.is=T, row.names=1, header=T)
			smry[rownames(stat),paste0(st,"@",sampleName)] <- stat[,1]
		}
	}
	
	smry <- t(as.matrix(smry))
	
	write.csv(smry,paste0(preprocessDataSummaryPath,"visium_UMI_Gene.csv"),quote=F)
	
	
	fg.df1 <- data.frame(x=smry[,3],y=smry[,4],z=smry[,2])
	fg.df1 <- fg.df1[order(fg.df1[,2]),]
	
	library(ggplot2)
	p1 <- ggplot(fg.df1, aes(x=x, y=y)) + 
		#geom_point(size=2, color="#f19670", alpha=0.5)+
		geom_point(aes(size=z), color="#f19670", alpha=0.5)+
		#geom_vline(xintercept=1000, color = "green",linewidth=0.8)+
		#geom_hline(yintercept=500, color = "green",linewidth=0.8)+
		xlab("Median UMI counts per spot")+
		ylab("Median gene counts per spot")+
		guides(size=guide_legend(title="Spot number"))+
		theme_bw()+ 
		theme(
			plot.background = element_blank(),
			panel.grid = element_blank(),
			axis.title = element_text(size=13,colour = "black"),
			axis.text = element_text(size=12,colour = "black"),
			legend.position = c(0.8, 0.3)
		)			
	ggsave(paste0(preprocessDataSummaryPath,"visium_UMI_Gene.png"), p1, width = 10, height = 9.5, dpi=400, units = "cm")
	
	
	fg.df1 <- data.frame(x=smry[,2])
	
	library(ggplot2)
	p1 <- ggplot(fg.df1,aes(x=x)) + 
		geom_histogram( colour="grey1", fill="#9abf88", alpha=0.5)+
		geom_density(color="grey3",linewidth=0.6)+
		xlab("# Spots detected under tissue")+
		ylab("# Sample")+
		theme_bw()+ 
		theme(
			panel.background = element_blank(),
			panel.grid = element_blank(),
			axis.title = element_text(size=14,colour = "black"),
			axis.text = element_text(size=13,colour = "black"),
			legend.position="none"
		)	
	ggsave(paste0(preprocessDataSummaryPath,"visium_SpotNum.png"), p1, width = 12, height = 12, dpi=400, units = "cm")

}
