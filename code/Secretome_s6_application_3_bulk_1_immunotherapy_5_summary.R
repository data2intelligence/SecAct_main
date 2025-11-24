source("Secretome_s0_path.R")

inputPath <- paste0(dataAppPath,"Immunotherapy/")
outputPath <- paste0(applicationPath,"Immunotherapy/")

smy <- read.csv(paste0(outputPath,"prediction_act.csv"),as.is=T,row.names=1,header=T)
smy_l <- smy[, setdiff(colnames(smy),notGenomeWide_alt) ]
smy_l.median <- apply(smy_l,1,function(x) median(x,na.rm=T))

stat <- data.frame()
for(n in n_downsampling)
{
	for(i in 1:n_downsampling_rep)
	{
		smy_s <- read.csv(paste0(outputPath,"downsampling/",n,"/prediction_act_",i,".csv"),as.is=T,row.names=1,header=T)
		smy_s.median <- apply(smy_s,1,function(x) median(x,na.rm=T))
		
		fg.df <- data.frame(
			x=smy_l.median,
			y=smy_s.median
		)
		
		cor_res <- cor.test(fg.df[,1],fg.df[,2])
		rv <- round(cor_res$estimate,2)
		pv <- signif(cor_res$p.value,2)
		
		stat[i,as.character(n)] <- rv
		
		if(i==2&n%in%c(5,10,20,100,1000))
		{
			library(ggplot2)
			p1 <- ggplot(fg.df,aes(x=x, y=y)) + 
				geom_point(alpha=0.3, size=0.6, color="skyblue")+
				annotate("text", x = 1, y=-0.7, label = paste0("r = ",rv))+
				xlab("Median z (Genome-wide)")+
				ylab("Median z (Downsampling)")+
				ggtitle(paste0("Gene coverage: ",n))+
				theme_classic()+ 
				theme(
					plot.background = element_blank(),
					panel.grid = element_blank(),
					plot.title = element_text(size = 10),
					legend.position = "none",
					legend.title = element_blank(),
					legend.key.size = unit(0.7, 'lines')
				)			
			ggsave(paste0(outputPath,"prediction_genome_downsample_compare_",n,".png"), p1, width = 6.2, height = 6.2, dpi=500, units = "cm")
		}
	}
}


fg.df <- reshape2::melt(as.matrix(stat))
fg.df[,2] <- as.character(fg.df[,2])
fg.df[,2] <- factor(fg.df[,2], levels=as.character(n_downsampling) )

library(ggplot2)
p1 <- ggplot(fg.df,aes(x=Var2,y=value)) + 
	geom_hline(yintercept=0, color = "grey", linewidth=0.8, linetype="dashed")+
	geom_boxplot( aes(color=Var2), alpha=0, width=0.5, outlier.shape = NA)+
	geom_jitter(color="darkgrey",alpha=0.3, size=1, width=0.1)+
	scale_y_continuous(breaks = seq(0, 1, by = 0.2))+
	ylab("Correlation r")+
	xlab("Gene coverage")+
	theme_classic()+ 
	theme(
		panel.background = element_blank(),
		panel.grid = element_blank(),
		axis.title = element_text(colour = "black"),
		axis.text = element_text(colour = "black"),
		axis.text.x = element_text(size=10),
		legend.position="none"
	)
	
ggsave(paste0(outputPath,"prediction_genome_downsample_compare_summary.png"), p1, width = 26, height = 6, dpi=500, units = "cm")



