transform_barcode_to_spotID <- function(barcodes)
{
	barcode_v1 <- read.csv(paste0("../_raw/ST/BRCA_10x_Datasets/Version1.0.0_Breast.Cancer_rep1/spatial/tissue_positions_list.csv"),as.is=TRUE,row.names=1,header=FALSE)
	barcode_v2 <- read.csv(gzfile(paste0("../_raw/ST/HB_2025_Munter/GSM8657788_H1227_tissue_positions.csv.gz")),as.is=TRUE,row.names=1,header=TRUE)				
	
	if(all(barcodes %in% rownames(barcode_v1)))
	{
		barcode <- barcode_v1
	}else if(all(barcodes %in% rownames(barcode_v2))){
		barcode <- barcode_v2
	}else{
		stop("Wrong barcode version!")
	}
	
	barcode[["comb"]] <- paste0(barcode[,2],"x",barcode[,3])
	
	barcode[barcodes,"comb"]
}

transferSymbol <- function(x, dataPath=NCBIPath)
{
	alias2symbol <- read.csv(gzfile(paste0(dataPath,"NCBI_20251008_gene_result_alias2symbol.csv.gz")))
	alias2symbol[is.na(alias2symbol[,"Alias"]),"Alias"] <- "NA"
	
	x[x%in%alias2symbol[,1]] <- alias2symbol[
		match(
			x[x%in%alias2symbol[,1]],
			alias2symbol[,1]
		), 2]
	
	x
}

rm_duplicates <- function(mat)
{
  gene_count <- table(rownames(mat))
  gene_dupl <- names(gene_count)[gene_count>1]

  if(length(gene_dupl) > 0){
    gene_unique <- names(gene_count)[gene_count==1]
    gene_unique_index <- which(rownames(mat)%in%gene_unique)

    gene_dupl_index <- c()
    for(gene in gene_dupl)
    {
      gene_dupl_index_gene <- which(rownames(mat)%in%gene)
      mat_dupl_gene <- mat[gene_dupl_index_gene,,drop=FALSE]
      dupl_sum <- Matrix::rowSums(mat_dupl_gene)
      max_flag <- which(dupl_sum==max(dupl_sum))
      gene_dupl_index <- c(gene_dupl_index,gene_dupl_index_gene[max_flag[1]])
    }

    mat <- mat[sort(c(gene_unique_index,gene_dupl_index)),,drop=FALSE]
  }

  return(mat)
}


sweep_sparse <- function(m, margin, stats, fun)
{
  f <- match.fun(fun)

  if(margin==1)
  {
    idx <- m@i + 1
  }else{
    if(class(m)[1]=="dgCMatrix")
    {
      idx <- rep(1:m@Dim[2], diff(m@p))
    }else{
      idx <- x@j + 1
    }
  }

  m@x <- f(m@x, stats[idx])
  m
}


calWeights_old <- function(SpotIDs, r, diag0=TRUE)
{
	d <- matrix(Inf,ncol=length(SpotIDs),nrow=length(SpotIDs))
	colnames(d) <- SpotIDs
	rownames(d) <- SpotIDs
	
	for(i in 1:ncol(d))
	{
		xy <- rownames(d)[i]
		x <- as.numeric(unlist(strsplit(xy,"x")))[1]
		y <- as.numeric(unlist(strsplit(xy,"x")))[2]
		xm <- r-1
		ym <- 2*(r-1)+1
		x_y <- expand.grid((x-xm):(x+xm),(y-ym):(y+ym))
		x_y_d <- cbind(x_y,d=sqrt( (0.5*sqrt(3)*(x_y[,1]-x))^2 + (0.5*(x_y[,2]-y))^2) )
		rownames(x_y_d) <- paste0(x_y[,1],"x",x_y[,2])
		x_y_d <- x_y_d[rownames(x_y_d)%in%rownames(d),]
			
		d[xy,rownames(x_y_d)] <- x_y_d[,3]
		d[rownames(x_y_d),xy] <- x_y_d[,3]
	}
	
	W <- exp(-d^2/2)
	
	if(diag0==TRUE) diag(W) <- 0
	
	W <- W[,colSums(W)>0] # remove spot island
	W <- W[rowSums(W)>0,] # remove spot island
	
	W
}


calWeights <- function(SpotIDs, radius=200, sigma=100, diagAsZero=TRUE)
{
	spotCoordinates <- t(matrix(as.numeric(unlist(strsplit(SpotIDs,"x"))),nrow=2))
	rownames(spotCoordinates) <- SpotIDs
	colnames(spotCoordinates) <- c("array_row","array_col")
	
	# transform array ID to coordinates (um)
	spotCoordinates[,1] <- spotCoordinates[,1] * 0.5 * sqrt(3) * 100
	spotCoordinates[,2] <- spotCoordinates[,2] * 0.5 * 100
	
	nn_result <- RANN::nn2(spotCoordinates, searchtype="radius", radius=radius, k=nrow(spotCoordinates))

	neighbor_indices <- nn_result$nn.idx
	neighbor_distances <- nn_result$nn.dists
	
	i <- rep(1:nrow(neighbor_indices), each=ncol(neighbor_indices)) # row indices (cell index)
	j <- as.vector(t(neighbor_indices)) 
	x <- as.vector(t(neighbor_distances))
	
	valid <- x<=radius & x>0
	i <- i[valid]          # Keep only valid indices
	j <- j[valid]          # Valid neighbor indices
	x <- x[valid]          # Valid distances
	
	# transform distance to weight
	x <- exp( -x^2 / (2*sigma^2) ) 
	
	# Create the sparse matrix using the 'i', 'j', and 'x' vectors
	library(Matrix)
	W <- sparseMatrix(i=i, j=j, x=x, dims=c(nrow(neighbor_indices), nrow(neighbor_indices)), repr="T")
	rownames(W) <- rownames(spotCoordinates)
	colnames(W) <- rownames(spotCoordinates)
	
	if(diagAsZero==FALSE) diag(W) <- 1
	
	W
}


spatialCrossCorrelation <- function(mat, W)
{
	W <- W[,colSums(W)>0] # remove spot island
	W <- W[rowSums(W)>0,] # remove spot island
	
	mat <- mat[,colnames(W)]
	
	N <- ncol(mat)
	
	for (i in 1:nrow(mat))
	{
		x <- mat[i, ]
		
		dx <- x - mean(x, na.rm=TRUE)
		stdx <- sum(dx^2, na.rm=TRUE)
		stdx <- sqrt(stdx/N)
		dx <- dx/stdx
		
		dx -> mat[i, ]
	}
		
	XW <- as.matrix(mat %*% W)
	
  	XWX <- tcrossprod(XW, mat)
	
	m <- XWX/sum(W)
	
	m
}


spatialAutoCorrelation <- function(mat, W, permuteNo=1000)
{
	W <- W[,colSums(W)>0] # remove spot island
	W <- W[rowSums(W)>0,] # remove spot island
	
	mat <- mat[,colnames(W)]
	
	N <- ncol(mat)
	
	for (i in 1:nrow(mat))
	{
		x <- mat[i, ]
		
		dx <- x - mean(x, na.rm=TRUE)
		stdx <- sum(dx^2, na.rm=TRUE)
		stdx <- sqrt(stdx/N)
		dx <- dx/stdx
		
		dx -> mat[i, ]
	}
	
	set.seed(123456)
	MoranI_permute_mat <- data.frame()
	for(i in 1:permuteNo)
	{	
		randomOrder <- sample(1:N)
		
		XW <- as.matrix(mat[,randomOrder] %*% W)
  		MoranI_permute_mat[rownames(mat),i] <- rowSums(XW * mat[,randomOrder])
	}
	
	XW <- as.matrix(mat %*% W)
	MoranI_permute_mat[rownames(mat), permuteNo+1] <- rowSums(XW * mat)
	
	MoranI_permute_mat <- MoranI_permute_mat/sum(W)
	
	p.Moran_I <- MoranI_permute_mat[,1+permuteNo]
	p.Moran_Z <- apply(MoranI_permute_mat, 1, function(x) (x[permuteNo+1]-mean(x[1:permuteNo]) ) / sd(x[1:permuteNo]) )
	p.Moran_P <- apply(MoranI_permute_mat, 1, function(x) (sum(x[1:permuteNo]>=x[permuteNo+1])+1) / (permuteNo+1) )
	p.Moran_Padj <- p.adjust(p.Moran_P, method="BH")
	
	data.frame(p.Moran_I, p.Moran_Z, p.Moran_P, p.Moran_Padj)
}



spatialAutoCorrelation.Analytical <- function(x, weight)
{
	x <- x[rownames(weight)]
	
	N <- length(x)
	W <- sum(weight)
	
	dx <- x - mean(x, na.rm=TRUE)
	varx <- sum(dx^2, na.rm=TRUE)
	stdx <- sqrt(varx/N)
	zx <- dx/stdx
	
	cv <- zx %o% zx
	cv <- sum(weight*cv, na.rm=TRUE)
	
	Moran_I <- cv/W
	
	# calculate EI
	EI <- (-1) / (N-1)
	
	# calculate VarI
	N2 <- N^2
	W2 <- W^2
	S1 <- sum( (weight + t(weight))^2 ) / 2
	S2 <- sum( (apply(weight, 1, sum) + apply(weight, 2, sum))^2 )
	S3 <- (sum(dx^4)/N) / (varx/N)^2
	S4 <- (N2-3*N+3)*S1 - N*S2 + 3*W2
	S5 <- (N2-N)*S1 - 2*N*S2 + 6*W2
	VarI <- (N*S4 - S3*S5) / ((N-1)*(N-2)*(N-3)*W2) - (EI)^2	
	
	SD <- sqrt(VarI)
	
	Moran_Z <- (Moran_I-EI) / SD	
	Moran_P <- pnorm(Moran_I, mean = EI, sd = SD, lower.tail=FALSE)
  	
  	c(morans.i=Moran_I, z.score=Moran_Z, p.value=Moran_P)
}

run_moran_listw <- function(W)
{
	listw <- spdep::mat2listw(W, row.names = NULL, style="M")
	listw
}

run_moran <- function(mat,listw)
{
	moran_bv.res <- spdep::moran_bv(
  		mat[1,], 
  		mat[2,], 
  		listw, nsim = 2)

	s <- moran_bv.res $ t0
	#p <- (sum(moran_bv.res $ t >= moran_bv.res $ t0)+1)/1001
	
	#c(s,p)
	s
}

run_moran2 <- function(i,fixedj,mat,listw)
{
	moran_bv.res <- spdep::moran_bv(
  		mat[i,], 
  		mat[fixedj,], 
  		listw, nsim = 2)

	s <- moran_bv.res $ t0
	s
}

run_moran_permute <- function(i,mat,listw)
{
	mat_temp <- mat
	colnames(mat_temp) <- sample(colnames(mat),ncol(mat))
	mat_temp <- mat_temp[,colnames(mat)]
	
	mvalue <- run_moran(mat_temp,listw)
					
	mvalue
}

run_moran_permute2 <- function(i,LRdb_rand.m.filter,mat,listw)
{
	mat_temp <- mat[c(LRdb_rand.m.filter[i,1],LRdb_rand.m.filter[i,2]),]
	
	mvalue <- run_moran(mat_temp,listw)
					
	mvalue
}



read.Xena <- function(cancer, dataPath=TCGAPath)
{
	if(cancer%in%c("GBM","LGG","SKCM"))
	{
		tpm1 <- as.matrix(read.csv(gzfile(paste0(dataPath,cancer,".Primary.TPM.gz")),sep="\t",header=T,row.names=1))
		
		if(cancer%in%c("GBM","LGG"))
		{
			tpm2 <- as.matrix(read.csv(gzfile(paste0(dataPath,cancer,".Recurrent.TPM.gz")),sep="\t",header=T,row.names=1))
		}else{
			tpm2 <- as.matrix(read.csv(gzfile(paste0(dataPath,cancer,".Metastatic.TPM.gz")),sep="\t",header=T,row.names=1))
		}
		
		tpm2 <- tpm2[,!grepl(".Normal",colnames(tpm2),fixed=T)]
		tpm <- merge(tpm1,tpm2,by="row.names",all.x=TRUE)
		rownames(tpm) <- tpm[,1]
		tpm <- as.matrix(tpm[,-1])
		tpm[is.na(tpm)] <- 0		
	}else if(cancer=="LAML"){
		tpm <- as.matrix(read.csv(gzfile(paste0(dataPath,cancer,".Blood.TPM.gz")),sep="\t",header=T,row.names=1))
	}else{
		tpm <- as.matrix(read.csv(gzfile(paste0(dataPath,cancer,".Primary.TPM.gz")),sep="\t",header=T,row.names=1))
	}
	
	colnames(tpm) <- gsub(".Normal","",colnames(tpm),fixed = TRUE)
	
	tpm
}


#firebrowse
read.gene.exp <- function(cancer)
{
	gene_exp_path <- paste0("/data/rub2/data/TCGA_firebrowse/mRNAseq_gene_exp_symbol/")

	gene_exp_link <- paste0(gene_exp_path,cancer)
	gene_exp <- read.csv(gene_exp_link,as.is=T,row.names=1,sep="\t")
	
	gene_exp
}

extract.samples <- function(data_matrix,type)
{
	sample_type <- as.numeric(substr(colnames(data_matrix),14,15))
	if(type=="T")
	{
		TS <- sample_type>=1&sample_type<=9
		new_matrix <- data_matrix[,TS]
	}
	if(type=="N")
	{
		NS <- sample_type>=10&sample_type<=19
		new_matrix <- data_matrix[,NS,drop=F]
	}
	return(new_matrix)
}

extract.paired.samples <- function(data_matrix)
{
	sample_type <- as.numeric(substr(colnames(data_matrix),14,15))
	TS_IDs <- colnames(data_matrix)[sample_type==1]
	NS_IDs <- colnames(data_matrix)[sample_type==11]
	TS_matched_IDs <- TS_IDs[substr(TS_IDs,1,12)%in%substr(NS_IDs,1,12)]
	NS_matched_IDs <- NS_IDs[substr(NS_IDs,1,12)%in%substr(TS_IDs,1,12)]
	data_matrix <- data_matrix[,c(NS_matched_IDs,TS_matched_IDs)]
	return(data_matrix)
}

filterMat <- function(mat,cutoff=0.05)
{
	top5p <- round(ncol(mat)*cutoff)
	mat <- mat[rowSums(mat>0)>top5p,]
	mat
}

filter.counts <- function(counts, sample_percent_cutoff = 0.1)
{
	sample_number_cutoff <- round( ncol(counts)*sample_percent_cutoff )
	keep <- rowSums( counts > 1 ) >= sample_number_cutoff
	counts <- counts[keep,]
	return(counts)
}

cor.act.exp <- function(X,Y)
{
	olp <- intersect(rownames(X),rownames(Y))
	X_olp <- X[olp,,drop=F]
	Y_olp <- Y[olp,,drop=F]

	cc_corr_r <- WGCNA::cor(X_olp,Y_olp,use="pairwise.complete.obs")	
	
	olp <- intersect(rownames(cc_corr_r),rownames(Y_olp))
	
	mat1 <- cc_corr_r[olp,]
	mat2 <- Y_olp[olp,]
	
	smyMat <- mapply(function(x, y) cor(x, y), 
                    split(mat1, row(mat1)), 
                    split(mat2, row(mat2)))
	names(smyMat) <- olp

	smyMat
}

cor.two.matrix.same.rows <- function(X,Y)
{
	olp <- intersect(rownames(X),rownames(Y))
	X_olp <- X[olp,,drop=F]
	Y_olp <- Y[olp,,drop=F]
	
	smyMat <- mapply(function(x, y) cor(x, y), 
                    split(X_olp, row(X_olp)), 
                    split(Y_olp, row(Y_olp)))
	names(smyMat) <- olp

	smyMat
}

expand_rows <- function(mat) {
  new_rows <- lapply(1:nrow(mat), function(i) {
    names <- strsplit(rownames(mat)[i], "\\|")[[1]]
    do.call(rbind, replicate(length(names), mat[i, , drop = FALSE], simplify = FALSE)) |>
      `rownames<-`(names)
  })
  do.call(rbind, new_rows)
}

two_col_heatmap_with_arrow_legends <- function(v1, v2,
                                               left_title = "Expression",
                                               right_title = "Activity",
                                               mid_title = "Cell Type",
                                               title = NULL,
                                               gap_x = 3,
                                               barheight_pt = 80,
                                               row_label_size = 3.2,
                                               row_label_color_other = "grey20",
                                               row_label_color_highlight = "black",
                                               band_pad = 0.10,
                                               band_fill = "white",
                                               band_alpha = 1,
                                               arrow_color = "black",
                                               arrow_bg_color = "white",
                                               arrow_lwd = 1.1,
                                               arrow_bg_lwd = 3.5,
                                               arrow_head_len_pt = 8,
                                               label_inset = 0) {
  stopifnot(length(v1) == length(v2))
  n <- length(v1)
  if (is.null(names(v1)) || is.null(names(v2))) {
    names(v1) <- names(v2) <- seq_len(n)
  }
  stopifnot(identical(names(v1), names(v2)))
  lab <- names(v1)

  lims1 <- range(v1, na.rm = TRUE); if (diff(lims1) == 0) lims1 <- lims1 + c(-0.5, 0.5)
  lims2 <- range(v2, na.rm = TRUE); if (diff(lims2) == 0) lims2 <- lims2 + c(-0.5, 0.5)

  df1 <- data.frame(x = 1, index = seq_len(n), value = v1, label = lab)
  df2 <- data.frame(x = gap_x, index = seq_len(n), value = v2, label = lab)

  y1 <- which.max(v1); y2 <- which.max(v2)

  half_w <- 0.45
  arrow_df <- data.frame(x = 1 + half_w, y = y1,
                         xend = gap_x - half_w, yend = y2)

  band_xmin <- 1 + half_w + band_pad
  band_xmax <- gap_x - half_w - band_pad
  if (band_xmax <= band_xmin) { band_xmin <- 1 + 0.55; band_xmax <- gap_x - 0.55 }
  band_df <- data.frame(xmin = band_xmin, xmax = band_xmax,
                        ymin = 0.5, ymax = n + 0.5)

  mid_x <- (1 + gap_x) / 2
  base_labels <- data.frame(x = mid_x, index = seq_len(n), label = lab)
  df_other <- base_labels[!(base_labels$index %in% c(y1, y2)), ]
  df_start <- base_labels[base_labels$index == y1, , drop = FALSE]
  df_end   <- base_labels[base_labels$index == y2, , drop = FALSE]
  #df_start$x <- band_xmin + label_inset
  #df_end$x   <- band_xmax - label_inset

  library(ggnewscale)
  p_main <- ggplot() +
    # Heatmaps
    geom_tile(data = df1, aes(x, index, fill = value),
              width = 0.9, height = 0.9, color = "white") +
    scale_fill_gradient(name = left_title,
                        low = "grey80", high = "red", limits = lims1, guide = "none") +
    new_scale_fill() +
    geom_tile(data = df2, aes(x, index, fill = value),
              width = 0.9, height = 0.9, color = "white") +
    scale_fill_gradient(name = right_title,
                        low = "grey80", high = "blue", limits = lims2, guide = "none") +
    # Outline maxima
    geom_tile(data = subset(df1, index == y1), aes(x, index),
              fill = NA, color = "black") +
    geom_tile(data = subset(df2, index == y2), aes(x, index),
              fill = NA, color = "black") +
    # Middle band
    geom_rect(data = band_df,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              inherit.aes = FALSE, fill = band_fill,
              alpha = band_alpha, colour = NA) +
    # Labels
    geom_text(data = df_other, aes(x = x, y = index, label = label),
              hjust = 0.5, size = row_label_size,
              color = row_label_color_other) +
    geom_text(data = df_start, aes(x = x, y = index, label = label),
              hjust = 0.5, size = row_label_size, fontface = "bold",
              color = row_label_color_highlight) +
    geom_text(data = df_end, aes(x = x, y = index, label = label),
              hjust = 0.5, size = row_label_size, fontface = "bold",
              color = row_label_color_highlight) +
    # Arrow
    geom_segment(data = arrow_df, aes(x, y, xend = xend, yend = yend),
                 colour = arrow_bg_color, linewidth = arrow_bg_lwd,
                 lineend = "round", inherit.aes = FALSE) +
    geom_segment(data = arrow_df, aes(x, y, xend = xend, yend = yend),
                 colour = arrow_color, linewidth = arrow_lwd,
                 lineend = "round",
                 arrow = arrow(type = "closed",
                               length = unit(arrow_head_len_pt, "pt")),
                 inherit.aes = FALSE) +
    # Axis on top
    scale_x_continuous(breaks = c(1, mid_x, gap_x),
                       labels = c(left_title, mid_title, right_title),
                       position = "top") +
    expand_limits(x = c(0.3, gap_x + 0.7)) +
    coord_cartesian(expand = FALSE, clip = "off") +
    labs(x = NULL, y = NULL, title = title) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.text = element_text(colour = "black"),
			axis.title = element_text(colour = "black"),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none")

  # Legends
  p_leg_left <- ggplot(df1, aes(1, 1, fill = value)) +
    geom_raster() +
    scale_fill_gradient(name = "Exp.",
                        low = "grey80", high = "red", limits = lims1,
                        guide = guide_colorbar(title.position = "top", title.hjust = 0,
                                               direction = "vertical",
                                               barheight = unit(barheight_pt, "pt"),
                                               label.position = "right",
                                               label.theme = element_text(hjust = 1))) +
    theme_void()

  p_leg_right <- ggplot(df2, aes(1, 1, fill = value)) +
    geom_raster() +
    scale_fill_gradient(name = "Act.",
                        low = "grey80", high = "blue", limits = lims2,
                        guide = guide_colorbar(title.position = "top", title.hjust = 1,
                                               direction = "vertical",
                                               barheight = unit(barheight_pt, "pt"),
                                               label.position = "left")) +
    theme_void()

  leg_left  <- cowplot::get_legend(p_leg_left)
  leg_right <- cowplot::get_legend(p_leg_right)

  cowplot::plot_grid(leg_left, p_main, leg_right,
                     nrow = 1, align = "h",
                     rel_widths = c(0.22, 1, 0.22))
}



transferMouseToHuman <- function(mat)
{
	mgi_ann = read.table("/data/rub2/data/HOM_MouseHumanSequence.rpt.gz", sep = "\t", header = TRUE)
	mgi_ann = mgi_ann[, c(1,2,4)]
	colnames(mgi_ann) = c("homoloid", "org", "symbol")
	mgi_ann$org = gsub(", laboratory", "", mgi_ann$org)
	
	mouseGenes <- rownames(mat)
	for(i in 1:length(mouseGenes))
	{
		id <- mgi_ann[mgi_ann[,3]==mouseGenes[i],1]
		if(length(id)==0)
		{
			mouseGenes[i] <- "humanGeneNotExist"
		}else{
			hGenes <- mgi_ann[mgi_ann[,1]%in%id&mgi_ann[,2]=="human",3]
			if(length(hGenes)==0)
			{
				mouseGenes[i] <- "humanGeneNotExist"
			}else{
				if(length(hGenes)==1)
				{
					mouseGenes[i] <- hGenes
				}else{
					if(toupper(mouseGenes[i])%in%hGenes)
					{
						mouseGenes[i] <- toupper(mouseGenes[i])
					}else{
						mouseGenes[i] <- "humanGeneMultiple"
					}
				}
			}
		}
	}
	rownames(mat) <- mouseGenes
	mat <- mat[!rownames(mat)%in%c("humanGeneNotExist","humanGeneMultiple"),,drop=F]
	
	mat
}

transferSymbolFromMouseToHuman <- function(mouseGenes)
{
	mgi_ann = read.table("/data/rub2/data/HOM_MouseHumanSequence.rpt.gz", sep = "\t", header = TRUE)
	mgi_ann = mgi_ann[, c(1,2,4)]
	colnames(mgi_ann) = c("homoloid", "org", "symbol")
	mgi_ann$org = gsub(", laboratory", "", mgi_ann$org)
	
	for(i in 1:length(mouseGenes))
	{
		id <- mgi_ann[mgi_ann[,3]==mouseGenes[i],1]
		if(length(id)==0)
		{
			mouseGenes[i] <- "humanGeneNotExist"
		}else{
			hGenes <- mgi_ann[mgi_ann[,1]%in%id&mgi_ann[,2]=="human",3]
			if(length(hGenes)==0)
			{
				mouseGenes[i] <- "humanGeneNotExist"
			}else{
				if(length(hGenes)==1)
				{
					mouseGenes[i] <- hGenes
				}else{
					if(toupper(mouseGenes[i])%in%hGenes)
					{
						mouseGenes[i] <- toupper(mouseGenes[i])
					}else{
						mouseGenes[i] <- "humanGeneMultiple"
					}
				}
			}
		}
	}
	
	mouseGenes
}

extract_experimental_varified_LR <- function(dataPath)
{
	LRdb <- data.frame()
	
	for(f in c("Cytokine","Chemokine","Growth_Factor"))
	{
		data_list <- rio::import_list(paste0(dataPath,"Receptor_",f,".xlsx"))
	
		for(i in names(data_list))
		{
			temp <- data_list[[i]]
			
			for(j in 1:nrow(temp))
			{	
				items_1 <- unlist( strsplit(temp[j,2],"+",fixed=T) )
				items_1 <- unlist( strsplit(items_1,",",fixed=T) )
				items_1 <- gsub(" ","",items_1,fixed=T)
				items_1 <- gsub("-","",items_1,fixed=T)
				items_1 <- unique(items_1)
				
				items_2 <- unlist( strsplit(temp[j,3],"+",fixed=T) )
				items_2 <- unlist( strsplit(items_2,",",fixed=T) )
				items_2 <- gsub(" ","",items_2,fixed=T)
				items_2 <- gsub("-","",items_2,fixed=T)
				items_2 <- unique(items_2)
				
				items12 <- as.matrix(expand.grid(items_1,items_2))
				
				LRdb[paste0(items12[,1],"_",items12[,2]),"Ligand"] <- items12[,1]
				LRdb[paste0(items12[,1],"_",items12[,2]),"Receptor"] <- items12[,2]
			}
		}
	}
	
	LRdb
}





CoxPH_best_separation <- function(X, Y, margin)
{
  # part 1: continuous regression
  errflag = F
  
  coxph.fit = tryCatch(
    coxph(Y~., data=X),
    error = function(e) errflag <<- T,
    warning = function(w) errflag <<- T)

  if(errflag) return (NA)
  
  n_r = nrow(X)
  n_c = ncol(X)
  
  arr_result = summary(coxph.fit)$coef[n_c,]
  
  if(is.na(arr_result["z"])) return (NA)
  
  # no need to find optimal threshold
  if(is.null(margin)) return (arr_result)
  
  # part 2: find the optimal threshold
  arr = X[, n_c]
  vthres = sort(arr)
  
  # these are missing values, not NULL not existing values
  zscore_opt = thres_opt = NA
  
  for(i in (margin+1):(n_r-margin))
  {
    X[, n_c] = as.numeric(arr >= vthres[i])
    
    errflag = F
    coxph.fit = tryCatch(
      coxph(Y~., data=X),
      error = function(e) errflag <<- T,
      warning = function(w) errflag <<- T)
    
    if(errflag) next
    
    z = summary(coxph.fit)$coef[n_c, "z"]
    if(is.na(z)) next
    
    if (is.na(zscore_opt)){
      zscore_opt = z
      thres_opt= vthres[i]
    
    }else if(arr_result['z'] > 0){
      if(z > zscore_opt){
        zscore_opt = z
        thres_opt = vthres[i]
      }
      
    }else{ # arr_result['z'] <= 0
      if(z < zscore_opt){
        zscore_opt = z
        thres_opt = vthres[i]
      }
    }
  }
  
  arr_result['thres.opt'] = thres_opt
  arr_result['z.opt'] = zscore_opt
  
  return (arr_result)
}

run_CoxPH_best_separation <- function(data,survival,margin)
{
	# align matrix names
	common = Reduce(intersect, list(rownames(data),rownames(survival)))
	sprintf("%s samples", length(common))
	
	data = data[common,,drop=F]
	survival = survival[common,,drop=F]
	
	# stop at low death rate
	death_rate = sum(survival[,2])/dim(survival)[1]
	if(length(death_rate) < 0.1) q()
	
	# split up survival and background
	surv = Surv(survival[,1], survival[,2])
	
	if(dim(survival)[2] > 2){
	  B = survival[,3:dim(survival)[2], drop=F]
	}else{
	  B = survival[,c(), drop=F]
	}
	
	# build up regression data space
	B = cbind(B, rep(0, dim(data)[1]))
	B = as.data.frame(B)
	N_B = ncol(B)
	colnames(B)[N_B] = "pivot"
	
	# iterate over features
	features = colnames(data)
	N = length(features)
	
	result = NULL
	
	step = round(max(N/100,1))
	for (i in 1:N)
	{
	  # progress report
	  if((i %% step) == 0){
	    sprintf("%s", round(100 * i/N, 2))
	  }
	  
	  fid = features[i]
	  #if(!(tail(strsplit(fid, '@')[[1]], n=1) %in% c('DDX3Y', 'FIBP', 'FCMR'))) next
	  
	  # part 1: overall regression
	  arr = B[,N_B] = data[,i]
	  
	  arr_result = CoxPH_best_separation(B, surv, margin)
	  
	  if(sum(is.na(arr_result)) > 0){
	    warning(paste0('Jump with failed continuous regression ', fid))
	    next
	  }
	  
	  if(is.null(result)){
	    result = matrix(nrow = N, ncol=length(arr_result))
	    colnames(result) = names(arr_result)
	    rownames(result) = features
	  }
	  
	  result[i,] = arr_result
	}
	
	mean.value = colMeans(data)
	result = cbind(result, mean.value)
	
	N = rep(length(common), dim(result)[1])
	result = cbind(result, N)
	
	list(data,survival,result)
}

