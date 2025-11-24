source("Secretome_s0_path.R")

library(rio) 
xlsx <- import_list(paste0(dataValPath,"/blocking/IL11_Widjaja2024_TableS5.xlsx")) 

s0 <- as.data.frame(xlsx[[1]]) 
s0 <- s0[s0[,2]>0,]
s0[,3] <- as.numeric(s0[,3])
s0 <- s0[order(s0[,2],decreasing=T),]
s0 <- s0[!duplicated(s0[,8]),]
rownames(s0) <- s0[,8]
s1 <- s0[,3,drop=F]

s0 <- as.data.frame(xlsx[[2]]) 
s0 <- s0[s0[,2]>0,]
s0[,3] <- as.numeric(s0[,3])
s0 <- s0[order(s0[,2],decreasing=T),]
s0 <- s0[!duplicated(s0[,8]),]
rownames(s0) <- s0[,8]
s2 <- s0[,3,drop=F]

s0 <- as.data.frame(xlsx[[3]]) 
s0 <- s0[s0[,2]>0,]
s0[,3] <- as.numeric(s0[,3])
s0 <- s0[order(s0[,2],decreasing=T),]
s0 <- s0[!duplicated(s0[,8]),]
rownames(s0) <- s0[,8]
s3 <- s0[,3,drop=F]

olp <- intersect(rownames(s1), intersect(rownames(s2),rownames(s3)))

f3 <- cbind(vWat=s1[olp,], cbind(liver=s2[olp,], gastro=s3[olp,]) )
rownames(f3) <- olp

f3 <- transferMouseToHuman(f3)

rownames(f3) <- transferSymbol(rownames(f3))
f3 <- rm_duplicates(f3)

write.table(f3,paste0(dataValPath,"blocking/IL11_Widjaja2024.diff"),quote=FALSE,sep="\t")

