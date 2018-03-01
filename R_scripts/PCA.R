# PCA on binary matrices of TEs overlapping DNase or H3K27ac by sample

#load("R_datasets/TE_DNase_peaks.RData")
#load("R_datasets/TE_H3K27ac_peaks.RData")

# Configure matrices
TE_DNase_binary = TE_DNase_peaks[,1:60]
TE_DNase_binary[,8:60] = ifelse(TE_DNase_binary[,8:60]>0, 1, 0)
TE_DNase_binary[which(TE_DNase_binary$chromosome == "chrY"),as.vector(metadata[which(metadata$chrY == "No"),]$Sample)] = NA
TE_DNase_binary = TE_DNase_binary[which(apply(TE_DNase_binary[,8:60],1,var) != 0),]
rownames(TE_DNase_binary) = apply(TE_DNase_binary,1,function(x) paste(x[1],x[2],x[3],x[4],x[5],x[6],x[7],sep="_"))

TE_H3K27ac_binary = TE_H3K27ac_peaks[,1:105]
TE_H3K27ac_binary[,8:105] = ifelse(TE_H3K27ac_binary[,8:105]>0, 1, 0)
TE_H3K27ac_binary[which(TE_H3K27ac_binary$chromosome == "chrY"),as.vector(metadata[which(metadata$chrY == "No"),]$Sample)] = NA
TE_H3K27ac_binary = TE_H3K27ac_binary[which(apply(TE_H3K27ac_binary[,8:105],1,var) != 0),]
rownames(TE_H3K27ac_binary) = apply(TE_H3K27ac_binary,1,function(x) paste(x[1],x[2],x[3],x[4],x[5],x[6],x[7],sep="_"))

# Compute PCA
DNase_pca = prcomp(t(na.omit(TE_DNase_binary[,8:60])),scale=TRUE,center=TRUE)
H3K27ac_pca = prcomp(t(na.omit(TE_H3K27ac_binary[,8:105])),scale=TRUE,center=TRUE)
