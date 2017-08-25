# PCA on binary matrix of TEs in enhancer states x sample
# See 11/3/2016, 11/4/2016, 2/6/2017, 2/9/2017

# Matrix of TEs
enhancer_matrix = rbind(read.table("chromHMM/enhancer_matrix.txt",sep='\t',header=TRUE),read.table("chromHMM/enhancer_matrix_other.txt",sep='\t',header=TRUE))
enhancer_matrix = enhancer_matrix[which(apply(enhancer_matrix[,8:134],1,function(x) sum(x > 0)) > 0),]
enhancer_matrix[which(enhancer_matrix$chromosome == "chrY"),which(colnames(enhancer_matrix) %in% c("E116", "E117", "E123", "E124", "E126", "E127"))] = NA
rownames(enhancer_matrix) = apply(enhancer_matrix,1,function(x) paste(x[1],x[2],x[3],x[4],x[5],x[6],x[7],sep="_"))

# Compute PCA
enhancer_matrix_pca = prcomp(t(na.omit(enhancer_matrix[,8:134])),scale=TRUE,center=TRUE)