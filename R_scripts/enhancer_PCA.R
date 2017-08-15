# PCA on binary matrix of TEs in enhancer states x sample
# See 11/3/2016, 11/4/2016, 2/6/2017, 2/9/2017

# Matrix of TEs
enhancer_matrix = read.table("enhancer_matrix.txt",sep='\t',header=TRUE)

# Matrix of other TEs
enhancer_matrix_other = read.table("other/enhancer_matrix_other.txt",sep='\t',header=TRUE)

# Combine
enhancer_matrix = rbind(enhancer_matrix,enhancer_matrix_other)
rm(enhancer_matrix_other)
enhancer_matrix = enhancer_matrix[which(apply(enhancer_matrix[,8:134],1,function(x) sum(x > 0)) > 0),]
enhancer_matrix[which(enhancer_matrix$chromosome == "chrY"),which(colnames(enhancer_matrix) %in% c("E116", "E117", "E123", "E124", "E126", "E127"))] = NA
rownames(enhancer_matrix) = apply(enhancer_matrix,1,function(x) paste(x[1],x[2],x[3],x[4],x[5],x[6],x[7],sep="_"))

# Compute PCA
enhancer_matrix_TEother_pca = prcomp(t(na.omit(enhancer_matrix[,8:134])),scale=TRUE,center=TRUE)

# Additional analyses
# enhancer_matrix_dist = cor(enhancer_matrix[,8:134],use="pairwise.complete.obs")
# enhancer_matrix_hclust = hclust(as.dist(enhancer_matrix_dist),method="complete")
# plot(enhancer_matrix_hclust,hang=-1)
# summary(enhancer_matrix_pca) #To get the proportion and cumulative proportion of variance explained by each PC
