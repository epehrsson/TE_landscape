# PCA on binary matrix of TEs in enhancer states x sample
# See 11/3/2016, 11/4/2016, 2/6/2017, 2/9/2017

# Matrix of TEs
enhancer_matrix = read.table("chromHMM/enhancer_matrix.txt",sep='\t')
colnames(enhancer_matrix) = c(TE_coordinates[c(1:4,6,5,7)],"Sample")
enhancer_matrix = dcast(enhancer_matrix,chromosome+start+stop+subfamily+family+class+strand~Sample,length)
rownames(enhancer_matrix) = apply(enhancer_matrix,1,function(x) paste(x[1],x[2],x[3],x[4],x[5],x[6],x[7],sep="_"))
enhancer_matrix = enhancer_matrix[which(enhancer_matrix$chromosome != "chrY"),8:134]
enhancer_matrix[enhancer_matrix > 0] = 1

# Compute PCA
enhancer_matrix_pca = prcomp(t(enhancer_matrix),scale=TRUE,center=TRUE)