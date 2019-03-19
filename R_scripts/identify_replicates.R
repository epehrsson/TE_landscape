library(plyr)
library(reshape2)
library(ggplot2)

# Load matrices
load("raw_data/correlation_matrices/cor_H3K27me3_13.RData")
h3k27me3 = markcor

load("raw_data/correlation_matrices/cor_H3K36me3_4.RData")
h3k36me3 = markcor

load("raw_data/correlation_matrices/cor_H3K4me1_7.RData")
h3k4me1 = markcor

load("raw_data/correlation_matrices/cor_H3K4me3_1.RData")
h3k4me3 = markcor

load("raw_data/correlation_matrices/cor_H3K9me3_9.RData")
h3k9me3 = markcor

transform_matrix = function(matrix){
  matrix[lower.tri(matrix)] = NA
  matrix = melt(as.matrix(matrix))
  colnames(matrix) = c("Sample 1","Sample 2","Correlation")
  matrix = matrix[which(matrix$`Sample 1` != matrix$`Sample 2` & !is.na(matrix$Correlation)),]
  matrix = matrix[order(matrix$Correlation,decreasing = TRUE),]
  matrix$Rank = 1:nrow(matrix)
  return(matrix)
}

marks = ls(pattern="h3")
mark_cor = lapply(marks,function(x) transform_matrix(get(x)))
names(mark_cor) = marks
mark_cor = ldply(mark_cor)
colnames(mark_cor)[1] = "Mark"

ggplot(mark_cor,aes(x=Correlation)) + geom_histogram() + facet_wrap(~Mark)

mark_cor_ranked = ddply(mark_cor,.(`Sample 1`,`Sample 2`),summarise,Rank=mean(Rank))
mark_cor_ranked = mark_cor_ranked[order(mark_cor_ranked$Rank),]

replicates = mark_cor_ranked[c(1:2,5:7,13:15,17,19,21,23:24,31,33,38,40,46:47,59),]
rm(list=c("mark_cor_ranked","mark_cor","markcor","h3k27me3","h3k36me3","h3k4me1","h3k4me3","h3k9me3","marks","transform_matrix"))
