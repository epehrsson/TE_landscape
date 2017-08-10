# Histone modification ChIP-seq, DNase, and methylation over TEs in each chromHMM state
# See 12/13/2016, 3/8/2017, 8/7/2017

histones = read.table("compare_marks/ChIP_histone/TEother_average.txt",sep='\t')
colnames(histones) = c("Bin","Level","Bases","Mark","State")
histones$Bin = factor(histones$Bin,levels=unique(histones$Bin))
histones$Mark = factor(histones$Mark,levels=levels(histones$Mark)[c(5:6,8,3:4,7,2,1,9)])