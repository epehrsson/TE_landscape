# Average methylation level across CpGs by sample and feature
# See 5/25/2016, 9/16/2016, 12/15/2016, 5/11/2017, 5/12/2017, 7/21/2017, 8/7/2017

library(reshape2)
library(plyr)

# Import average methylation of CpGs (all, in TEs, not in TEs)
CpG_meth_average = rbind(read.table("WGBS/all_CpG_Meth_average_new.txt",header=TRUE),read.table("WGBS/CpG_noTE_Meth_average.txt",header=TRUE),read.table("WGBS/CpG_TE_Meth_average.txt",header=TRUE))
rownames(CpG_meth_average) = c("All","noTE","TE")
CpG_meth_average = melt(as.matrix(CpG_meth_average))
colnames(CpG_meth_average) = c("Cohort","Sample","Methylation")
CpG_meth_average = CpG_meth_average[order(CpG_meth_average$Cohort,CpG_meth_average$Methylation),]
CpG_meth_average$Sample = factor(CpG_meth_average$Sample,levels=CpG_meth_average[which(CpG_meth_average$Cohort == "All"),]$Sample)

# Import average methylation of CpGs in each Refseq genic feature
CpG_feature_meth_average = ldply(list.files(path = "WGBS/Refseq_features/", pattern = "average.txt", full.names = TRUE),function(x) read.table(x,header=TRUE))
feature_names = list.files(path = "WGBS/Refseq_features/", pattern = "average.txt")
feature_names = gsub("CpG_refseq_","",feature_names)
feature_names = gsub("_merge_noTE_Meth_average.txt","",feature_names)
feature_names = gsub("_noTE_Meth_average.txt","",feature_names)
rownames(CpG_feature_meth_average) = feature_names
CpG_feature_meth_average = melt(as.matrix(CpG_feature_meth_average))
colnames(CpG_feature_meth_average) = c("Cohort","Sample","Methylation")

# Compare all/TE CpGs to other features
CpG_feature_meth_average = rbind(droplevels(CpG_meth_average[which(CpG_meth_average$Cohort != "noTE"),]),CpG_feature_meth_average)
CpG_feature_meth_average = CpG_feature_meth_average[order(CpG_feature_meth_average$Cohort,CpG_feature_meth_average$Methylation),]