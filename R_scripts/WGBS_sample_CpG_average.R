# Average methylation level across CpGs by sample (all, TEs, not TEs)
# See 5/25/2016, 9/16/2016, 12/15/2016, 5/11/2017, 5/12/2017, 7/21/2017, 8/7/2017

# Import all averages 
CpG_meth_average = rbind(read.table("WGBS/all_CpG_Meth_average_new.txt",header=TRUE),read.table("WGBS/CpG_noTE_Meth_average.txt",header=TRUE),read.table("WGBS/CpG_TE_Meth_average.txt",header=TRUE))
rownames(CpG_meth_average) = c("All","noTE","TE")
CpG_meth_average = melt(as.matrix(CpG_meth_average))
colnames(CpG_meth_average) = c("Cohort","Sample","Methylation")
CpG_meth_average = CpG_meth_average[order(CpG_meth_average$Cohort,CpG_meth_average$Methylation),]
CpG_meth_average$Sample = factor(CpG_meth_average$Sample,levels=CpG_meth_average[which(CpG_meth_average$Cohort == "All"),]$Sample)

# Import average methylation of CpGs in each Refseq genic feature
CpG_feature_meth_average = ldply(list.files(path = "WGBS/Refseq_features/", pattern = "average.txt", full.names = TRUE),function(x) read.table(x,header=TRUE))
rownames(CpG_feature_meth_average) = as.vector(feature_contribution$Cohort[c(1:11,13:19)])
CpG_feature_meth_average = melt(as.matrix(CpG_feature_meth_average))
colnames(CpG_feature_meth_average) = c("Cohort","Sample","Methylation")
CpG_feature_meth_average = CpG_feature_meth_average[order(CpG_feature_meth_average$Cohort,CpG_feature_meth_average$Methylation),]

# Averages across samples
by(CpG_meth_average,CpG_meth_average$Cohort,function(x) mean(x$Methylation))
by(CpG_meth_average,CpG_meth_average$Cohort,function(x) sd(x$Methylation))
by(CpG_meth_average,CpG_meth_average$Cohort,function(x) mean(x$Methylation)-2*sd(x$Methylation))
ddply(CpG_feature_meth_average,"Cohort",summarise,Mean=mean(Methylation),SD=sd(Methylation))