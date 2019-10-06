# Global distribution of CpG methylation values
# See 9/20/2016
# NEEDS TO BE UPDATED

all_cpg_meth_distribution = read.table("all_CpG_Meth_distribution.txt",sep='\t',header=TRUE,row.names=1)
all_cpg_meth_distribution$Bins = c(0,rep(seq(0.05,1,0.05),times=rep(5,20)))
all_cpg_meth_distribution = aggregate(data=all_cpg_meth_distribution,.~Bins,FUN=sum)
rownames(all_cpg_meth_distribution) = all_cpg_meth_distribution$Bins
all_cpg_meth_distribution = all_cpg_meth_distribution[,2:38]
all_cpg_meth_distribution = all_cpg_meth_distribution/colSums(all_cpg_meth_distribution)
all_cpg_meth_distribution_long = melt(as.matrix(all_cpg_meth_distribution))
colnames(all_cpg_meth_distribution_long) = c("Bin","Sample","CpGs")
