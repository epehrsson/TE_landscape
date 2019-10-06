# Calculates density for CpGs overlapping TEs
# 9/20/2016

cpg = read.table("~/TE_landscape/WGBS/CpG_TE_Meth.bed",sep='\t')
samples = as.vector(read.table("~/TE_landscape/sample_lists/WGBS_samples.txt")$V1)
colnames(cpg) = c("chr","start","stop",samples)
cpg = melt(cpg[,4:40])
colnames(cpg) = c("Sample","Methylation")
cpg = cpg[which(cpg$Methylation != -1),]
save(cpg,"~/TE_landscape/R_datasets/cpg_density.RData")

# Old code: Global methylation in bins of 0.01
#TE_landscape/WGBS/all_CpG_Meth_distribution.txt
#python ../bin/TE_landscape/methylation_distribution_multi.py all_CpG_Meth.bed all_CpG_Meth_distribution.txt WGBS_samples.txt