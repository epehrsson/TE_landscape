cpg = read.table("~/TE_landscape/WGBS/CpG_TE_Meth.bed",sep='\t')
samples = as.vector(read.table("~/TE_landscape/sample_lists/WGBS_samples.txt")$V1)
colnames(cpg) = c("chr","start","stop",samples)
cpg = melt(cpg[,4:40])
colnames(cpg) = c("Sample","Methylation")
cpg = cpg[which(cpg$Methylation != -1),]
save(cpg,"~/TE_landscape/R_scripts/cpg_density.RData")