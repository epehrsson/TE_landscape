# Load matrix of methylation level per TE

# Methylation level for all TEs
TE_meth_average = read.table("WGBS/TE_CpG_Meth_new_average.txt",sep='\t',header=TRUE)

# Number of CpGs per TE
TE_CpG_count = read.table("WGBS/TE_CpG_count.txt",sep='\t')
colnames(TE_CpG_count) = c(TE_coordinates[c(1:4,6,5,7)],"CpGs")
TE_CpG_count$CpGs = TE_CpG_count$CpGs/2

# Methylation level for TEs with at least one CpG
TE_meth_average = merge(TE_meth_average,TE_CpG_count,by=TE_coordinates[c(1:4,6,5,7)])
rm(TE_CpG_count)

# Number of samples TE is in each methylation state (IMR90, no IMR90)
TE_meth_average$Hypomethylated = apply(TE_meth_average,1,function(x) sum(as.numeric(x[8:44]) < 0.3,na.rm=TRUE))
TE_meth_average$Hypermethylated = apply(TE_meth_average,1,function(x) sum(as.numeric(x[8:44]) > 0.7,na.rm=TRUE))
TE_meth_average$Intermediate = apply(TE_meth_average,1,function(x) sum(as.numeric(x[8:44]) <= 0.7 & as.numeric(x[8:44]) >= 0.3,na.rm=TRUE))
TE_meth_average$Missing = apply(TE_meth_average,1,function(x) sum(is.na(as.numeric(x[8:44]))))
TE_meth_average$Hypomethylated_noIMR90 = apply(TE_meth_average,1,function(x) sum(as.numeric(x[c(8:17,19:44)]) < 0.3,na.rm=TRUE))
TE_meth_average$Hypermethylated_noIMR90 = apply(TE_meth_average,1,function(x) sum(as.numeric(x[c(8:17,19:44)]) > 0.7,na.rm=TRUE))
TE_meth_average$Intermediate_noIMR90 = apply(TE_meth_average,1,function(x) sum(as.numeric(x[c(8:17,19:44)]) <= 0.7 & as.numeric(x[c(8:17,19:44)]) >= 0.3,na.rm=TRUE))
TE_meth_average$Missing_noIMR90 = apply(TE_meth_average,1,function(x) sum(is.na(as.numeric(x[c(8:17,19:44)]))))

# Adding Unconfident class
TE_meth_average$class_update = convert_class(TE_meth_average$class)

# Adding number of states
TE_meth_average$States = apply(TE_meth_average[,46:49],1,function(x) sum(x > 0))
TE_meth_average$States_noIMR90 = apply(TE_meth_average[,50:53],1,function(x) sum(x > 0))

# Write out for WGBS state
write.table(melt(TE_meth_average[,1:44],id.vars = TE_coordinates),file="WGBS/TE_WGBS_state.txt",sep='\t',quote=FALSE,row.names=FALSE)

save(TE_meth_average,file="R_datasets/TE_meth_average.RData")
