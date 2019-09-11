# Creates a matrix of the methylation level of each hg19 TE (row) in each sample (column) ("TE_meth_average")
# For all samples with WGBS data (n=37), including only TEs that overlap a CpG

# Load average methylation level for all TEs
TE_meth_average = read.table("WGBS/TE_CpG_Meth_new_average.txt",sep='\t',header=TRUE)

# Load number of CpGs per TE
TE_CpG_count = read.table("WGBS/TE_CpG_count.txt",sep='\t',col.names=c(TE_coordinates[c(1:4,6,5,7)],"CpGs"))
TE_CpG_count$CpGs = TE_CpG_count$CpGs/2

# Remove TEs that do not overlap any CpG
TE_meth_average = merge(TE_meth_average,TE_CpG_count,by=TE_coordinates[c(1:4,6,5,7)])
rm(TE_CpG_count)

# Number of samples each TE is in each methylation state
# With and without IMR90 (E017)
TE_meth_average$Hypomethylated = apply(TE_meth_average,1,function(x) sum(as.numeric(x[8:44]) < 0.3,na.rm=TRUE))
TE_meth_average$Hypermethylated = apply(TE_meth_average,1,function(x) sum(as.numeric(x[8:44]) > 0.7,na.rm=TRUE))
TE_meth_average$Intermediate = apply(TE_meth_average,1,function(x) sum(as.numeric(x[8:44]) <= 0.7 & as.numeric(x[8:44]) >= 0.3,na.rm=TRUE))
TE_meth_average$Missing = apply(TE_meth_average,1,function(x) sum(is.na(as.numeric(x[8:44]))))
TE_meth_average$Hypomethylated_noIMR90 = apply(TE_meth_average,1,function(x) sum(as.numeric(x[c(8:17,19:44)]) < 0.3,na.rm=TRUE))
TE_meth_average$Hypermethylated_noIMR90 = apply(TE_meth_average,1,function(x) sum(as.numeric(x[c(8:17,19:44)]) > 0.7,na.rm=TRUE))
TE_meth_average$Intermediate_noIMR90 = apply(TE_meth_average,1,function(x) sum(as.numeric(x[c(8:17,19:44)]) <= 0.7 & as.numeric(x[c(8:17,19:44)]) >= 0.3,na.rm=TRUE))
TE_meth_average$Missing_noIMR90 = apply(TE_meth_average,1,function(x) sum(is.na(as.numeric(x[c(8:17,19:44)]))))

# Update class assignments
TE_meth_average$class_update = convert_class(TE_meth_average$class)

# Number of unique states each TE is annotated with across all samples
# With and without IMR90
TE_meth_average$States = apply(TE_meth_average[,46:49],1,function(x) sum(x > 0))
TE_meth_average$States_noIMR90 = apply(TE_meth_average[,50:53],1,function(x) sum(x > 0))

# Write out matrix
write.table(melt(TE_meth_average[,1:44],id.vars = TE_coordinates),file="WGBS/TE_WGBS_state.txt",sep='\t',quote=FALSE,row.names=FALSE)

# Save dataframe
save(TE_meth_average,file="R_datasets/TE_meth_average.RData")