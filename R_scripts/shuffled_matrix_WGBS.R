# Creates matrices of the methylation level of each shuffled hg19 TE (row) in each sample (column) ("shuffled_WGBS_average")
# For 10 iterations of shuffled TEs

# Load matrices of the average methylation level for shuffled TEs
shuffled_WGBS_average = lapply(seq(1,10,1),function(x) read.table(paste("WGBS/shuffled/rmsk_TE_shuffle_",x,"_Meth_average.txt",sep=""),sep='\t',header=TRUE))
print("Loaded WGBS matrices")

# Load the number of CpGs per shuffled TE
shuffled_WGBS_CpG = lapply(seq(1,10,1),function(x) read.table(paste("WGBS/shuffled/TE_CpG_count_",x,".txt",sep=""),sep='\t', col.names=c(TE_coordinates[c(1:4,6,5,7)],"CpGs")))
print("Loaded CpG matrices")

# Remove TEs that do not overlap any CpG
for (i in 1:10){
  shuffled_WGBS_average[[i]] = merge(shuffled_WGBS_average[[i]],shuffled_WGBS_CpG[[i]],by=TE_coordinates[c(1:4,6,5,7)])
  shuffled_WGBS_average[[i]]$CpGs = shuffled_WGBS_average[[i]]$CpGs/2
}
rm(shuffled_WGBS_CpG)
print("Combined WGBS and CpG matrices")

# Number of samples each TE is in each methylation state
shuffled_WGBS_average = lapply(shuffled_WGBS_average,function(y) transform(y,Hypomethylated = apply(y,1,function(x) sum(as.numeric(x[8:44]) < 0.3,na.rm=TRUE)),
                                                                           Hypermethylated = apply(y,1,function(x) sum(as.numeric(x[8:44]) > 0.7,na.rm=TRUE)),
                                                                           Intermediate = apply(y,1,function(x) sum(as.numeric(x[8:44]) <= 0.7 & as.numeric(x[8:44]) >= 0.3,na.rm=TRUE)),
                                                                           Missing = apply(y,1,function(x) sum(is.na(x[8:44])))))
print("Calculating methylation states")

# Number of unique methylation states each TE is annotated with across all samples
shuffled_WGBS_average = lapply(shuffled_WGBS_average,function(x) transform(x,States = apply(x[,46:49],1,function(y) sum(y > 0))))
print("Counting total states")

# Save dataframes
save(shuffled_WGBS_average,file="R_datasets/shuffled_WGBS.RData")