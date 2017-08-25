# Load shuffled TE potential matrices, WGBS

# Methylation level for shuffled TEs
shuffled_WGBS_average = lapply(list.files(path="WGBS/shuffled/",pattern="_Meth_average.txt",full.names = TRUE),function(x) read.table(x,sep='\t',header=TRUE))
print("Loaded WGBS matrices")

# Number of CpGs per TE
shuffled_WGBS_CpG = lapply(list.files(path="WGBS/shuffled/",pattern="TE_CpG_count_",full.names = TRUE),function(x) read.table(x,sep='\t'))
shuffled_WGBS_CpG = lapply(shuffled_WGBS_CpG, setNames, nm =c("chromosome","start","stop","subfamily","class","family","strand","CpGs"))
print("Loaded CpG matrices")

# Methylation level for TEs with at least one CpG
for (i in 1:10){
  shuffled_WGBS_average[[i]] = merge(shuffled_WGBS_average[[i]],shuffled_WGBS_CpG[[i]],by=c("chromosome","start","stop","subfamily","class","family","strand"))
}
rm(shuffled_WGBS_CpG)
print("Combined WGBS and CpG matrices")

# Number of samples TE is in each methylation state (IMR90, no IMR90)
shuffled_WGBS_average = lapply(shuffled_WGBS_average,function(y) transform(y,Hypomethylated = apply(y,1,function(x) sum(x[8:44] < 0.3,na.rm=TRUE)),Hypermethylated = apply(y,1,function(x) sum(x[8:44] > 0.7,na.rm=TRUE)),Intermediate = apply(y,1,function(x) sum(x[8:44] <= 0.7 & x[8:44] >= 0.3,na.rm=TRUE)),Missing = apply(y,1,function(x) sum(is.na(x[8:44]))),Hypomethylated_noIMR90 = apply(y,1,function(x) sum(x[c(8:17,19:44)] < 0.3,na.rm=TRUE)),Hypermethylated_noIMR90 = apply(y,1,function(x) sum(x[c(8:17,19:44)] > 0.7,na.rm=TRUE)),Intermediate_noIMR90 = apply(y,1,function(x) sum(x[c(8:17,19:44)] <= 0.7 & x[c(8:17,19:44)] >= 0.3,na.rm=TRUE)),Missing_noIMR90 = apply(y,1,function(x) sum(is.na(x[c(8:17,19:44)])))))
print("Calculating methylation states")

save(shuffled_WGBS_average,file="R_scripts/shuffled_WGBS.RData")