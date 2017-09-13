# Load matrix of methylation level per TE
# See 6/2/2016, 9/17/2016, 12/15/2016, 2/6/2017, 7/21/2017, 7/22/2017

# Methylation level for all TEs
TE_meth_average = read.table("WGBS/TE_CpG_Meth_new_average.txt",sep='\t',header=TRUE)

# Number of CpGs per TE
TE_CpG_count = read.table("WGBS/TE_CpG_count.txt",sep='\t')
colnames(TE_CpG_count) = c("chromosome","start","stop","subfamily","class","family","strand","CpGs")

# Methylation level for TEs with at least one CpG
TE_meth_average = merge(TE_meth_average,TE_CpG_count,by=c("chromosome","start","stop","subfamily","class","family","strand"))
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
TE_meth_average$class_update = TE_meth_average$class
TE_meth_average$class_update = factor(TE_meth_average$class_update,levels=c("DNA","LINE","LTR","SINE","SVA","Other"))
TE_meth_average[which(TE_meth_average$class == "Other"),]$class_update = "SVA"
TE_meth_average[which(TE_meth_average$class %in% c("DNA?","LINE?","LTR?","SINE?","Unknown","Unknown?","RC")),]$class_update = "Other"

save(TE_meth_average,file="R_datasets/TE_meth_average.RData")
