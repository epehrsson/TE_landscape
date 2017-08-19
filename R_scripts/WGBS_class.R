# CpG methylation state by class
# See 5/9/2016, 6/2/2016, 6/27/2016, 8/24/2016, 9/7/2016, 9/17/2016, 9/28/2016, 11/27/2016, 12/13/2016, 12/15/2016, 1/13/2017, 2/6/2017, 3/2/2017, 3/3/2017, 5/14/2017, 5/18/2017, 7/21/2017, 7/24/2017, 8/1/2017

# Proportion of class CpGs in methylation state by sample
class_CpG_meth = read.table("WGBS/class_CpG_Meth_states.txt",sep='\t')
class_CpG_meth = class_CpG_meth[which(class_CpG_meth$V1 != 41),]
class_CpG_meth$Sample = rep(as.vector(metadata[which(!is.na(metadata$WGBS)),]$Sample),8)
class_CpG_meth = class_CpG_meth[,c(7,6,2:5)]
colnames(class_CpG_meth)[3:6] = meth_states
colnames(class_CpG_meth)[2] = "class"
class_CpG_meth[is.na(class_CpG_meth)] = 0
class_CpG_meth = class_CpG_meth[which(class_CpG_meth$class %in% c("DNA","LINE","LTR","SINE","Other","Unconfident_RC")),]
class_CpG_meth$class = factor(c(rep("LINE",37),rep("SINE",37),rep("LTR",37),rep("DNA",37),rep("SVA",37),rep("Other",37)),levels=c("DNA","LINE","LTR","SINE","SVA","Other"))
class_CpG_meth[,3:6] = class_CpG_meth[,3:6]/rowSums(class_CpG_meth[,3:6])