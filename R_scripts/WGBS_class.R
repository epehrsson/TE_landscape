# CpG methylation state by class
# See 5/9/2016, 6/2/2016, 6/27/2016, 8/24/2016, 9/7/2016, 9/17/2016, 9/28/2016, 11/27/2016, 12/13/2016, 12/15/2016, 1/13/2017, 2/6/2017, 3/2/2017, 3/3/2017, 5/14/2017, 5/18/2017, 7/21/2017, 7/24/2017, 8/1/2017

# Number of CpGs per class
TE_class_CpG_count = read.table("WGBS/class/TE_CpG_class.txt",sep='\t',row.names=1)
colnames(TE_class_CpG_count) = c("CpGs")
TE_class_CpG_count$TEs = rmsk_TEother_stats_class[match(rownames(TE_class_CpG_count),rmsk_TEother_stats_class$Class),]$Count
TE_class_CpG_count$TEs_wCpG = as.data.frame(table(TE_CpG_count$class))[match(rownames(TE_class_CpG_count),as.data.frame(table(TE_CpG_count$class))$Var1),]$Freq
TE_class_CpG_count[13,2] = sum(TE_class_CpG_count$TEs[c(2,4,6:8,10)])
TE_class_CpG_count[13,3] = sum(TE_class_CpG_count$TEs_wCpG[c(2,4,6:8,10)])

# CpGs per TE, number/proportion of TEs with CpGs
TE_class_CpG_count$Mean = TE_class_CpG_count$CpGs/TE_class_CpG_count$TEs
TE_class_CpG_count$Mean_wCpG = TE_class_CpG_count$CpGs/TE_class_CpG_count$TEs_wCpG
TE_class_CpG_count$TEs_wCpG_per = TE_class_CpG_count$TEs_wCpG/TE_class_CpG_count$TEs

# Proportion of CpGs (all, overlapping TEs) in each class
TE_class_CpG_count$Percent_TE_CpGs = TE_class_CpG_count$CpGs/28373958
TE_class_CpG_count$Percent_all_CpGs = TE_class_CpG_count$CpGs/56434896

# Proportion of class CpGs in methylation state by sample
class_CpG_meth = read.table("WGBS/class_CpG_Meth_states.txt",sep='\t')
class_CpG_meth = class_CpG_meth[which(class_CpG_meth$V1 != 41),]
class_CpG_meth$Sample = rep(EID_metadata_meth$Sample,7)
class_CpG_meth = class_CpG_meth[,c(7,6,2:5)]
colnames(class_CpG_meth)[3:6] = meth_states
colnames(class_CpG_meth)[2] = "class"
class_CpG_meth[is.na(class_CpG_meth)] = 0
class_CpG_meth[,3:6] = class_CpG_meth[,3:6]/rowSums(class_CpG_meth[,3:6])

# Comparison of CpG % hypomethylated by class
kruskal.test(class_CpG_meth$Hypomethylated~class_CpG_meth$class)
pairwise.wilcox.test(class_CpG_meth$Hypomethylated,class_CpG_meth$class,p.adj="bonf")
kruskal.test(class_CpG_meth[which(class_CpG_meth$Sample != "E017"),]$Hypomethylated~class_CpG_meth[which(class_CpG_meth$Sample != "E017"),]$class)
pairwise.wilcox.test(class_CpG_meth[which(class_CpG_meth$Sample != "E017"),]$Hypomethylated,class_CpG_meth[which(class_CpG_meth$Sample != "E017"),]$class,p.adj="bonf")