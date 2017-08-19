# Proportion of each TE class in state
# See 4/27/2016, 5/20/2016, 6/27/2016, 9/17/2016, 2/3/2017, 2/6/2017, 3/2/2017, 5/18/2017, 6/5/2017, 7/4/2017
# See 5/9/2016, 6/2/2016, 6/27/2016, 8/24/2016, 9/7/2016, 9/17/2016, 9/28/2016, 11/27/2016, 12/13/2016, 12/15/2016, 1/13/2017, 2/6/2017, 3/2/2017, 3/3/2017, 5/14/2017, 5/18/2017, 7/21/2017, 7/24/2017, 8/1/2017

source("R_scripts/contribution.R")
source("R_scripts/TE_class_stats.R")

# chromHMM
contribution_class = contribution[c(4:7,11,17),1:16]
rownames(contribution_class)[5:6] = c("SVA","Other")
contribution_class = contribution_class[,2:16]/contribution_class$Total
contribution_class_long = melt(as.matrix(contribution_class))
colnames(contribution_class_long) = c("Class","State","Proportion")

# WGBS
# Proportion of CpGs in TEs in each state in each class
class_CpG_meth_contribution = read.table("WGBS/class_CpG_Meth_states.txt",sep='\t')
class_CpG_meth_contribution = class_CpG_meth_contribution[which(class_CpG_meth_contribution$V1 != 41),]
class_CpG_meth_contribution = class_CpG_meth_contribution[,c(6,2:5)]
colnames(class_CpG_meth_contribution) = c("class",meth_states)
class_CpG_meth_contribution[is.na(class_CpG_meth_contribution)] = 0
class_CpG_meth_contribution = aggregate(data=class_CpG_meth_contribution,.~class,sum)[c(1:4,6,8),]
class_CpG_meth_contribution$class = factor(c("DNA","LINE","LTR","SVA","SINE","Other"),levels=c("DNA","LINE","LTR","SINE","Other","SVA"))
class_CpG_meth_contribution[,2:5] = class_CpG_meth_contribution[,2:5]/(rmsk_TE_class[match(class_CpG_meth_contribution$class,rmsk_TE_class$class_update),]$CpGs*37)
class_CpG_meth_contribution_long = melt(class_CpG_meth_contribution,id.vars="class")
colnames(class_CpG_meth_contribution_long) = c("Class","State","Proportion")

# Proportion of CpGs in TEs in each state in each class, no IMR90
class_CpG_meth_contribution_noIMR90 = read.table("WGBS/class_CpG_Meth_states.txt",sep='\t')
class_CpG_meth_contribution_noIMR90 = class_CpG_meth_contribution_noIMR90[which(!(class_CpG_meth_contribution_noIMR90$V1 %in% c(11,41))),]
class_CpG_meth_contribution_noIMR90 = class_CpG_meth_contribution_noIMR90[,c(6,2:5)]
colnames(class_CpG_meth_contribution_noIMR90) = c("class",meth_states)
class_CpG_meth_contribution_noIMR90[is.na(class_CpG_meth_contribution_noIMR90)] = 0
class_CpG_meth_contribution_noIMR90 = aggregate(data=class_CpG_meth_contribution_noIMR90,.~class,sum)[c(1:4,6,8),]
class_CpG_meth_contribution_noIMR90$class = factor(c("DNA","LINE","LTR","SVA","SINE","Other"),levels=c("DNA","LINE","LTR","SINE","Other","SVA"))
class_CpG_meth_contribution_noIMR90[,2:5] = class_CpG_meth_contribution_noIMR90[,2:5]/(rmsk_TE_class[match(class_CpG_meth_contribution_noIMR90$class,rmsk_TE_class$class_update),]$CpGs*36)

# DNase
# Contribution of DNase peaks overlapping TEs by class
TE_DNase_class = read.table("DNase/class_DNase_sample.txt",sep='\t')
colnames(TE_DNase_class) = c("class","Sample","Overlap")
TE_DNase_class = colSums(dcast(TE_DNase_class,Sample~class)[,c(2:5,7,9)])
names(TE_DNase_class)[c(4,6)] = c("SVA","Other")
TE_DNase_class = melt(as.matrix(TE_DNase_class/rmsk_TE_class[match(names(TE_DNase_class),rmsk_TE_class$class_update),]$DNase_total_width))
colnames(TE_DNase_class) = c("Class","State","Proportion")
TE_DNase_class$State = rep("DNase",6)

# H3K27ac
# Contribution of H3K27ac peaks overlapping TEs by class
TE_H3K27ac_class = read.table("H3K27ac/class_H3K27ac_sample.txt",sep='\t')
colnames(TE_H3K27ac_class) = c("class","Sample","Overlap")
TE_H3K27ac_class = colSums(dcast(TE_H3K27ac_class,Sample~class)[,c(2:5,7,9)])
names(TE_H3K27ac_class)[c(4,6)] = c("SVA","Other")
TE_H3K27ac_class = melt(as.matrix(TE_H3K27ac_class/rmsk_TE_class[match(names(TE_H3K27ac_class),rmsk_TE_class$class_update),]$H3K27ac_total_width))
colnames(TE_H3K27ac_class) = c("Class","State","Proportion")
TE_H3K27ac_class$State = rep("H3K27ac",6)