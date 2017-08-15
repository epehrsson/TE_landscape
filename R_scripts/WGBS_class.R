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

# Contribution
# Proportion of CpGs in TEs in each state in each class
class_CpG_meth_contribution = read.table("WGBS/class_CpG_Meth_states.txt",sep='\t')
class_CpG_meth_contribution = class_CpG_meth_contribution[which(class_CpG_meth_contribution$V1 != 41),]
class_CpG_meth_contribution = class_CpG_meth_contribution[,c(6,2:5)]
colnames(class_CpG_meth_contribution) = c("class",meth_states)
class_CpG_meth_contribution[is.na(class_CpG_meth_contribution)] = 0
class_CpG_meth_contribution = aggregate(data=class_CpG_meth_contribution,.~class,sum)
class_CpG_meth_contribution$class = factor(class_CpG_meth_contribution$class,levels=c("Total",levels(class_CpG_meth_contribution$class)))
class_CpG_meth_contribution = rbind(class_CpG_meth_contribution,c("Total",colSums(TE_CpG_meth)))
class_CpG_meth_contribution[,2:5] = apply(class_CpG_meth_contribution[,2:5],2,function(x) as.numeric(x))
class_CpG_meth_contribution[,2:5] = apply(class_CpG_meth_contribution[,2:5],2,function(x) x/x[8])
class_CpG_meth_contribution = class_CpG_meth_contribution[1:7,]
class_CpG_meth_contribution$All = TE_class_CpG_count[as.vector(class_CpG_meth_contribution$class),]$Percent_TE_CpGs

# Proportion of CpGs in TEs in each state in each class, no IMR90
class_CpG_meth_contribution_noIMR90 = read.table("WGBS/class_CpG_Meth_states.txt",sep='\t')
class_CpG_meth_contribution_noIMR90 = class_CpG_meth_contribution_noIMR90[which(!(class_CpG_meth_contribution_noIMR90$V1 %in% c(11,41))),]
class_CpG_meth_contribution_noIMR90 = class_CpG_meth_contribution_noIMR90[,c(6,2:5)]
colnames(class_CpG_meth_contribution_noIMR90) = c("class",meth_states)
class_CpG_meth_contribution_noIMR90[is.na(class_CpG_meth_contribution_noIMR90)] = 0
class_CpG_meth_contribution_noIMR90 = aggregate(data=class_CpG_meth_contribution_noIMR90,.~class,sum)
class_CpG_meth_contribution_noIMR90$class = factor(class_CpG_meth_contribution_noIMR90$class,levels=c("Total",levels(class_CpG_meth_contribution_noIMR90$class)))
class_CpG_meth_contribution_noIMR90 = rbind(class_CpG_meth_contribution_noIMR90,c("Total",colSums(TE_CpG_meth[which(rownames(TE_CpG_meth) != "E017"),])))
class_CpG_meth_contribution_noIMR90[,2:5] = apply(class_CpG_meth_contribution_noIMR90[,2:5],2,function(x) as.numeric(x)/as.numeric(x)[8])
class_CpG_meth_contribution_noIMR90 = class_CpG_meth_contribution_noIMR90[1:7,]
class_CpG_meth_contribution_noIMR90$All = TE_class_CpG_count[as.vector(class_CpG_meth_contribution_noIMR90$class),]$Percent_TE_CpGs

# Potential
# Cumulative distribution of methylation states and statistics, by class
TE_meth_average_class = by(TE_meth_average,TE_meth_average$class_update,function(x) sample_distribution(x,c(46:49),37))
TE_meth_average_class_cum = by(TE_meth_average,TE_meth_average$class_update,function(x) cumulative_distribution(x,c(46:49),37))
TE_meth_average_class_cum_long = rbind(melt(as.matrix(TE_meth_average_class_cum$DNA)),melt(as.matrix(TE_meth_average_class_cum$LINE)),melt(as.matrix(TE_meth_average_class_cum$LTR)),melt(as.matrix(TE_meth_average_class_cum$SINE)),melt(as.matrix(TE_meth_average_class_cum$Other)),melt(as.matrix(TE_meth_average_class_cum$RC)),melt(as.matrix(TE_meth_average_class_cum$Unconfident)))
TE_meth_average_class_cum_long$Class = c(rep("DNA",148),rep("LINE",148),rep("LTR",148),rep("SINE",148),rep("Other",148),rep("RC",148),rep("Unconfident",148))
colnames(TE_meth_average_class_cum_long)[1:3] = c("Samples","Methylation","Proportion")
TE_meth_average_class_cum_long$Group = paste(TE_meth_average_class_cum_long$Methylation,TE_meth_average_class_cum_long$Class,sep="_")
TE_meth_average_class_stats = ldply(TE_meth_average_class,function(y) as.data.frame(t(rbind(apply(y[,2:5],2,function(x) sum(x[2:38])/(sum(x)/100)),apply(y[,2:5],2,function(x) sum(as.numeric(x)*seq(0,37))/sum(x))/37,apply(y[2:38,2:5],2,function(x) sum(as.numeric(x)*seq(1,37))/sum(x))/37))))
colnames(TE_meth_average_class_stats) = c("Class","Proportion_ever","Samples_avg_all","Samples_avg_ever")
TE_meth_average_class_stats$Class = factor(TE_meth_average_class_stats$Class,levels=c("DNA","LINE","LTR","SINE","Other","RC","Unconfident"))
TE_meth_average_class_stats$State = factor(rep(meth_states[c(1,3,2,4)],7),levels=meth_states)
TE_meth_average_class_stats[,2:4] = apply(TE_meth_average_class_stats[,2:4],2,function(x) as.numeric(x))

# Cumulative distribution of methylation states and statistics, no IMR90, by class
TE_meth_average_class_noIMR90 = by(TE_meth_average,TE_meth_average$class_update,function(x) sample_distribution(x,c(50:53),36))
TE_meth_average_class_noIMR90_cum = by(TE_meth_average,TE_meth_average$class_update,function(x) cumulative_distribution(x,c(50:53),36))
TE_meth_average_class_noIMR90_cum_long = rbind(melt(as.matrix(TE_meth_average_class_noIMR90_cum$DNA)),melt(as.matrix(TE_meth_average_class_noIMR90_cum$LINE)),melt(as.matrix(TE_meth_average_class_noIMR90_cum$LTR)),melt(as.matrix(TE_meth_average_class_noIMR90_cum$SINE)),melt(as.matrix(TE_meth_average_class_noIMR90_cum$Other)),melt(as.matrix(TE_meth_average_class_noIMR90_cum$RC)),melt(as.matrix(TE_meth_average_class_noIMR90_cum$Unconfident)))
TE_meth_average_class_noIMR90_cum_long$Class = c(rep("DNA",144),rep("LINE",144),rep("LTR",144),rep("SINE",144),rep("Other",144),rep("RC",144),rep("Unconfident",144))
colnames(TE_meth_average_class_noIMR90_cum_long)[1:3] = c("Samples","Methylation","Proportion")
TE_meth_average_class_noIMR90_cum_long$Group = paste(TE_meth_average_class_noIMR90_cum_long$Methylation,TE_meth_average_class_noIMR90_cum_long$Class,sep="_")
TE_meth_average_class_noIMR90_stats = ldply(TE_meth_average_class_noIMR90,function(y) as.data.frame(t(rbind(apply(y[,2:5],2,function(x) sum(x[2:37])/(sum(x)/100)),apply(y[,2:5],2,function(x) sum(as.numeric(x)*seq(0,36))/sum(x))/36,apply(y[2:37,2:5],2,function(x) sum(as.numeric(x)*seq(1,36))/sum(x))/36))))
colnames(TE_meth_average_class_noIMR90_stats) = c("Class","Proportion_ever","Samples_avg_all","Samples_avg_ever")
TE_meth_average_class_noIMR90_stats$Class = factor(TE_meth_average_class_noIMR90_stats$Class,levels=c("DNA","LINE","LTR","SINE","Other","RC","Unconfident"))
TE_meth_average_class_noIMR90_stats$State = factor(rep(meth_states[c(1,3,2,4)],7),levels=meth_states)
TE_meth_average_class_noIMR90_stats[,2:4] = apply(TE_meth_average_class_noIMR90_stats[,2:4],2,function(x) as.numeric(x))