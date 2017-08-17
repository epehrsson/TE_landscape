# Proportion of each TE class in state
# See 4/27/2016, 5/20/2016, 6/27/2016, 9/17/2016, 2/3/2017, 2/6/2017, 3/2/2017, 5/18/2017, 6/5/2017, 7/4/2017
# See 5/9/2016, 6/2/2016, 6/27/2016, 8/24/2016, 9/7/2016, 9/17/2016, 9/28/2016, 11/27/2016, 12/13/2016, 12/15/2016, 1/13/2017, 2/6/2017, 3/2/2017, 3/3/2017, 5/14/2017, 5/18/2017, 7/21/2017, 7/24/2017, 8/1/2017

# chromHMM
contribution_class = contribution[c(4:7,11:12,16),]
contribution_class = contribution_class$Total/rowSums(contribution_class[,2:16])
contribution_class_long = melt(as.matrix(contribution_class))
colnames(contribution_class_long) = c("State","Class","Proportion")

# WGBS
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

# DNase
# Contribution of DNase peaks overlapping TEs by class
TE_DNase_class = read.table("DNase_peaks/class_DNase_sample.txt",sep='\t')
colnames(TE_DNase_class) = c("class","Sample","Overlap")
TE_DNase_class = dcast(TE_DNase_class,Sample~class)
TE_DNase_class = merge(TE_DNase_class,DNase_stats[,c(1,5)],by=c("Sample"))
TE_DNase_class[,2:9] = t(apply(TE_DNase_class,1,function(x) as.numeric(x[2:9])/as.numeric(x[9])))

# Proportion of TEs overlapping DNase peak by class by sample
TE_DNase_peaks_class = aggregate(data=TE_DNase_peaks[,c(8:60,62)],.~class_update,function(x) sum(x > 0))
rownames(TE_DNase_peaks_class) = TE_DNase_peaks_class$class_update
TE_DNase_peaks_class = TE_DNase_peaks_class[,2:54]
TE_DNase_peaks_class = TE_DNase_peaks_class/rmsk_TEother_class_pi[match(rownames(TE_DNase_peaks_class),rmsk_TEother_class_pi$Class),]$Elements

# H3K27ac
# Contribution of H3K27ac peaks overlapping TEs by class
TE_H3K27ac_class = read.table("H3K27ac/class_H3K27ac_sample.txt",sep='\t')
colnames(TE_H3K27ac_class) = c("class","Sample","Overlap")
TE_H3K27ac_class = dcast(TE_H3K27ac_class,Sample~class)
TE_H3K27ac_class = merge(TE_H3K27ac_class,H3K27ac_stats[,c(1,5)],by=c("Sample"))
TE_H3K27ac_class[,2:9] = t(apply(TE_H3K27ac_class,1,function(x) as.numeric(x[2:9])/as.numeric(x[9])))

# Proportion of TEs overlapping DNase peak by class by sample
TE_H3K27ac_peaks_class = aggregate(data=TE_H3K27ac_peaks[,c(8:105,107)],.~class_update,function(x) sum(x > 0))
rownames(TE_H3K27ac_peaks_class) = TE_H3K27ac_peaks_class$class_update
TE_H3K27ac_peaks_class = TE_H3K27ac_peaks_class[,2:99]
TE_H3K27ac_peaks_class = TE_H3K27ac_peaks_class/rmsk_TEother_class_pi[match(rownames(TE_H3K27ac_peaks_class),rmsk_TEother_class_pi$Class),]$Elements