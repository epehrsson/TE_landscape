# Proportion of each TE class in state (overall, by sample)
# See 4/27/2016, 5/20/2016, 6/27/2016, 9/17/2016, 2/3/2017, 2/6/2017, 3/2/2017, 5/18/2017, 6/5/2017, 7/4/2017
# See 5/9/2016, 6/2/2016, 6/27/2016, 8/24/2016, 9/7/2016, 9/17/2016, 9/28/2016, 11/27/2016, 12/13/2016, 12/15/2016, 1/13/2017, 2/6/2017, 3/2/2017, 3/3/2017, 5/14/2017, 5/18/2017, 7/21/2017, 7/24/2017, 8/1/2017

#source("R_scripts/TE_class_stats.R")
#source("R_scripts/proportion.R")

# chromHMM
# By sample
class_chromHMM = read.table("chromHMM/class/class_state_sample.txt",sep='\t')
colnames(class_chromHMM) = c("class","State","Sample","Bases_state_class")
class_chromHMM$class = convert_class(class_chromHMM$class)
class_chromHMM = merge(class_chromHMM,expand.grid(class=levels(class_chromHMM$class),Sample=levels(class_chromHMM$Sample),State=levels(class_chromHMM$State)),by=c("class","Sample","State"),all.y=TRUE)
class_chromHMM[is.na(class_chromHMM)] = 0
class_chromHMM$State = factor(class_chromHMM$State,chromHMM_states)

# Add total bases in class in sample
class_chromHMM = merge(class_chromHMM,aggregate(data=class_chromHMM,Bases_state_class~Sample+class,sum),by=c("Sample","class"))
colnames(class_chromHMM)[4:5] = c("Bases_state_class","Bases_class")
class_chromHMM$Proportion_class = class_chromHMM$Bases_state_class/class_chromHMM$Bases_class

# Add total bases in state in sample
class_chromHMM$Bases_state = apply(class_chromHMM,1,function(x) mnemonics_states_genome[x[3],x[1]])
class_chromHMM$Proportion_state = class_chromHMM$Bases_state_class/class_chromHMM$Bases_state

# Add total bases in sample
class_chromHMM$Bases_sample = apply(class_chromHMM,1,function(x) mnemonics_states_genome[1,x[1]])
class_chromHMM$Enrichment = log2((class_chromHMM$Bases_state_class/class_chromHMM$Bases_class)/(class_chromHMM$Bases_state/class_chromHMM$Bases_sample))

# Total
contribution_class = aggregate(data=class_chromHMM,Bases_state_class~State+class,sum)
contribution_class$Proportion = apply(contribution_class,1,function(x) as.numeric(x[3])/rmsk_TE_class[match(x[2],rmsk_TE_class$class_update),]$chromHMM_total_width)

# WGBS
# By sample
class_CpG_meth = read.table("WGBS/class_CpG_Meth_states.txt",sep='\t')
class_CpG_meth = class_CpG_meth[which(class_CpG_meth$V1 != 41),]
class_CpG_meth$Sample = rep(as.vector(metadata[which(!is.na(metadata$WGBS)),]$Sample),8)
class_CpG_meth = class_CpG_meth[,c(7,6,2:5)]
colnames(class_CpG_meth)[2:6] = c("class",meth_states)
class_CpG_meth[is.na(class_CpG_meth)] = 0
class_CpG_meth = class_CpG_meth[which(class_CpG_meth$class %in% c("DNA","LINE","LTR","SINE","Other","Unconfident_RC")),]
class_CpG_meth$class = convert_class(class_CpG_meth$class)

# Total
class_CpG_meth_total = aggregate(data=class_CpG_meth[,2:6],.~class,sum)
class_CpG_meth_total[,2:5] = class_CpG_meth_total[,2:5]/(rmsk_TE_class[match(class_CpG_meth_total$class,rmsk_TE_class$class_update),]$CpGs*37)
class_CpG_meth_total = melt(class_CpG_meth_total,id.vars="class")
colnames(class_CpG_meth_total) = c("class","State","Proportion")

# Total, no IMR90
class_CpG_meth_total_noIMR90 = aggregate(data=class_CpG_meth[which(class_CpG_meth$Sample != "E017"),2:6],.~class,sum)
class_CpG_meth_total_noIMR90[,2:5] = class_CpG_meth_total_noIMR90[,2:5]/(rmsk_TE_class[match(class_CpG_meth_total_noIMR90$class,rmsk_TE_class$class_update),]$CpGs*36)

# By sample - Add CpGs in class in sample
class_CpG_meth$CpGs_class = rowSums(class_CpG_meth[,3:6])
class_CpG_meth = melt(class_CpG_meth,id.vars=c("class","Sample","CpGs_class"))
colnames(class_CpG_meth)[4:5] = c("State","CpGs_state_class")
class_CpG_meth$State = factor(class_CpG_meth$State,meth_states)
class_CpG_meth$Proportion_class = class_CpG_meth$CpGs_state_class/class_CpG_meth$CpGs_class

# Add CpGs in state in sample
class_CpG_meth$CpGs_state = apply(class_CpG_meth,1,function(x) all_CpG_meth[x[2],x[4]])
class_CpG_meth$Proportion_state = class_CpG_meth$CpGs_state_class/class_CpG_meth$CpGs_state

# Add CpGs in sample
class_CpG_meth$CpGs_sample = rep(ALL_CPGS,dim(class_CpG_meth)[1])
class_CpG_meth$Enrichment = log2((class_CpG_meth$CpGs_state_class/class_CpG_meth$CpGs_class)/(class_CpG_meth$CpGs_state/class_CpG_meth$CpGs_sample))

# DNase
# By sample
TE_DNase_class = read.table("DNase/class_DNase_sample.txt",sep='\t')
colnames(TE_DNase_class) = c("class","Sample","Bases_state_class")
TE_DNase_class = TE_DNase_class[which(TE_DNase_class$class %in% c("DNA","LINE","LTR","SINE","Other","Unconfident_RC")),]
TE_DNase_class$class = convert_class(TE_DNase_class$class)
TE_DNase_class$State = rep("DNase",dim(TE_DNase_class)[1])

# Add bases in class in sample
TE_DNase_class$Bases_class = apply(TE_DNase_class,1,function(x) rmsk_TE_class[match(x[1],rmsk_TE_class$class_update),]$Total_length)
TE_DNase_class[which(metadata[match(TE_DNase_class$Sample,metadata$Sample),]$chrY == "No"),]$Bases_class = apply(TE_DNase_class[which(metadata[match(TE_DNase_class$Sample,metadata$Sample),]$chrY == "No"),],1,function(x) rmsk_TE_class[match(x[1],rmsk_TE_class$class_update),]$Total_length_noY)
TE_DNase_class$Proportion_class = TE_DNase_class$Bases_state_class/TE_DNase_class$Bases_class

# Add bases in state in sample
TE_DNase_class$Bases_state = apply(TE_DNase_class,1,function(x) DNase_stats[match(x[2],DNase_stats$Sample),]$Total_width)
TE_DNase_class$Proportion_state = TE_DNase_class$Bases_state_class/TE_DNase_class$Bases_state

# Add total bases in sample
TE_DNase_class$Bases_sample = apply(TE_DNase_class,1,function(x) mnemonics_states_genome[1,x[2]])
TE_DNase_class$Enrichment = log2((TE_DNase_class$Bases_state_class/TE_DNase_class$Bases_class)/(TE_DNase_class$Bases_state/TE_DNase_class$Bases_sample))

# Total
TE_DNase_class_total = colSums(dcast(TE_DNase_class[,1:3],Sample~class)[,2:7])
TE_DNase_class_total = melt(as.matrix(TE_DNase_class_total/rmsk_TE_class[match(names(TE_DNase_class_total),rmsk_TE_class$class_update),]$DNase_total_width))
colnames(TE_DNase_class_total) = c("class","State","Proportion")
TE_DNase_class_total$State = rep("DNase",6)

# H3K27ac
# By sample
TE_H3K27ac_class = read.table("H3K27ac/class_H3K27ac_sample.txt",sep='\t')
colnames(TE_H3K27ac_class) = c("class","Sample","Bases_state_class")
TE_H3K27ac_class = TE_H3K27ac_class[which(TE_H3K27ac_class$class %in% c("DNA","LINE","LTR","SINE","Other","Unconfident_RC")),]
TE_H3K27ac_class$class = convert_class(TE_H3K27ac_class$class)
TE_H3K27ac_class$State = rep("H3K27ac",dim(TE_H3K27ac_class)[1])

# Add bases in class in sample
TE_H3K27ac_class$Bases_class = apply(TE_H3K27ac_class,1,function(x) rmsk_TE_class[match(x[1],rmsk_TE_class$class_update),]$Total_length)
TE_H3K27ac_class[which(metadata[match(TE_H3K27ac_class$Sample,metadata$Sample),]$chrY == "No"),]$Bases_class = apply(TE_H3K27ac_class[which(metadata[match(TE_H3K27ac_class$Sample,metadata$Sample),]$chrY == "No"),],1,function(x) rmsk_TE_class[match(x[1],rmsk_TE_class$class_update),]$Total_length_noY)
TE_H3K27ac_class$Proportion_class = TE_H3K27ac_class$Bases_state_class/TE_H3K27ac_class$Bases_class

# Add bases in state in sample
TE_H3K27ac_class$Bases_state = apply(TE_H3K27ac_class,1,function(x) H3K27ac_stats[match(x[2],H3K27ac_stats$Sample),]$Total_width)
TE_H3K27ac_class$Proportion_state = TE_H3K27ac_class$Bases_state_class/TE_H3K27ac_class$Bases_state

# Add total bases in sample
TE_H3K27ac_class$Bases_sample = apply(TE_H3K27ac_class,1,function(x) mnemonics_states_genome[1,x[2]])
TE_H3K27ac_class$Enrichment = log2((TE_H3K27ac_class$Bases_state_class/TE_H3K27ac_class$Bases_class)/(TE_H3K27ac_class$Bases_state/TE_H3K27ac_class$Bases_sample))

# Total
TE_H3K27ac_class_total = colSums(dcast(TE_H3K27ac_class[,1:3],Sample~class)[,2:7])
TE_H3K27ac_class_total = melt(as.matrix(TE_H3K27ac_class_total/rmsk_TE_class[match(names(TE_H3K27ac_class_total),rmsk_TE_class$class_update),]$H3K27ac_total_width))
colnames(TE_H3K27ac_class_total) = c("class","State","Proportion")
TE_H3K27ac_class_total$State = rep("H3K27ac",6)