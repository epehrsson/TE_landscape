# Creates dataframes of the number of TEs in each epigenetic state in each number of samples, by class
# As well as the number of TEs in each state in at least one sample, by class

## chromHMM_TE_state_class: Number of TEs in each chromHMM state for each number of samples, by class
## TE_meth_average_class: Number of TEs in each methylation state for each number of samples, by class
## combine_stats_class: Proportion of TEs ever in each state and mean/SE proportion of samples in state, by class
## combine_boxplot_class: TEs in each state per sample, by class
## combine_potential_class: Number of TEs in each state in each number/proportion of samples, by class,
## plus the proportion of TEs at each number/proportion of samples that belongs to each TE class
## combine_class_ever: Proportion of TEs ever in each state that belongs to each TE class

# chromHMM

# Number of TEs in each chromHMM state for each number of samples (0-127), by class
chromHMM_TE_state_class = ddply(chromHMM_TE_state,~class_update,function(x) sample_distribution(x,c(8:22),sample_counts["All","chromHMM"]))

# Proportion of TEs ever in each state and mean/SE proportion of samples in state, for all TEs and those ever in the state, by class
chromHMM_TE_state_class_stats = ddply(chromHMM_TE_state_class,~class_update,function(x) potential_stats(x[,2:17],15,sample_counts["All","chromHMM"]))
chromHMM_TE_state_class_stats$State = factor(rep(chromHMM_states,6),levels=chromHMM_states)
chromHMM_TE_state_class_stats$class_update = factor(chromHMM_TE_state_class_stats$class_update,levels=c("DNA","LINE","LTR","SINE","SVA","Other"))
chromHMM_TE_state_class_stats[,2:4] = apply(chromHMM_TE_state_class_stats[,2:4],2,function(x) as.numeric(x))

# Number/proportion of TEs in each chromHMM state by class by sample
class_state_sample = read.table("chromHMM/class_state_sample_summit.txt",sep='\t',col.names=c("Class","Sample","State","Count"))

## Including missing combinations
class_state_sample = merge(class_state_sample,expand.grid(Class = levels(class_state_sample$Class),Sample = levels(class_state_sample$Sample),State = levels(class_state_sample$State)),by=c("Class","State","Sample"),all.y=TRUE)
class_state_sample[is.na(class_state_sample)] = 0
class_state_sample$Count = as.numeric(class_state_sample$Count)
class_state_sample$Proportion = ifelse(metadata[match(class_state_sample$Sample,metadata$Sample),]$chrY == "Yes",
                                       class_state_sample$Count/as.numeric(rmsk_TE_class[match(class_state_sample$Class,rmsk_TE_class$class_update),]$Count),
                                       class_state_sample$Count/as.numeric(rmsk_TE_class[match(class_state_sample$Class,rmsk_TE_class$class_update),]$Count_noY))
class_state_sample$State = factor(class_state_sample$State,levels=chromHMM_states)                                                                                            

# WGBS

# Number of TEs in each methylation state for each number of samples (0-37), by class
TE_meth_average_class = ddply(TE_meth_average,~class_update,function(x) sample_distribution(x,c(46:49),sample_counts["All","WGBS"]))

# Proportion of TEs ever in each state and mean/SE proportion of samples in state, for all TEs and those ever in the state, by class
TE_meth_average_class_stats = ddply(TE_meth_average_class,~class_update,function(x) potential_stats(x[,2:6],4,sample_counts["All","WGBS"]))
TE_meth_average_class_stats$State = factor(rep(meth_states[c(1,3,2,4)],6),levels=meth_states)
TE_meth_average_class_stats$class_update = factor(TE_meth_average_class_stats$class_update,levels=c("DNA","LINE","LTR","SINE","SVA","Other"))
TE_meth_average_class_stats[,2:4] = apply(TE_meth_average_class_stats[,2:4],2,function(x) as.numeric(x))

# Number/proportion of TEs in each methylation state by class by sample
TE_meth_average_state_class = melt(TE_meth_average[,c(1:44,54)],id.vars=c(TE_coordinates,"class_update"),
                                   variable.name="Sample",value.name="Methylation")
WGBS_sample_state_class = ddply(TE_meth_average_state_class,.(class_update,Sample),summarise,
                                Hypomethylated=sum(na.omit(Methylation) < 0.3),Intermediate=sum(na.omit(Methylation) <= 0.7 & na.omit(Methylation) >= 0.3),
                                Hypermethylated=sum(na.omit(Methylation) > 0.7),Missing=sum(is.na(Methylation)))
WGBS_sample_state_class = melt(WGBS_sample_state_class,id.vars=c("class_update","Sample"))
colnames(WGBS_sample_state_class) = c("Class","Sample","State","Count")
WGBS_sample_state_class$Proportion = WGBS_sample_state_class$Count/as.numeric(rmsk_TE_class[match(WGBS_sample_state_class$Class,rmsk_TE_class$class_update),]$Count)

# DHS

# Number of TEs overlapping a DHS peak summit for each number of samples (0-53), by class
colnames(TE_DNase_peaks)[62] = "DNase"
potential_TE_DNase_class = ddply(TE_DNase_peaks,~class_update,function(x) sample_distribution(x,62,sample_counts["All","DNase"]))
potential_TE_DNase_class[which(potential_TE_DNase_class$Samples == 0),3] =  rmsk_TE_class$Count - aggregate(data=potential_TE_DNase_class,DNase~class_update,sum)$DNase

# Proportion of TEs ever in each state and mean/SE proportion of samples in state, for all TEs and those ever in the state, by class
potential_TE_DNase_class_stats = ddply(potential_TE_DNase_class,~class_update,function(x) potential_stats(x[,2:3],1,sample_counts["All","DNase"]))
potential_TE_DNase_class_stats$State = rep("DNase",6)

# Number/proportion of TEs overlapping a DHS peak summit by class by sample
TE_DNase_peaks_class = aggregate(data=TE_DNase_peaks[,c(8:61)],.~class_update,function(x) sum(x > 0))
TE_DNase_peaks_class = melt(TE_DNase_peaks_class,id.var="class_update")
colnames(TE_DNase_peaks_class) = c("Class","Sample","Count")
TE_DNase_peaks_class$Proportion = TE_DNase_peaks_class$Count/ifelse(metadata[match(TE_DNase_peaks_class$Sample,metadata$Sample),]$chrY == "Yes",
                                                                    rmsk_TE_class[match(TE_DNase_peaks_class$Class,rmsk_TE_class$class_update),]$Count,
                                                                    rmsk_TE_class[match(TE_DNase_peaks_class$Class,rmsk_TE_class$class_update),]$Count_noY)
TE_DNase_peaks_class$State = rep("DNase",dim(TE_DNase_peaks_class)[1])

# Revert column name
colnames(TE_DNase_peaks)[62] = "Samples"

# H3K27ac

# Number of TEs overlapping an H3K27ac peak summit for each number of samples (0-98), by class
colnames(TE_H3K27ac_peaks)[107] = "H3K27ac"
potential_TE_H3K27ac_class = ddply(TE_H3K27ac_peaks,~class_update,function(x) sample_distribution(x,107,sample_counts["All","H3K27ac"]))
potential_TE_H3K27ac_class[which(potential_TE_H3K27ac_class$Samples == 0),3] =  rmsk_TE_class$Count - aggregate(data=potential_TE_H3K27ac_class,H3K27ac~class_update,sum)$H3K27ac

# Proportion of TEs ever in each state and mean/SE proportion of samples in state, for all TEs and those ever in the state, by class
potential_TE_H3K27ac_class_stats = ddply(potential_TE_H3K27ac_class,~class_update,function(x) potential_stats(x[,2:3],1,sample_counts["All","H3K27ac"]))
potential_TE_H3K27ac_class_stats$State = rep("H3K27ac",6)

# Number/proportion of TEs overlapping an H3K27ac peak summit by class by sample
TE_H3K27ac_peaks_class = aggregate(data=TE_H3K27ac_peaks[,c(8:106)],.~class_update,function(x) sum(x > 0))
TE_H3K27ac_peaks_class = melt(TE_H3K27ac_peaks_class,id.var="class_update")
colnames(TE_H3K27ac_peaks_class) = c("Class","Sample","Count")
TE_H3K27ac_peaks_class$Proportion = TE_H3K27ac_peaks_class$Count/ifelse(metadata[match(TE_H3K27ac_peaks_class$Sample,metadata$Sample),]$chrY == "Yes",
                                                                    rmsk_TE_class[match(TE_H3K27ac_peaks_class$Class,rmsk_TE_class$class_update),]$Count,
                                                                    rmsk_TE_class[match(TE_H3K27ac_peaks_class$Class,rmsk_TE_class$class_update),]$Count_noY)
TE_H3K27ac_peaks_class$State = rep("H3K27ac",dim(TE_H3K27ac_peaks_class)[1])

# Revert column name
colnames(TE_H3K27ac_peaks)[107] = "Samples"

# Expression

# Number of TEs RPKM >1 for each number of samples (0-56), by class
RNA_potential_class = ddply(RNA_TE,~class_update,function(x) sample_distribution(x,64,sample_counts["All","RNA"]))

# Proportion of TEs ever in each state and mean/SE proportion of samples in state, for all TEs and those ever in the state, by class
RNA_potential_class_stats = ddply(RNA_potential_class,~class_update,function(x) potential_stats(x[,2:3],1,sample_counts["All","RNA"]))
RNA_potential_class_stats$class_update = factor(RNA_potential_class_stats$class_update,levels=c("DNA","LINE","LTR","SINE","SVA","Other"))
RNA_potential_class_stats$State = rep("Expressed_samples",6)

# Number/proportion of TEs expressed RPKM > 1 by class by sample
RNA_RPKM_class = aggregate(data=RNA_TE[,c(8:63,67)],.~class_update,function(x) sum(x > 1))
RNA_RPKM_class = melt(RNA_RPKM_class,id.var="class_update")
colnames(RNA_RPKM_class) = c("Class","Sample","Count")
RNA_RPKM_class$Proportion = RNA_RPKM_class$Count/ifelse(metadata[match(RNA_RPKM_class$Sample,metadata$Sample),]$chrY == "Yes",
                                                                        rmsk_TE_class[match(RNA_RPKM_class$Class,rmsk_TE_class$class_update),]$Count,
                                                                        rmsk_TE_class[match(RNA_RPKM_class$Class,rmsk_TE_class$class_update),]$Count_noY)
RNA_RPKM_class$State = rep("Expressed_samples",dim(RNA_RPKM_class)[1])


# Combine dataframes for all techniques
## Proportion of TEs ever in each state and mean/SE proportion of samples in state, by class
combine_stats_class = rbind(chromHMM_TE_state_class_stats,TE_meth_average_class_stats,potential_TE_DNase_class_stats,potential_TE_H3K27ac_class_stats,RNA_potential_class_stats)
combine_stats_class$Group = factor(c(rep("chromHMM",90),rep("WGBS",24),rep("DNase",6),rep("H3K27ac",6),rep("Expressed_samples",6)),levels=c("chromHMM","WGBS","DNase","H3K27ac","Expressed_samples"))
colnames(combine_stats_class)[1] = "Class"

## TEs in each state per sample, by class
combine_boxplot_class = rbind(class_state_sample,WGBS_sample_state_class,TE_DNase_peaks_class,TE_H3K27ac_peaks_class,RNA_RPKM_class)

## Number of TEs in each state in each number/proportion of samples, by class
## As well as the number of all TEs in each state in each number/proportion of samples
## And the proportion of TEs at each number/proportion of samples that belongs to each TE class
combine_potential_class = rbind(melt(chromHMM_TE_state_class,id.vars=c("Samples","class_update")),melt(TE_meth_average_class,id.vars=c("Samples","class_update")),
                          melt(potential_TE_DNase_class,id.var=c("Samples","class_update")),melt(potential_TE_H3K27ac_class,id.var=c("Samples","class_update")),
                          melt(RNA_potential_class,id.var=c("Samples","class_update")))
colnames(combine_potential_class) = c("Samples","Class","State","Count")
combine_potential_class = ddply(combine_potential_class,.(State,Class),transform,Sample.Proportion = Samples/(length(Samples)-1))
combine_potential_class = merge(combine_potential_class,combine_potential,by=c("State","Samples","Sample.Proportion"),all.x=TRUE)
colnames(combine_potential_class)[5:6] = c("Count","Total")
combine_potential_class$TEs.Proportion = combine_potential_class$Count/combine_potential_class$Total

# Proportion of TEs ever in each state that belongs to each TE class
combine_class_ever = ddply(combine_potential_class,.(Class,State),summarise,Proportion=sum(Count[which(Samples > 0)])/sum(Total[which(Samples > 0)]))

rm(list=c("potential_TE_DNase_class","potential_TE_H3K27ac_class","RNA_potential_class",
          "chromHMM_TE_state_class_stats","TE_meth_average_class_stats","potential_TE_DNase_class_stats","potential_TE_H3K27ac_class_stats","RNA_potential_class_stats",
          "class_state_sample","TE_meth_average_state_class","WGBS_sample_state_class","TE_DNase_peaks_class","TE_H3K27ac_peaks_class","RNA_RPKM_class"))
