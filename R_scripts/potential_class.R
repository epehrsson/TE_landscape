# Potential statistics by class
# See 4/27/2016, 5/20/2016, 6/27/2016, 9/17/2016, 2/3/2017, 2/6/2017, 3/2/2017, 5/18/2017, 6/5/2017, 7/4/2017
# See 5/9/2016, 6/2/2016, 6/27/2016, 8/24/2016, 9/7/2016, 9/17/2016, 9/28/2016, 11/27/2016, 12/13/2016, 12/15/2016, 1/13/2017, 2/6/2017, 3/2/2017, 3/3/2017, 5/14/2017, 5/18/2017, 7/21/2017, 7/24/2017, 8/1/2017

# chromHMM

chromHMM_TE_state_class = ddply(chromHMM_TE_state,~class_update,function(x) sample_distribution(x,c(8:22),sample_counts["All","chromHMM"]))

chromHMM_TE_state_class_stats = ddply(chromHMM_TE_state_class,~class_update,function(x) potential_stats(x[,2:17],15,sample_counts["All","chromHMM"]))
chromHMM_TE_state_class_stats$State = factor(rep(chromHMM_states,6),levels=chromHMM_states)
chromHMM_TE_state_class_stats$class_update = factor(chromHMM_TE_state_class_stats$class_update,levels=c("DNA","LINE","LTR","SINE","SVA","Other"))
chromHMM_TE_state_class_stats[,2:4] = apply(chromHMM_TE_state_class_stats[,2:4],2,function(x) as.numeric(x))

# Proportion of TEs in each chromHMM state by class by sample
class_state_sample = read.table("chromHMM/class_state_sample.txt",sep='\t')
colnames(class_state_sample) = c("Class","State","Sample","Count")
class_state_sample = merge(class_state_sample,expand.grid(Class = levels(class_state_sample$Class),Sample = levels(class_state_sample$Sample),State = levels(class_state_sample$State)),by=c("Class","State","Sample"),all.y=TRUE)
class_state_sample[is.na(class_state_sample)] = 0
class_state_sample$Count = as.numeric(class_state_sample$Count)
class_state_sample$Proportion = ifelse(metadata[match(class_state_sample$Sample,metadata$Sample),]$chrY == "Yes",
                                       class_state_sample$Count/as.numeric(rmsk_TE_class[match(class_state_sample$Class,rmsk_TE_class$class_update),]$Count),
                                       class_state_sample$Count/as.numeric(rmsk_TE_class[match(class_state_sample$Class,rmsk_TE_class$class_update),]$Count_noY))
class_state_sample$State = factor(class_state_sample$State,levels=chromHMM_states)                                                                                            

# WGBS

TE_meth_average_class = ddply(TE_meth_average,~class_update,function(x) sample_distribution(x,c(46:49),sample_counts["All","WGBS"]))

TE_meth_average_class_stats = ddply(TE_meth_average_class,~class_update,function(x) potential_stats(x[,2:6],4,sample_counts["All","WGBS"]))
TE_meth_average_class_stats$State = factor(rep(meth_states[c(1,3,2,4)],6),levels=meth_states)
TE_meth_average_class_stats$class_update = factor(TE_meth_average_class_stats$class_update,levels=c("DNA","LINE","LTR","SINE","SVA","Other"))
TE_meth_average_class_stats[,2:4] = apply(TE_meth_average_class_stats[,2:4],2,function(x) as.numeric(x))

# Proportion of TEs in each methylation state by class by sample
TE_meth_average_state_class = ddply(TE_meth_average,~class_update,function(y) 
  as.data.frame(cbind(apply(y[,8:44],2,function(x) sum(na.omit(as.numeric(x)) < 0.3)/length(x)),
                      apply(y[,8:44],2,function(x) sum(na.omit(as.numeric(x)) <= 0.7 & na.omit(as.numeric(x)) >= 0.3)/length(x)),
                      apply(y[,8:44],2,function(x) sum(na.omit(as.numeric(x)) > 0.7)/length(x)),
                      apply(y[,8:44],2,function(x) sum(is.na(x))/length(x)))))
colnames(TE_meth_average_state_class)[2:5] = c("Hypomethylated","Intermediate","Hypermethylated","Missing")
TE_meth_average_state_class$Sample = rep(colnames(TE_meth_average)[8:44],6)
TE_meth_average_state_class = melt(TE_meth_average_state_class)
colnames(TE_meth_average_state_class) = c("Class","Sample","State","Proportion")

# DNase

colnames(TE_DNase_peaks)[62] = "DNase"
potential_TE_DNase_class = ddply(TE_DNase_peaks,~class_update,function(x) sample_distribution(x,62,sample_counts["All","DNase"]))
potential_TE_DNase_class[which(potential_TE_DNase_class$Samples == 0),3] =  rmsk_TE_class$Count - aggregate(data=potential_TE_DNase_class,DNase~class_update,sum)$DNase

potential_TE_DNase_class_stats = ddply(potential_TE_DNase_class,~class_update,function(x) potential_stats(x[,2:3],1,sample_counts["All","DNase"]))
potential_TE_DNase_class_stats$State = rep("DNase",6)

# Proportion of TEs overlapping DNase peak by class by sample
TE_DNase_peaks_class = aggregate(data=TE_DNase_peaks[,c(8:61)],.~class_update,function(x) sum(x > 0))
TE_DNase_peaks_class = melt(TE_DNase_peaks_class,id.var="class_update")
colnames(TE_DNase_peaks_class) = c("Class","Sample","Proportion")
TE_DNase_peaks_class$Proportion = TE_DNase_peaks_class$Count/ifelse(metadata[match(TE_DNase_peaks_class$Sample,metadata$Sample),]$chrY == "Yes",
                                                                    rmsk_TE_class[match(TE_DNase_peaks_class$Class,rmsk_TE_class$class_update),]$Count,
                                                                    rmsk_TE_class[match(TE_DNase_peaks_class$Class,rmsk_TE_class$class_update),]$Count_noY)
TE_DNase_peaks_class$State = rep("DNase",dim(TE_DNase_peaks_class)[1])

colnames(TE_DNase_peaks)[62] = "Samples"

# H3K27ac

colnames(TE_H3K27ac_peaks)[107] = "H3K27ac"
potential_TE_H3K27ac_class = ddply(TE_H3K27ac_peaks,~class_update,function(x) sample_distribution(x,107,sample_counts["All","H3K27ac"]))
potential_TE_H3K27ac_class[which(potential_TE_H3K27ac_class$Samples == 0),3] =  rmsk_TE_class$Count - aggregate(data=potential_TE_H3K27ac_class,H3K27ac~class_update,sum)$H3K27ac

potential_TE_H3K27ac_class_stats = ddply(potential_TE_H3K27ac_class,~class_update,function(x) potential_stats(x[,2:3],1,sample_counts["All","H3K27ac"]))
potential_TE_H3K27ac_class_stats$State = rep("H3K27ac",6)

# Proportion of TEs overlapping H3K27ac peak by class by sample
TE_H3K27ac_peaks_class = aggregate(data=TE_H3K27ac_peaks[,c(8:106)],.~class_update,function(x) sum(x > 0))
TE_H3K27ac_peaks_class = melt(TE_H3K27ac_peaks_class,id.var="class_update")
colnames(TE_H3K27ac_peaks_class) = c("Class","Sample","Proportion")
TE_H3K27ac_peaks_class$Proportion = TE_H3K27ac_peaks_class$Count/ifelse(metadata[match(TE_H3K27ac_peaks_class$Sample,metadata$Sample),]$chrY == "Yes",
                                                                    rmsk_TE_class[match(TE_H3K27ac_peaks_class$Class,rmsk_TE_class$class_update),]$Count,
                                                                    rmsk_TE_class[match(TE_H3K27ac_peaks_class$Class,rmsk_TE_class$class_update),]$Count_noY)
TE_H3K27ac_peaks_class$State = rep("H3K27ac",dim(TE_H3K27ac_peaks_class)[1])

colnames(TE_H3K27ac_peaks)[107] = "Samples"

# Expression

RNA_potential_class = ddply(RNA_TE_agnostic,~class_update,function(x) sample_distribution(x,61,sample_counts["All","RNA"]))

RNA_potential_class_stats = ddply(RNA_potential_class,~class_update,function(x) potential_stats(x[,2:3],1,sample_counts["All","RNA"]))
RNA_potential_class_stats$class_update = factor(RNA_potential_class_stats$class_update,levels=c("DNA","LINE","LTR","SINE","SVA","Other"))
RNA_potential_class_stats$State = rep("Expression",6)

# Proportion of TEs expressed RPKM > 1 by class by sample
RNA_RPKM_class = aggregate(data=RNA_TE_agnostic[,c(9:60,64)],.~class_update,function(x) sum(x > 1))
RNA_RPKM_class = melt(RNA_RPKM_class,id.var="class_update")
colnames(RNA_RPKM_class) = c("Class","Sample","Proportion")
RNA_RPKM_class$Proportion = RNA_RPKM_class$Count/ifelse(metadata[match(RNA_RPKM_class$Sample,metadata$Sample),]$chrY == "Yes",
                                                                        rmsk_TE_class[match(RNA_RPKM_class$Class,rmsk_TE_class$class_update),]$Count,
                                                                        rmsk_TE_class[match(RNA_RPKM_class$Class,rmsk_TE_class$class_update),]$Count_noY)
RNA_RPKM_class$State = rep("Expression",dim(RNA_RPKM_class)[1])
