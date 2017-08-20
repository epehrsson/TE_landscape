# Potential statistics by class
# See 4/27/2016, 5/20/2016, 6/27/2016, 9/17/2016, 2/3/2017, 2/6/2017, 3/2/2017, 5/18/2017, 6/5/2017, 7/4/2017
# See 5/9/2016, 6/2/2016, 6/27/2016, 8/24/2016, 9/7/2016, 9/17/2016, 9/28/2016, 11/27/2016, 12/13/2016, 12/15/2016, 1/13/2017, 2/6/2017, 3/2/2017, 3/3/2017, 5/14/2017, 5/18/2017, 7/21/2017, 7/24/2017, 8/1/2017

library(plyr)

source("R_scripts/TE_class_stats.R")

# chromHMM
load("R_scripts/chromHMM_TE_state.RData")

# Distribution of state, all TEs, by class
chromHMM_TE_state_class = by(chromHMM_TE_state,chromHMM_TE_state$class_update,function(x) sample_distribution(x,c(8:22),127))

# Statistics
chromHMM_TE_state_class_stats = ldply(chromHMM_TE_state_class,function(y) as.data.frame(t(rbind(apply(y[,2:16],2,function(x) sum(x[2:128])/(sum(x)/100)),apply(y[,2:16],2,function(x) sum(x*as.numeric(y$Samples))/sum(x))/1.27,apply(y[2:128,2:16],2,function(x) sum(x*as.numeric(y$Samples[2:128]))/sum(x))/1.27))))
colnames(chromHMM_TE_state_class_stats) = c("Class","Proportion_ever","Samples_avg_all","Samples_avg_ever")
chromHMM_TE_state_class_stats$Class = factor(chromHMM_TE_state_class_stats$Class,levels=c("DNA","LINE","LTR","SINE","SVA","Other"))
chromHMM_TE_state_class_stats$State = factor(rep(chromHMM_states,6),levels=chromHMM_states)
chromHMM_TE_state_class_stats[,2:4] = apply(chromHMM_TE_state_class_stats[,2:4],2,function(x) as.numeric(x))
                                                                                                    
# WGBS
load("R_scripts/TE_meth_average.RData")

# Cumulative distribution of methylation states and statistics, by class
TE_meth_average_class = by(TE_meth_average,TE_meth_average$class_update,function(x) sample_distribution(x,c(46:49),37))

TE_meth_average_class_stats = ldply(TE_meth_average_class,function(y) as.data.frame(t(rbind(apply(y[,2:5],2,function(x) sum(x[2:38])/(sum(x)/100)),apply(y[,2:5],2,function(x) sum(as.numeric(x)*seq(0,37))/sum(x))/0.37,apply(y[2:38,2:5],2,function(x) sum(as.numeric(x)*seq(1,37))/sum(x))/0.37))))
colnames(TE_meth_average_class_stats) = c("Class","Proportion_ever","Samples_avg_all","Samples_avg_ever")
TE_meth_average_class_stats$Class = factor(TE_meth_average_class_stats$Class,levels=c("DNA","LINE","LTR","SINE","SVA","Other"))
TE_meth_average_class_stats$State = factor(rep(meth_states[c(1,3,2,4)],6),levels=meth_states)
TE_meth_average_class_stats[,2:4] = apply(TE_meth_average_class_stats[,2:4],2,function(x) as.numeric(x))

# Cumulative distribution of methylation states and statistics, no IMR90, by class
TE_meth_average_class_noIMR90 = by(TE_meth_average,TE_meth_average$class_update,function(x) sample_distribution(x,c(50:53),36))

TE_meth_average_class_noIMR90_stats = ldply(TE_meth_average_class_noIMR90,function(y) as.data.frame(t(rbind(apply(y[,2:5],2,function(x) sum(x[2:37])/(sum(x)/100)),apply(y[,2:5],2,function(x) sum(as.numeric(x)*seq(0,36))/sum(x))/0.36,apply(y[2:37,2:5],2,function(x) sum(as.numeric(x)*seq(1,36))/sum(x))/0.36))))
colnames(TE_meth_average_class_noIMR90_stats) = c("Class","Proportion_ever","Samples_avg_all","Samples_avg_ever")
TE_meth_average_class_noIMR90_stats$Class = factor(TE_meth_average_class_noIMR90_stats$Class,levels=c("DNA","LINE","LTR","SINE","SVA","Other"))
TE_meth_average_class_noIMR90_stats$State = factor(rep(meth_states[c(1,3,2,4)],6),levels=meth_states)
TE_meth_average_class_noIMR90_stats[,2:4] = apply(TE_meth_average_class_noIMR90_stats[,2:4],2,function(x) as.numeric(x))

# DNase
load("R_scripts/TE_DNase_peaks.RData")

# Distribution of TEs overlapping DNase peaks, by class
potential_TE_DNase_class = t(as.matrix(table(TE_DNase_peaks[,61:62])))

# Adding TEs never overlapping DNase peak
potential_TE_DNase_class = rbind(rmsk_TE_class[match(colnames(potential_TE_DNase_class),rmsk_TE_class$class_update),]$Count - colSums(potential_TE_DNase_class),potential_TE_DNase_class)
rownames(potential_TE_DNase_class) = seq(0,53,1)
potential_TE_DNase_class = as.data.frame(potential_TE_DNase_class)
potential_TE_DNase_class$Total = rowSums(potential_TE_DNase_class)

# Statistics
potential_TE_DNase_class_stats = as.data.frame(t(rbind(apply(potential_TE_DNase_class[,1:6],2,function(x) sum(x[2:54])/(sum(x)/100)),apply(potential_TE_DNase_class[,1:6],2,function(x) sum(x*as.numeric(rownames(potential_TE_DNase_class)))/sum(x))/0.53,apply(potential_TE_DNase_class[2:54,1:6],2,function(x) sum(x*as.numeric(rownames(potential_TE_DNase_class)[2:54]))/sum(x))/0.53)))
colnames(potential_TE_DNase_class_stats) = c("Proportion_ever","Samples_avg_all","Samples_avg_ever")
potential_TE_DNase_class_stats$State = "DNase"
potential_TE_DNase_class_stats$Class = rownames(potential_TE_DNase_class_stats)

# Proportion of TEs overlapping DNase peak by class by sample
TE_DNase_peaks_class = aggregate(data=TE_DNase_peaks[,c(8:61)],.~class_update,function(x) sum(x > 0))
rownames(TE_DNase_peaks_class) = TE_DNase_peaks_class$class_update
TE_DNase_peaks_class = TE_DNase_peaks_class[,2:54]
TE_DNase_peaks_class = TE_DNase_peaks_class/rmsk_TE_class[match(rownames(TE_DNase_peaks_class),rmsk_TE_class$class_update),]$Count

# H3K27ac
load("R_scripts/TE_H3K27ac_peaks.RData")

# Distribution of TEs overlapping H3K27ac peaks, by class
potential_TE_H3K27ac_class = t(as.matrix(table(TE_H3K27ac_peaks[,106:107])))

# Adding TEs never overlapping H3K27ac peak
potential_TE_H3K27ac_class = rbind(rmsk_TE_class[match(colnames(potential_TE_H3K27ac_class),rmsk_TE_class$class_update),]$Count - colSums(potential_TE_H3K27ac_class),potential_TE_H3K27ac_class)
rownames(potential_TE_H3K27ac_class) = seq(0,98,1)
potential_TE_H3K27ac_class = as.data.frame(potential_TE_H3K27ac_class)
potential_TE_H3K27ac_class$Total = rowSums(potential_TE_H3K27ac_class)

# Statistics
potential_TE_H3K27ac_class_stats = as.data.frame(t(rbind(apply(potential_TE_H3K27ac_class[,1:6],2,function(x) sum(x[2:99])/(sum(x)/100)),apply(potential_TE_H3K27ac_class[,1:6],2,function(x) sum(x*as.numeric(rownames(potential_TE_H3K27ac_class)))/sum(x))/0.98,apply(potential_TE_H3K27ac_class[2:99,1:6],2,function(x) sum(x*as.numeric(rownames(potential_TE_H3K27ac_class)[2:99]))/sum(x))/0.98)))
colnames(potential_TE_H3K27ac_class_stats) = c("Proportion_ever","Samples_avg_all","Samples_avg_ever")
potential_TE_H3K27ac_class_stats$State = "H3K27ac"
potential_TE_H3K27ac_class_stats$Class = rownames(potential_TE_H3K27ac_class_stats)

# Proportion of TEs overlapping DNase peak by class by sample
TE_H3K27ac_peaks_class = aggregate(data=TE_H3K27ac_peaks[,c(8:106)],.~class_update,function(x) sum(x > 0))
rownames(TE_H3K27ac_peaks_class) = TE_H3K27ac_peaks_class$class_update
TE_H3K27ac_peaks_class = TE_H3K27ac_peaks_class[,2:99]
TE_H3K27ac_peaks_class = TE_H3K27ac_peaks_class/rmsk_TE_class[match(rownames(TE_H3K27ac_peaks_class),rmsk_TE_class$class_update),]$Count