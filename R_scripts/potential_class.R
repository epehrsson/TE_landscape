# Potential statistics by class
# See 4/27/2016, 5/20/2016, 6/27/2016, 9/17/2016, 2/3/2017, 2/6/2017, 3/2/2017, 5/18/2017, 6/5/2017, 7/4/2017
# See 5/9/2016, 6/2/2016, 6/27/2016, 8/24/2016, 9/7/2016, 9/17/2016, 9/28/2016, 11/27/2016, 12/13/2016, 12/15/2016, 1/13/2017, 2/6/2017, 3/2/2017, 3/3/2017, 5/14/2017, 5/18/2017, 7/21/2017, 7/24/2017, 8/1/2017

# chromHMM
# Distribution of state, all TEs, by class
potential_TEother_state_class = by(potential_TEother_state,potential_TEother_state$class,function(x) sample_distribution(x,c(8:22),127))
potential_TEother_state_class$Unconfident = potential_TEother_state_class[[5]] + potential_TEother_state_class[[6]] + potential_TEother_state_class[[7]] + potential_TEother_state_class[[10]] + potential_TEother_state_class[[11]] + potential_TEother_state_class[[12]]
potential_TEother_state_class$Unconfident$Samples = potential_TEother_state_class$Unknown$Samples

# Cumulative distribution of state, all TEs, by class
potential_TEother_state_class_cum = by(potential_TEother_state,potential_TEother_state$class,function(x) cumulative_distribution(x,c(8:22),127))
potential_TEother_state_class_cum_long = rbind(melt(as.matrix(potential_TEother_state_class_cum$DNA)),melt(as.matrix(potential_TEother_state_class_cum$LINE)),melt(as.matrix(potential_TEother_state_class_cum$LTR)),melt(as.matrix(potential_TEother_state_class_cum$SINE)),melt(as.matrix(potential_TEother_state_class_cum$DNA?)),melt(as.matrix(potential_TEother_state_class_cum$LINE?)),melt(as.matrix(potential_TEother_state_class_cum$LTR?)),melt(as.matrix(potential_TEother_state_class_cum$SINE?)),melt(as.matrix(potential_TEother_state_class_cum$Unknown)),melt(as.matrix(potential_TEother_state_class_cum$Unknown?)),melt(as.matrix(potential_TEother_state_class_cum$Other)),melt(as.matrix(potential_TEother_state_class_cum$RC)))
potential_TEother_state_class_cum_long$Class = c(rep("DNA",1905),rep("LINE",1905),rep("LTR",1905),rep("SINE",1905),rep("DNA?",1905),rep("LINE?",1905),rep("LTR?",1905),rep("SINE?",1905),rep("Unknown",1905),rep("Unknown?",1905),rep("Ohter",1905),rep("RC",1905))
colnames(potential_TEother_state_class_cum_long)[1:3] = c("Samples","State","Proportion")
potential_TEother_state_class_cum_long$Group = paste(potential_TEother_state_class_cum_long$State,potential_TEother_state_class_cum_long$Class,sep="_")

# Average number of samples all TEs are in a state, number of all TEs ever in state, by class
potential_TEother_state_class_stats = lapply(potential_TEother_state_class,function(y) rbind(apply(y,2,function(x) sum(x[2:128])/(sum(x)/100)),apply(y,2,function(x) sum(x*as.numeric(y$Samples))/sum(x))/1.27,apply(y[2:128,],2,function(x) sum(x*as.numeric(y$Samples[2:128]))/sum(x))/1.27))

# Number of all TEs ever in state by class
potential_TEother_state_class_stats_ever = ldply(potential_TEother_state_class_stats[c(1:4,8:9,13)],function(x) as.data.frame(cbind(chromHMM_states,x[1,2:16])))
colnames(potential_TEother_state_class_stats_ever) = c("Class","State","Proportion")
potential_TEother_state_class_stats_ever$State = factor(potential_TEother_state_class_stats_ever$State,levels=chromHMM_states)
potential_TEother_state_class_stats_ever$Proportion = as.numeric(as.character(potential_TEother_state_class_stats_ever$Proportion))
potential_TEother_state_class_stats_ever$Class = factor(potential_TEother_state_class_stats_ever$Class,levels=c("DNA","LINE","LTR","SINE","Other","RC","Unconfident"))

# Average number of samples all TEs are in a state by class
potential_TEother_state_class_stats_avg = ldply(potential_TEother_state_class_stats[c(1:4,8:9,13)], function(x) as.data.frame(cbind(chromHMM_states,x[2,2:16])))
colnames(potential_TEother_state_class_stats_avg) = c("Class","State","Samples")
potential_TEother_state_class_stats_avg$Samples = as.numeric(as.character(potential_TEother_state_class_stats_avg$Samples))
potential_TEother_state_class_stats_avg$State = factor(potential_TEother_state_class_stats_avg$State,levels=chromHMM_states[c(1:3,6:7,10:12,4:5,8:9,13:15)])
potential_TEother_state_class_stats_avg$Class = factor(potential_TEother_state_class_stats_avg$Class,levels=c("DNA","LINE","LTR","SINE","Other","RC","Unconfident"))

# Average number of samples all TEs in a state remain in the state by class
potential_TEother_state_class_stats_cum = ldply(potential_TEother_state_class_stats[c(1:4,8:9,13)],function(x) as.data.frame(cbind(chromHMM_states,x[3,2:16])))
colnames(potential_TEother_state_class_stats_cum) = c("Class","State","Samples")
potential_TEother_state_class_stats_cum$State = factor(potential_TEother_state_class_stats_cum$State,levels=chromHMM_states)
potential_TEother_state_class_stats_cum$Samples = as.numeric(as.character(potential_TEother_state_class_stats_cum$Samples))
potential_TEother_state_class_stats_cum$Class = factor(potential_TEother_state_class_stats_cum$Class,levels=c("DNA","LINE","LTR","SINE","Other","RC","Unconfident"))

# WGBS
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

# DNase
# Distribution of TEs overlapping DNase peaks, by class
potential_TEother_DNase_class = as.matrix(table(TE_DNase_peaks[,61:62]))

# Adding TEs never overlapping DNase peak
potential_TEother_DNase_class = rbind(rmsk_TEother_class_pi[match(colnames(potential_TEother_DNase_class),rmsk_TEother_class_pi$Class),]$Elements - colSums(potential_TEother_DNase_class),potential_TEother_DNase_class)
rownames(potential_TEother_DNase_class) = seq(0,53,1)
potential_TEother_DNase_class = as.data.frame(potential_TEother_DNase_class)
potential_TEother_DNase_class$Total = rowSums(potential_TEother_DNase_class)

# Statistics
potential_TEother_DNase_class_stats = rbind(apply(potential_TEother_DNase_class,2,function(x) sum(x[2:54])/(sum(x)/100)),apply(potential_TEother_DNase_class,2,function(x) sum(x*as.numeric(rownames(potential_TEother_DNase_class)))/sum(x))/0.53,apply(potential_TEother_DNase_class[2:54,],2,function(x) sum(x*as.numeric(rownames(potential_TEother_DNase_class)[2:54]))/sum(x))/0.53)

# H3K27ac
# Distribution of TEs overlapping H3K27ac peaks, by class
potential_TEother_H3K27ac_class = as.matrix(table(TE_H3K27ac_peaks[,106:107]))

# Adding TEs never overlapping H3K27ac peak
potential_TEother_H3K27ac_class = rbind(rmsk_TEother_class_pi[match(colnames(potential_TEother_H3K27ac_class),rmsk_TEother_class_pi$Class),]$Elements - 
                                          colSums(potential_TEother_H3K27ac_class),potential_TEother_H3K27ac_class)
rownames(potential_TEother_H3K27ac_class) = seq(0,98,1)
potential_TEother_H3K27ac_class = as.data.frame(potential_TEother_H3K27ac_class)
potential_TEother_H3K27ac_class$Total = rowSums(potential_TEother_H3K27ac_class)

# Statistics
potential_TEother_H3K27ac_class_stats = rbind(apply(potential_TEother_H3K27ac_class,2,function(x) sum(x[2:99])/(sum(x)/100)),apply(potential_TEother_H3K27ac_class,2,function(x) sum(x*as.numeric(rownames(potential_TEother_H3K27ac_class)))/sum(x))/0.98,apply(potential_TEother_H3K27ac_class[2:99,],2,function(x) sum(x*as.numeric(rownames(potential_TEother_H3K27ac_class)[2:99]))/sum(x))/0.98)