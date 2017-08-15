# ChromHMM analysis by class
# See 4/27/2016, 5/20/2016, 6/27/2016, 9/17/2016, 2/3/2017, 2/6/2017, 3/2/2017, 5/18/2017

# Potential
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

# Contribution
# Contribution of each TE class to each chromHMM state (with Unconfident)
contribution_class = as.data.frame(t(contribution[c(4:7,11:12,16),]))
contribution_class = contribution_class/rowSums(contribution_class)
contribution_class_long = melt(as.matrix(contribution_class))
colnames(contribution_class_long) = c("State","Class","Proportion")
