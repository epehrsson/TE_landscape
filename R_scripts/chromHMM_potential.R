# Potential for chromHMM
# See 4/25/2016, 4/26/2016, 4/27/2016, 5/3/2016, 5/11/2016, 5/12/2016, 6/27/2016, 8/30/2016, 9/20/2016, 9/28/2016, 2/3/2017, 2/6/2017, 5/10/2017, 6/14/2017, 8/1/2017, 8/2/2017

# Load chromHMM matrix
load("xx")

# All samples
# TEs ever in composite states
sum(apply(potential_TEother_state[8:15],1,function(x) sum(x > 0)) > 0)
sum(apply(potential_TEother_state[16:22],1,function(x) sum(x > 0)) > 0)
sum(apply(potential_TEother_state[c(8:10,13:14)],1,function(x) sum(x > 0)) > 0)
sum(apply(potential_TEother_state[11:12],1,function(x) sum(x > 0)) > 0)
sum(apply(potential_TEother_state[17:19],1,function(x) sum(x > 0)) > 0)
sum(apply(potential_TEother_state[c(16,20,21)],1,function(x) sum(x > 0)) > 0)

# Distribution
potential_TEother_state_dist = sample_distribution(potential_TEother_state,c(8:22),127)

# Cumulative distribution
potential_TEother_state_cum = cumulative_distribution(potential_TEother_state,c(8:22),127)
potential_TEother_state_cum_long = melt(as.matrix(potential_TEother_state_cum))
colnames(potential_TEother_state_cum_long) = c("Samples","State","Proportion")

# Stats
potential_TEother_state_dist_stats = rbind(apply(potential_TEother_state_dist,2,function(x) sum(x[2:128])/44307.88),apply(potential_TEother_state_dist,2,function(x) sum(x*as.numeric(potential_TEother_state_dist$Samples))/sum(x))/1.27,apply(potential_TEother_state_dist[2:128,],2,function(x) sum(x*as.numeric(potential_TEother_state_dist$Samples[2:128]))/sum(x))/1.27)
potential_TEother_state_dist_stats_ever = as.data.frame(cbind(chromHMM_states,potential_TEother_state_dist_stats[1,2:16]))
colnames(potential_TEother_state_dist_stats_ever) = c("State","Proportion")
potential_TEother_state_dist_stats_ever$State = factor(potential_TEother_state_dist_stats_ever$State,levels=chromHMM_states)
potential_TEother_state_dist_stats_ever$Proportion = as.numeric(as.character(potential_TEother_state_dist_stats_ever$Proportion))
potential_TEother_state_dist_stats_avg = as.data.frame(cbind(chromHMM_states,potential_TEother_state_dist_stats[2,2:16]))
colnames(potential_TEother_state_dist_stats_avg) = c("State","Samples")
potential_TEother_state_dist_stats_avg$Samples = as.numeric(as.character(potential_TEother_state_dist_stats_avg$Samples))
potential_TEother_state_dist_stats_avg$State = factor(potential_TEother_state_dist_stats_avg$State,levels=chromHMM_states[c(1:3,6:7,10:12,4:5,8:9,13:15)])

# Number of TEs in each chromHMM state by sample
state_sample_count = read.table("state_sample_counts.txt",sep='\t')
colnames(state_sample_count) = c("State","Sample","Count")
state_sample_count[1905,] = c("3_TxFlnk","E002",0)
state_sample_count$Count = as.numeric(state_sample_count$Count)
state_sample_count$Proportion = state_sample_count$Count/4430788
state_sample_count[which(state_sample_count$Sample %in% c("E116","E117","E123","E124","E126","E127")),]$Proportion = state_sample_count[which(state_sample_count$Sample %in% c("E116","E117","E123","E124","E126","E127")),]$Count/(4430788-31580)

# No cancer cell lines
# Load chromHMM matrix
load("xx")

# Distribution
potential_TEother_state_dist_noCancer = sample_distribution(potential_TEother_state_noCancer,c(8:22),121)

# Cumulative distribution
potential_TEother_state_noCancer_cum = cumulative_distribution(potential_TEother_state_noCancer,c(8:22),121)
potential_TEother_state_noCancer_cum_long = melt(as.matrix(potential_TEother_state_noCancer_cum))
colnames(potential_TEother_state_noCancer_cum_long) = c("Samples","State","Proportion")

# Stats
potential_TEother_state_dist_noCancer_stats = rbind(apply(potential_TEother_state_dist_noCancer,2,function(x) sum(x[2:122])/44307.88),apply(potential_TEother_state_dist_noCancer,2,function(x) sum(x*as.numeric(potential_TEother_state_dist_noCancer$Samples))/sum(x))/1.21,apply(potential_TEother_state_dist_noCancer[2:122,],2,function(x) sum(x*as.numeric(potential_TEother_state_dist_noCancer$Samples[2:122]))/sum(x))/1.21)
potential_TEother_state_dist_noCancer_stats_ever = as.data.frame(cbind(chromHMM_states,potential_TEother_state_dist_noCancer_stats[1,2:16]))
colnames(potential_TEother_state_dist_noCancer_stats_ever) = c("State","Proportion")
potential_TEother_state_dist_noCancer_stats_ever$State = factor(potential_TEother_state_dist_noCancer_stats_ever$State,levels=chromHMM_states)
potential_TEother_state_dist_noCancer_stats_ever$Proportion = as.numeric(as.character(potential_TEother_state_dist_noCancer_stats_ever$Proportion))
potential_TEother_state_dist_noCancer_stats_avg = as.data.frame(cbind(chromHMM_states,potential_TEother_state_dist_noCancer_stats[2,2:16]))
colnames(potential_TEother_state_dist_noCancer_stats_avg) = c("State","Samples")
potential_TEother_state_dist_noCancer_stats_avg$Samples = as.numeric(as.character(potential_TEother_state_dist_noCancer_stats_avg$Samples))
potential_TEother_state_dist_noCancer_stats_avg$State = factor(potential_TEother_state_dist_noCancer_stats_avg$State,levels=chromHMM_states[c(1:3,6:7,10:12,4:5,8:9,13:15)])