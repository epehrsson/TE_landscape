# Average number of samples in state by Refseq feature overlap
# Currently with IMR90 for methylation
# See 9/20/2016, 9/28/2016, 12/14/2016, 2/9/2017, 6/6/2017, 6/14/17, 7/24/17, 8/1/2017

source("R_scripts/TE_correlation.R")
source("R_scripts/potential.R")

feature_state_mean = apply(rmsk_TE_measure[,14:21],2,function(y) colMeans(rmsk_TE_measure[which(y != "NA"),c(23:37,40:43,48:49)],na.rm=TRUE))
rownames(feature_state_mean)[1:15] = chromHMM_states
feature_state_mean[1:15,] = feature_state_mean[1:15,]/1.27
feature_state_mean[16:19,] = feature_state_mean[16:19,]/0.37
feature_state_mean[20,] = feature_state_mean[20,]/0.53
feature_state_mean[21,] = feature_state_mean[21,]/0.98

combine_stats = rbind(chromHMM_TE_state_dist_stats,TE_meth_average_category_stats,TE_DNase_potential_stats,TE_H3K27ac_potential_stats)
feature_state_mean = as.data.frame(cbind(combine_stats[match(rownames(feature_state_mean),combine_stats$State),]$Samples_avg_all,feature_state_mean))
colnames(feature_state_mean)[1] = "All"