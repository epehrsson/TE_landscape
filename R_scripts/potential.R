# Potential 
# See 4/25/2016, 4/26/2016, 4/27/2016, 5/3/2016, 5/11/2016, 5/12/2016, 6/2/2016, 6/27/2016, 7/11/2016, 8/25/2016, 8/26/2016, 8/30/2016, 9/16/2016, 9/20/2016, 9/28/2016, 12/15/2016, 2/3/2017, 2/6/2017, 5/10/2017, 6/14/2017, 7/22/2017, 8/1/2017, 8/2/2017, 8/3/2017

# chromHMM potential
source("~/TE_landscape/R_scripts/chromHMM_matrix.R")

# All samples
# Distribution
chromHMM_TE_state_dist = sample_distribution(chromHMM_TE_state,c(8:22),127)

# Cumulative distribution
chromHMM_TE_state_cum = cumulative_distribution(chromHMM_TE_state,c(8:22),127)
chromHMM_TE_state_cum_long = melt(as.matrix(chromHMM_TE_state_cum))
colnames(chromHMM_TE_state_cum_long) = c("Samples","State","Proportion")

# Stats
chromHMM_TE_state_dist_stats = as.data.frame(t(rbind(apply(chromHMM_TE_state_dist,2,function(x) sum(x[2:128])/44307.88),apply(chromHMM_TE_state_dist,2,function(x) sum(x*as.numeric(chromHMM_TE_state_dist$Samples))/sum(x))/1.27,apply(chromHMM_TE_state_dist[2:128,],2,function(x) sum(x*as.numeric(chromHMM_TE_state_dist$Samples[2:128]))/sum(x))/1.27)))
colnames(chromHMM_TE_state_dist_stats)[1:3] = c("Proportion_ever","Samples_avg_all","Samples_avg_ever")
chromHMM_TE_state_dist_stats$State = factor(chromHMM_TE_state_dist_stats$State,levels=chromHMM_states)

# No cancer cell lines
# Distribution
chromHMM_TE_state_dist_noCancer = sample_distribution(chromHMM_TE_state_noCancer,c(8:22),121)

# Cumulative distribution
chromHMM_TE_state_noCancer_cum = cumulative_distribution(chromHMM_TE_state_noCancer,c(8:22),121)
chromHMM_TE_state_noCancer_cum_long = melt(as.matrix(chromHMM_TE_state_noCancer_cum))
colnames(chromHMM_TE_state_noCancer_cum_long) = c("Samples","State","Proportion")

# Stats
chromHMM_TE_state_noCancer_dist_stats = as.data.frame(t(rbind(apply(chromHMM_TE_state_noCancer_dist,2,function(x) sum(x[2:122])/44307.88),apply(chromHMM_TE_state_noCancer_dist,2,function(x) sum(x*as.numeric(chromHMM_TE_state_noCancer_dist$Samples))/sum(x))/1.21,apply(chromHMM_TE_state_noCancer_dist[2:122,],2,function(x) sum(x*as.numeric(chromHMM_TE_state_noCancer_dist$Samples[2:122]))/sum(x))/1.21)))
colnames(chromHMM_TE_state_noCancer_dist_stats)[1:3] = c("Proportion_ever","Samples_avg_all","Samples_avg_ever")
chromHMM_TE_state_noCancer_dist_stats$State = factor(chromHMM_TE_state_noCancer_dist_stats$State,levels=chromHMM_states)

# Number of TEs in each chromHMM state by sample
state_sample_count = read.table("state_sample_counts.txt",sep='\t')
colnames(state_sample_count) = c("State","Sample","Count")
state_sample_count[1905,] = c("3_TxFlnk","E002",0)
state_sample_count$Count = as.numeric(state_sample_count$Count)
state_sample_count$Proportion = state_sample_count$Count/4430788
state_sample_count[which(state_sample_count$Sample %in% c("E116","E117","E123","E124","E126","E127")),]$Proportion = state_sample_count[which(state_sample_count$Sample %in% c("E116","E117","E123","E124","E126","E127")),]$Count/(4430788-31580)

# WGBS potential
source("~/TE_landscape/R_scripts/WGBS_TE_avg_methylation.R")

# Cumulative distribution of methylation states and statistics
TE_meth_average_category = sample_distribution(TE_meth_average,c(46:49),37)
TE_meth_average_category_cum = cumulative_distribution(TE_meth_average,c(46:49),37)
TE_meth_average_category_cum_long = melt(as.matrix(TE_meth_average_category_cum))
colnames(TE_meth_average_category_cum_long) = c("Samples","Methylation","Proportion")
TE_meth_average_category_stats = as.data.frame(t(rbind(apply(TE_meth_average_category[,2:5],2,function(x) sum(x[2:38])/3200428),apply(TE_meth_average_category[,2:5],2,function(x) sum(as.numeric(x)*seq(0,37))/sum(x))/37,apply(TE_meth_average_category[2:38,2:5],2,function(x) sum(as.numeric(x)*seq(1,37))/sum(x))/37)))
TE_meth_average_category_stats$State = rownames(TE_meth_average_category_stats)
colnames(TE_meth_average_category_stats)[1:3] = c("Proportion_ever","Samples_avg_all","Samples_avg_ever")
TE_meth_average_category_stats$State = factor(TE_meth_average_category_stats$State,levels=as.vector(TE_meth_average_category_stats$State)[c(4,2,3,1)])

# Cumulative distribution of methylation states and statistics, no IMR90
TE_meth_average_noIMR90_category = sample_distribution(TE_meth_average,c(50:53),36)
TE_meth_average_noIMR90_category_cum = cumulative_distribution(TE_meth_average,c(50:53),36)
TE_meth_average_noIMR90_category_cum_long = melt(as.matrix(TE_meth_average_noIMR90_category_cum))
colnames(TE_meth_average_noIMR90_category_cum_long) = c("Samples","Methylation","Proportion")
TE_meth_average_noIMR90_category_stats = as.data.frame(t(rbind(apply(TE_meth_average_noIMR90_category[,2:5],2,function(x) sum(x[2:37])/3200428),apply(TE_meth_average_noIMR90_category[,2:5],2,function(x) sum(as.numeric(x)*seq(0,36))/sum(x))/36,apply(TE_meth_average_noIMR90_category[2:37,2:5],2,function(x) sum(as.numeric(x)*seq(1,36))/sum(x))/36)))
TE_meth_average_noIMR90_category_stats$State = factor(c("Hypomethylated","Hypermethylated","Intermediate","Missing"),levels=c("Missing","Hypermethylated","Intermediate","Hypomethylated"))
colnames(TE_meth_average_noIMR90_category_stats)[1:3] = c("Proportion_ever","Samples_avg_all","Samples_avg_ever")

# DNase potential
source("~/TE_landscape/R_scripts/DNase_peaks.R")

# Distribution of TEs overlapping DNase peaks (needs matrix)
TE_DNase_potential = sample_distribution(TE_DNase_peaks,61,53)
TE_DNase_potential[1,2] = 4430788-1840984
colnames(TE_DNase_potential)[2] = "DNase"

# Cumulative distribution of TEs overlapping DNase peaks
TE_DNase_potential_cum = as.data.frame(cumsum(TE_DNase_potential[which(TE_DNase_potential$Samples %in% seq(1,53)),]$DNase)/sum(TE_DNase_potential[which(TE_DNase_potential$Samples %in% seq(1,53)),]$DNase))
rownames(TE_DNase_potential_cum) = seq(1:53)
TE_DNase_potential_cum$Samples = rownames(TE_DNase_potential_cum)
colnames(TE_DNase_potential_cum)[1] = "Proportion"

# Statistics
TE_DNase_potential_stats = as.data.frame(t(rbind(sum(TE_DNase_potential$DNase[2:54])/sum(TE_DNase_potential$DNase),(sum(as.numeric(TE_DNase_potential$DNase)*seq(0,53))/sum(TE_DNase_potential$DNase))/53,(sum(as.numeric(TE_DNase_potential$DNase[2:54])*seq(1,53))/sum(TE_DNase_potential$DNase))/53)))
colnames(TE_DNase_potential_stats)[1:3] = c("Proportion_ever","Samples_avg_all","Samples_avg_ever")
TE_DNase_potential_stats$State = "yes"

# Proportion of TEs overlapping Dnase peaks, by sample
TE_DNase_peaks_sample = as.data.frame(apply(TE_DNase_peaks[,8:60],2,function(x) sum(x > 0)))
colnames(TE_DNase_peaks_sample) = "DNase"
TE_DNase_peaks_sample$Proportion = TE_DNase_peaks_sample$DNase/4430788
TE_DNase_peaks_sample[which(rownames(TE_DNase_peaks_sample) %in% c("E116","E117","E123","E124","E126","E127")),]$Proportion = TE_DNase_peaks_sample[which(rownames(TE_DNase_peaks_sample) %in% c("E116","E117","E123","E124","E126","E127")),]$DNase/4399208

# H3K27ac potential
source("~/TE_landscape/R_scripts/H3K27ac_peaks.R")

# Distribution of TEs overlapping H3K27ac peaks
TE_H3K27ac_potential = sample_distribution(TE_H3K27ac_peaks,106,98)
TE_H3K27ac_potential[1,2] = 4430788-1801527
colnames(TE_H3K27ac_potential)[2] = "H3K27ac"

# Cumulative distribution of TEs overlapping H3K27ac peaks
TE_H3K27ac_potential_cum = as.data.frame(cumsum(TE_H3K27ac_potential[which(TE_H3K27ac_potential$Samples %in% seq(1,98)),]$H3K27ac)/sum(TE_H3K27ac_potential[which(TE_H3K27ac_potential$Samples %in% seq(1,98)),]$H3K27ac))
rownames(TE_H3K27ac_potential_cum) = seq(1:98)
TE_H3K27ac_potential_cum$Samples = rownames(TE_H3K27ac_potential_cum)
colnames(TE_H3K27ac_potential_cum)[1] = "Proportion"

# Statistics
TE_H3K27ac_potential_stats = as.data.frame(t(rbind(sum(TE_H3K27ac_potential$H3K27ac[2:99])/sum(TE_H3K27ac_potential$H3K27ac),(sum(as.numeric(TE_H3K27ac_potential$H3K27ac)*seq(0,98))/sum(TE_H3K27ac_potential$H3K27ac))/98,(sum(as.numeric(TE_H3K27ac_potential$H3K27ac[2:99])*seq(1,98))/sum(TE_H3K27ac_potential$H3K27ac))/98)))
colnames(TE_H3K27ac_potential_stats)[1:3] = c("Proportion_ever","Samples_avg_all","Samples_avg_ever")
TE_H3K27ac_potential_stats$State = "yes"

# Proportion of TEs overlapping H3K27ac peaks, by sample
TE_H3K27ac_peaks_sample = as.data.frame(apply(TE_H3K27ac_peaks[,8:105],2,function(x) sum(x > 0)))
colnames(TE_H3K27ac_peaks_sample) = "H3K27ac"
TE_H3K27ac_peaks_sample$Proportion = TE_H3K27ac_peaks_sample$H3K27ac/4430788
TE_H3K27ac_peaks_sample[which(rownames(TE_H3K27ac_peaks_sample) %in% c("E116","E117","E123","E124","E126","E127")),]$Proportion = TE_H3K27ac_peaks_sample[which(rownames(TE_H3K27ac_peaks_sample) %in% c("E116","E117","E123","E124","E126","E127")),]$H3K27ac/4399208
