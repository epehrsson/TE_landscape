# Potential for promoters

load("R_datasets/promoter_matrices.RData")

# chromHMM potential
chromHMM_promoter_state_dist = sample_distribution(chromHMM_promoter_state,c(5:19),127)
chromHMM_promoter_state_cum = cumulative_distribution(chromHMM_promoter_state,c(5:19),127)
chromHMM_promoter_state_cum$State = factor(rep(chromHMM_states,each=127),levels=chromHMM_states)
chromHMM_promoter_state_dist_stats = potential_stats(chromHMM_promoter_state_dist,15,127)
chromHMM_promoter_state_dist_stats$State = factor(chromHMM_states,levels=chromHMM_states)

# WGBS potential
# Cumulative distribution of methylation states and statistics
promoter_meth_average_category = sample_distribution(promoter_meth_average,c(42:45),37)
promoter_meth_average_category_cum = cumulative_distribution(promoter_meth_average,c(42:45),37)
promoter_meth_average_category_stats = potential_stats(promoter_meth_average_category,4,37)
promoter_meth_average_category_stats$State = factor(rownames(promoter_meth_average_category_stats),levels=meth_states)

# Cumulative distribution of methylation states and statistics, no IMR90
promoter_meth_average_noIMR90_category = sample_distribution(promoter_meth_average,c(46:49),36)
promoter_meth_average_noIMR90_category_cum = cumulative_distribution(promoter_meth_average,c(46:49),36)
promoter_meth_average_noIMR90_category_stats = potential_stats(promoter_meth_average_noIMR90_category,4,36)
promoter_meth_average_noIMR90_category_stats$State = factor(c("Hypomethylated","Hypermethylated","Intermediate","Missing"),levels=meth_states)

# DNase potential
promoter_DNase_potential = sample_distribution(promoter_DNase_peaks,58,53)
colnames(promoter_DNase_potential)[2] = "DNase"
promoter_DNase_potential_cum = cumulative_distribution(promoter_DNase_peaks,58,53)
promoter_DNase_potential_cum$State = rep("DNase",53)
promoter_DNase_potential_stats = potential_stats(promoter_DNase_potential,1,53)
promoter_DNase_potential_stats$State = "DNase"

# H3K27ac potential
promoter_H3K27ac_potential = sample_distribution(promoter_H3K27ac_peaks,103,98)
colnames(promoter_H3K27ac_potential)[2] = "H3K27ac"
promoter_H3K27ac_potential_cum = cumulative_distribution(promoter_H3K27ac_peaks,103,98)
promoter_H3K27ac_potential_cum$State = rep("H3K27ac",98)
promoter_H3K27ac_potential_stats = potential_stats(promoter_H3K27ac_potential,1,98)
promoter_H3K27ac_potential_stats$State = "H3K27ac"

# Promoters in state by sample
# Number of promoters in each chromHMM state by sample
promoter_state_sample_count = read.table("chromHMM/Refseq_promoters/promoter_state_sample_count.txt",sep='\t')
colnames(promoter_state_sample_count) = c("State","Sample","Count")
promoter_state_sample_count$Proportion = promoter_state_sample_count$Count/34750
promoter_state_sample_count[which(promoter_state_sample_count$Sample %in% c("E116","E117","E123","E124","E126","E127")),]$Proportion = promoter_state_sample_count[which(promoter_state_sample_count$Sample %in% c("E116","E117","E123","E124","E126","E127")),]$Count/34568
promoter_state_sample_count$State = factor(promoter_state_sample_count$State,levels=chromHMM_states)

# Number of promoters in each WGBS state by sample
promoter_meth_average_state = as.data.frame(cbind(apply(promoter_meth_average[,7:41],2,function(x) sum(na.omit(x) < 0.3)/length(x)),apply(promoter_meth_average[,7:41],2,function(x)  sum(na.omit(x) <= 0.7 & na.omit(x) >= 0.3)/length(x)),apply(promoter_meth_average[,7:41],2,function(x) sum(na.omit(x) > 0.7)/length(x)),apply(promoter_meth_average[,7:41],2,function(x) sum(is.na(x))/length(x))))
colnames(promoter_meth_average_state) = c("Hypomethylated","Intermediate","Hypermethylated","Missing")
promoter_meth_average_state = promoter_meth_average_state[order(promoter_meth_average_state$Hypermethylated + promoter_meth_average_state$Missing),]

promoter_meth_average_state_long = melt(as.matrix(promoter_meth_average_state))
colnames(promoter_meth_average_state_long) = c("Sample","State","Proportion")

# Proportion of promoters overlapping Dnase peaks, by sample
promoter_DNase_peaks_sample = as.data.frame(apply(promoter_DNase_peaks[,5:57],2,function(x) sum(x > 0)))
colnames(promoter_DNase_peaks_sample) = "Count"
promoter_DNase_peaks_sample$Proportion = promoter_DNase_peaks_sample$Count/34750
promoter_DNase_peaks_sample[which(rownames(promoter_DNase_peaks_sample) %in% c("E116","E117","E123","E124","E126","E127")),]$Proportion = promoter_DNase_peaks_sample[which(rownames(promoter_DNase_peaks_sample) %in% c("E116","E117","E123","E124","E126","E127")),]$Count/34568
promoter_DNase_peaks_sample$Sample = rownames(promoter_DNase_peaks_sample)
promoter_DNase_peaks_sample$State = rep("DNase",53)

# Proportion of promoters overlapping H3K27ac peaks, by sample
promoter_H3K27ac_peaks_sample = as.data.frame(apply(promoter_H3K27ac_peaks[,5:102],2,function(x) sum(x > 0)))
colnames(promoter_H3K27ac_peaks_sample) = "Count"
promoter_H3K27ac_peaks_sample$Proportion = promoter_H3K27ac_peaks_sample$Count/34750
promoter_H3K27ac_peaks_sample[which(rownames(promoter_H3K27ac_peaks_sample) %in% c("E116","E117","E123","E124","E126","E127")),]$Proportion = promoter_H3K27ac_peaks_sample[which(rownames(promoter_H3K27ac_peaks_sample) %in% c("E116","E117","E123","E124","E126","E127")),]$Count/34568
promoter_H3K27ac_peaks_sample$Sample = rownames(promoter_H3K27ac_peaks_sample)
promoter_H3K27ac_peaks_sample$State = rep("H3K27ac",98)