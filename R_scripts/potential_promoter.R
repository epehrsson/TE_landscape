# Potential for promoters

load("R_datasets/promoter_matrices.RData")

# Unique promoters
NUM_PROMOTER = dim(chromHMM_promoter_state)[1]
NUM_PROMOTER_noY = dim(chromHMM_promoter_state[which(chromHMM_promoter_state$chromosome != "chrY"),])[1]

# chromHMM potential
chromHMM_promoter_state_dist = sample_distribution(chromHMM_promoter_state,c(5:19),sample_counts["All","chromHMM"])
chromHMM_promoter_state_cum = cumulative_distribution(chromHMM_promoter_state,c(5:19),sample_counts["All","chromHMM"])
chromHMM_promoter_state_cum$State = factor(rep(chromHMM_states,each=sample_counts["All","chromHMM"]),levels=chromHMM_states)
chromHMM_promoter_state_dist_stats = potential_stats(chromHMM_promoter_state_dist,15,sample_counts["All","chromHMM"])
chromHMM_promoter_state_dist_stats$State = factor(chromHMM_states,levels=chromHMM_states)

# WGBS potential
# Cumulative distribution of methylation states and statistics
promoter_meth_average_category = sample_distribution(promoter_meth_average,c(42:45),sample_counts["All","WGBS"])
promoter_meth_average_category_cum = cumulative_distribution(promoter_meth_average,c(42:45),sample_counts["All","WGBS"])
promoter_meth_average_category_stats = potential_stats(promoter_meth_average_category,4,sample_counts["All","WGBS"])
promoter_meth_average_category_stats$State = factor(rownames(promoter_meth_average_category_stats),levels=meth_states)

# DNase potential
promoter_DNase_potential = sample_distribution(promoter_DNase_peaks,58,sample_counts["All","DNase"])
colnames(promoter_DNase_potential)[2] = "DNase"
promoter_DNase_potential_cum = cumulative_distribution(promoter_DNase_peaks,58,sample_counts["All","DNase"])
promoter_DNase_potential_cum$State = rep("DNase",sample_counts["All","DNase"])
promoter_DNase_potential_stats = potential_stats(promoter_DNase_potential,1,sample_counts["All","DNase"])
promoter_DNase_potential_stats$State = "DNase"

# H3K27ac potential
promoter_H3K27ac_potential = sample_distribution(promoter_H3K27ac_peaks,103,sample_counts["All","H3K27ac"])
colnames(promoter_H3K27ac_potential)[2] = "H3K27ac"
promoter_H3K27ac_potential_cum = cumulative_distribution(promoter_H3K27ac_peaks,103,sample_counts["All","H3K27ac"])
promoter_H3K27ac_potential_cum$State = rep("H3K27ac",sample_counts["All","H3K27ac"])
promoter_H3K27ac_potential_stats = potential_stats(promoter_H3K27ac_potential,1,sample_counts["All","H3K27ac"])
promoter_H3K27ac_potential_stats$State = "H3K27ac"

# Promoters in state by sample
# Number of promoters in each chromHMM state by sample
promoter_state_sample_count = read.table("chromHMM/Refseq_promoters/promoter_state_sample_count.txt",sep='\t')
colnames(promoter_state_sample_count) = c("State","Sample","Count")
promoter_state_sample_count$Proportion = ifelse(metadata[match(promoter_state_sample_count$Sample,metadata$Sample),]$chrY == "Yes",promoter_state_sample_count$Count/NUM_PROMOTER,promoter_state_sample_count$Count/NUM_PROMOTER_noY)
promoter_state_sample_count$State = factor(promoter_state_sample_count$State,levels=chromHMM_states)

# Number of promoters in each WGBS state by sample
promoter_meth_average_state = as.data.frame(cbind(apply(promoter_meth_average[,7:41],2,function(x) sum(as.numeric(na.omit(x)) < 0.3)/length(x)),
                                                  apply(promoter_meth_average[,7:41],2,function(x)  sum(as.numeric(na.omit(x)) <= 0.7 & as.numeric(na.omit(x)) >= 0.3)/length(x)),
                                                  apply(promoter_meth_average[,7:41],2,function(x) sum(as.numeric(na.omit(x)) > 0.7)/length(x)),
                                                  apply(promoter_meth_average[,7:41],2,function(x) sum(is.na(x))/length(x))))
colnames(promoter_meth_average_state) = c("Hypomethylated","Intermediate","Hypermethylated","Missing")
promoter_meth_average_state = promoter_meth_average_state[order(promoter_meth_average_state$Hypermethylated + promoter_meth_average_state$Missing),]

promoter_meth_average_state_long = melt(as.matrix(promoter_meth_average_state))
colnames(promoter_meth_average_state_long) = c("Sample","State","Proportion")

# Proportion of promoters overlapping Dnase peaks, by sample
promoter_DNase_peaks_sample = as.data.frame(apply(promoter_DNase_peaks[,5:57],2,function(x) sum(x > 0)))
colnames(promoter_DNase_peaks_sample) = "Count"
promoter_DNase_peaks_sample$Proportion = ifelse(metadata[match(rownames(promoter_DNase_peaks_sample),metadata$Sample),]$chrY == "Yes",promoter_DNase_peaks_sample$Count/NUM_PROMOTER,promoter_DNase_peaks_sample$Count/NUM_PROMOTER_noY)
promoter_DNase_peaks_sample$Sample = rownames(promoter_DNase_peaks_sample)
promoter_DNase_peaks_sample$State = rep("DNase",sample_counts["All","DNase"])

# Proportion of promoters overlapping H3K27ac peaks, by sample
promoter_H3K27ac_peaks_sample = as.data.frame(apply(promoter_H3K27ac_peaks[,5:102],2,function(x) sum(x > 0)))
colnames(promoter_H3K27ac_peaks_sample) = "Count"
promoter_H3K27ac_peaks_sample$Proportion = ifelse(metadata[match(rownames(promoter_H3K27ac_peaks_sample),metadata$Sample),]$chrY == "Yes",promoter_H3K27ac_peaks_sample$Count/NUM_PROMOTER,promoter_H3K27ac_peaks_sample$Count/NUM_PROMOTER_noY)
promoter_H3K27ac_peaks_sample$Sample = rownames(promoter_H3K27ac_peaks_sample)
promoter_H3K27ac_peaks_sample$State = rep("H3K27ac",sample_counts["All","H3K27ac"])