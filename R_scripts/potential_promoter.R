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
promoter_state_sample_count = read.table("chromHMM/promoters_state_sample_counts_summit.txt",sep='\t')
colnames(promoter_state_sample_count) = c("Sample","State","Count")
promoter_state_sample_count[1905,] = c("E002","3_TxFlnk",0)
promoter_state_sample_count$Count = as.numeric(promoter_state_sample_count$Count)
promoter_state_sample_count$Proportion = ifelse(metadata[match(promoter_state_sample_count$Sample,metadata$Sample),]$chrY == "Yes",promoter_state_sample_count$Count/NUM_PROMOTER,promoter_state_sample_count$Count/NUM_PROMOTER_noY)
promoter_state_sample_count$State = factor(promoter_state_sample_count$State,levels=chromHMM_states)

# Number of promoters in each WGBS state by sample
promoter_meth_average_state = melt(promoter_meth_average[,5:41])
colnames(promoter_meth_average_state) = c("Sample","Methylation")
promoter_meth_average_state = ddply(promoter_meth_average_state,.(Sample),summarise,
                                    Hypomethylated=sum(na.omit(Methylation) < 0.3),Intermediate=sum(na.omit(Methylation) <= 0.7 & na.omit(Methylation) >= 0.3),
                                    Hypermethylated=sum(na.omit(Methylation) > 0.7),Missing=sum(is.na(Methylation)))
promoter_meth_average_state = melt(promoter_meth_average_state,id.vars="Sample")
colnames(promoter_meth_average_state) = c("Sample","State","Count")
promoter_meth_average_state$Proportion = promoter_meth_average_state$Count/NUM_PROMOTER_noY

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

# Combine
combine_boxplot_prom = rbind(promoter_state_sample_count,promoter_meth_average_state,promoter_DNase_peaks_sample,promoter_H3K27ac_peaks_sample)
combine_stats_prom = rbind(chromHMM_promoter_state_dist_stats,promoter_meth_average_category_stats,promoter_DNase_potential_stats,promoter_H3K27ac_potential_stats)
combine_stats_prom$Group = factor(c(rep("chromHMM",15),rep("WGBS",4),"DNase","H3K27ac"),levels=c("chromHMM","WGBS","DNase","H3K27ac"))
combine_cum_prom = rbind(chromHMM_promoter_state_cum,promoter_meth_average_category_cum,promoter_DNase_potential_cum,promoter_H3K27ac_potential_cum)
combine_potential_prom = rbind(melt(chromHMM_promoter_state_dist,id.vars="Samples"),melt(promoter_meth_average_category,id.vars="Samples"),
                               melt(promoter_DNase_potential,id.var="Samples"),melt(promoter_H3K27ac_potential,id.var="Samples"))
colnames(combine_potential_prom) = c("Samples","State","Count")
combine_potential_prom = ddply(combine_potential_prom,.(State),transform,Sample.Proportion = Samples/(length(Samples)-1))

rm(list=c("chromHMM_promoter_state_dist","promoter_meth_average_category","promoter_DNase_potential","promoter_H3K27ac_potential",
          "chromHMM_promoter_state","promoter_meth_average","promoter_DNase_peaks","promoter_H3K27ac_peaks",
          "promoter_state_sample_count","promoter_meth_average_state","promoter_DNase_peaks_sample","promoter_H3K27ac_peaks_sample",
          "chromHMM_promoter_state_dist_stats","promoter_meth_average_category_stats","promoter_DNase_potential_stats","promoter_H3K27ac_potential_stats",
          "chromHMM_promoter_state_cum","promoter_meth_average_category_cum","promoter_DNase_potential_cum","promoter_H3K27ac_potential_cum"))
