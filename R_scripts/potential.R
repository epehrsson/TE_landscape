# Potential 
# See 4/25/2016, 4/26/2016, 4/27/2016, 5/3/2016, 5/11/2016, 5/12/2016, 6/2/2016, 6/27/2016, 7/11/2016, 8/25/2016, 8/26/2016, 8/30/2016, 9/16/2016, 9/20/2016, 9/28/2016, 12/15/2016, 2/3/2017, 2/6/2017, 5/10/2017, 6/14/2017, 7/22/2017, 8/1/2017, 8/2/2017, 8/3/2017

# chromHMM potential
#load("R_datasets/chromHMM_TE_state.RData")

# All samples
chromHMM_TE_state_dist = sample_distribution(chromHMM_TE_state,c(8:22),127)
chromHMM_TE_state_cum = cumulative_distribution(chromHMM_TE_state,c(8:22),127)
chromHMM_TE_state_cum$State = factor(rep(chromHMM_states,each=127),levels=chromHMM_states)
chromHMM_TE_state_dist_stats = potential_stats(chromHMM_TE_state_dist,15,127)
chromHMM_TE_state_dist_stats$State = factor(chromHMM_states,levels=chromHMM_states)

# No cancer cell lines, IMR90
chromHMM_TE_state_dist_noCancer = sample_distribution(chromHMM_TE_state_noCancer,c(8:22),121)
chromHMM_TE_state_noCancer_cum = cumulative_distribution(chromHMM_TE_state_noCancer,c(8:22),121)
chromHMM_TE_state_noCancer_cum$State = factor(rep(chromHMM_states,each=121),levels=chromHMM_states)
chromHMM_TE_state_dist_noCancer_stats = potential_stats(chromHMM_TE_state_dist_noCancer,15,121)
chromHMM_TE_state_dist_noCancer_stats$State = factor(chromHMM_states,levels=chromHMM_states)

# Number of TEs in each chromHMM state by sample
state_sample_count = read.table("chromHMM/subfamily/state_sample_counts.txt",sep='\t')
colnames(state_sample_count) = c("State","Sample","Count")
state_sample_count[1905,] = c("3_TxFlnk","E002",0)
state_sample_count$Count = as.numeric(state_sample_count$Count)
state_sample_count$Proportion = state_sample_count$Count/4430788
state_sample_count[which(metadata[match(state_sample_count$Sample,metadata$Sample),]$chrY == "No"),]$Proportion = state_sample_count[which(metadata[match(state_sample_count$Sample,metadata$Sample),]$chrY == "No"),]$Count/(4430788-31580)
state_sample_count$State = factor(state_sample_count$State,levels=chromHMM_states)

# WGBS potential
#load("R_datasets/TE_meth_average.RData")

# Cumulative distribution of methylation states and statistics
TE_meth_average_category = sample_distribution(TE_meth_average,c(46:49),37)
TE_meth_average_category_cum = cumulative_distribution(TE_meth_average,c(46:49),37)
TE_meth_average_category_stats = potential_stats(TE_meth_average_category,4,37)
TE_meth_average_category_stats$State = factor(rownames(TE_meth_average_category_stats),levels=meth_states)

# Cumulative distribution of methylation states and statistics, no IMR90
TE_meth_average_noIMR90_category = sample_distribution(TE_meth_average,c(50:53),36)
TE_meth_average_noIMR90_category_cum = cumulative_distribution(TE_meth_average,c(50:53),36)
TE_meth_average_noIMR90_category_stats = potential_stats(TE_meth_average_noIMR90_category,4,36)
TE_meth_average_noIMR90_category_stats$State = factor(c("Hypomethylated","Hypermethylated","Intermediate","Missing"),levels=meth_states)

# Number of TEs in each WGBS state by sample
TE_meth_average_state = as.data.frame(cbind(apply(TE_meth_average[,8:44],2,function(x) sum(na.omit(x) < 0.3)/length(x)),apply(TE_meth_average[,8:44],2,function(x)  sum(na.omit(x) <= 0.7 & na.omit(x) >= 0.3)/length(x)),apply(TE_meth_average[,8:44],2,function(x) sum(na.omit(x) > 0.7)/length(x)),apply(TE_meth_average[,8:44],2,function(x) sum(is.na(x))/length(x))))
colnames(TE_meth_average_state) = c("Hypomethylated","Intermediate","Hypermethylated","Missing")
TE_meth_average_state = TE_meth_average_state[order(TE_meth_average_state$Hypermethylated + TE_meth_average_state$Missing),]

TE_meth_average_state_long = melt(as.matrix(TE_meth_average_state))
colnames(TE_meth_average_state_long) = c("Sample","State","Proportion")

# DNase potential
#load("R_datasets/TE_DNase_peaks.RData")

# Distribution of TEs overlapping DNase peaks
TE_DNase_potential = sample_distribution(TE_DNase_peaks,62,53)
TE_DNase_potential[1,2] = 4430788-1840984
colnames(TE_DNase_potential)[2] = "DNase"
TE_DNase_potential_cum = cumulative_distribution(TE_DNase_peaks,62,53)
TE_DNase_potential_cum$State = rep("DNase",53)
TE_DNase_potential_stats = potential_stats(TE_DNase_potential,1,53)
TE_DNase_potential_stats$State = "DNase"

# No cancer cell lines/IMR90
TE_DNase_potential_noCancer = sample_distribution(TE_DNase_peaks,63,48)
TE_DNase_potential_noCancer[1,2] = 4430788-1840984
colnames(TE_DNase_potential_noCancer)[2] = "DNase"
TE_DNase_potential_noCancer_cum = cumulative_distribution(TE_DNase_peaks,63,48)
TE_DNase_potential_noCancer_cum$State = rep("DNase",48)
TE_DNase_potential_noCancer_stats = potential_stats(TE_DNase_potential_noCancer,1,48) #Check
TE_DNase_potential_noCancer_stats$State = "DNase"

# Proportion of TEs overlapping Dnase peaks, by sample
TE_DNase_peaks_sample = as.data.frame(apply(TE_DNase_peaks[,8:60],2,function(x) sum(x > 0)))
colnames(TE_DNase_peaks_sample) = "Count"
TE_DNase_peaks_sample$Proportion = TE_DNase_peaks_sample$Count/4430788
TE_DNase_peaks_sample[which(metadata[match(TE_DNase_peaks_sample$Sample,metadata$Sample),]$chrY == "No"),]$Proportion = TE_DNase_peaks_sample[which(metadata[match(TE_DNase_peaks_sample$Sample,metadata$Sample),]$chrY == "No"),]$Count/4399208
TE_DNase_peaks_sample$Sample = rownames(TE_DNase_peaks_sample)
TE_DNase_peaks_sample$State = rep("DNase",53)

# H3K27ac potential
#load("R_datasets/TE_H3K27ac_peaks.RData")

# Distribution of TEs overlapping H3K27ac peaks
TE_H3K27ac_potential = sample_distribution(TE_H3K27ac_peaks,107,98)
TE_H3K27ac_potential[1,2] = 4430788-1801527
colnames(TE_H3K27ac_potential)[2] = "H3K27ac"
TE_H3K27ac_potential_cum = cumulative_distribution(TE_H3K27ac_peaks,107,98)
TE_H3K27ac_potential_cum$State = rep("H3K27ac",98)
TE_H3K27ac_potential_stats = potential_stats(TE_H3K27ac_potential,1,98)
TE_H3K27ac_potential_stats$State = "H3K27ac"

# No cancer cell lines/IMR90
TE_H3K27ac_potential_noCancer = sample_distribution(TE_H3K27ac_peaks,108,92)
TE_H3K27ac_potential_noCancer[1,2] = 4430788-1801527
colnames(TE_H3K27ac_potential_noCancer)[2] = "H3K27ac"
TE_H3K27ac_potential_noCancer_cum = cumulative_distribution(TE_H3K27ac_peaks,108,92)
TE_H3K27ac_potential_noCancer_cum$State = rep("H3K27ac",92)
TE_H3K27ac_potential_noCancer_stats = potential_stats(TE_H3K27ac_potential_noCancer,1,92) #Check
TE_H3K27ac_potential_noCancer_stats$State = "H3K27ac"

# Proportion of TEs overlapping H3K27ac peaks, by sample
TE_H3K27ac_peaks_sample = as.data.frame(apply(TE_H3K27ac_peaks[,8:105],2,function(x) sum(x > 0)))
colnames(TE_H3K27ac_peaks_sample) = "Count"
TE_H3K27ac_peaks_sample$Proportion = TE_H3K27ac_peaks_sample$Count/4430788
TE_H3K27ac_peaks_sample[which(metadata[match(TE_H3K27ac_peaks_sample$Sample,metadata$Sample),]$chrY == "No"),]$Proportion = TE_H3K27ac_peaks_sample[which(metadata[match(TE_H3K27ac_peaks_sample$Sample,metadata$Sample),]$chrY == "No"),]$Count/4399208
TE_H3K27ac_peaks_sample$Sample = rownames(TE_H3K27ac_peaks_sample)
TE_H3K27ac_peaks_sample$State = rep("H3K27ac",98)

# RNA-seq potential
#load("R_datasets/rna.RData")

# Distribution of TEs with RPKM >1
RNA_potential = sample_distribution(RNA_TE_agnostic,61,52)
colnames(RNA_potential)[2] = "Expression"
RNA_potential_cum = cumulative_distribution(RNA_TE_agnostic,61,52)
RNA_potential_cum$State = rep("Expression",52)
RNA_potential_stats = potential_stats(RNA_potential,1,52)
RNA_potential_stats$State = "Expression"

# Distribution of TEs with RPKM >1, no cancer cell lines/IMR90
RNA_potential_noCancer = sample_distribution(RNA_TE_agnostic,62,48)
colnames(RNA_potential_noCancer)[2] = "Expression"
RNA_potential_noCancer_cum = cumulative_distribution(RNA_TE_agnostic,62,48)
RNA_potential_noCancer_cum$State = rep("Expression",48)
RNA_potential_noCancer_stats = potential_stats(RNA_potential_noCancer,1,48)
RNA_potential_noCancer_stats$State = "Expression"

# Proportion of TEs with RPKM >1, by sample
RNA_RPKM_sample = as.data.frame(apply(RNA_TE_agnostic[,9:60],2,function(x) sum(x > 1)))
colnames(RNA_RPKM_sample) = "Count"
RNA_RPKM_sample$Proportion = RNA_RPKM_sample$Count/4430788
RNA_RPKM_sample[which(metadata[match(RNA_RPKM_sample$Sample,metadata$Sample),]$chrY == "No"),]$Proportion = RNA_RPKM_sample[which(metadata[match(RNA_RPKM_sample$Sample,metadata$Sample),]$chrY == "No"),]$Count/4399208
RNA_RPKM_sample$Sample = rownames(RNA_RPKM_sample)
RNA_RPKM_sample$State = rep("Expression",52)
