# Potential 
# See 4/25/2016, 4/26/2016, 4/27/2016, 5/3/2016, 5/11/2016, 5/12/2016, 6/2/2016, 6/27/2016, 7/11/2016, 8/25/2016, 8/26/2016, 8/30/2016, 9/16/2016, 9/20/2016, 9/28/2016, 12/15/2016, 2/3/2017, 2/6/2017, 5/10/2017, 6/14/2017, 7/22/2017, 8/1/2017, 8/2/2017, 8/3/2017

# chromHMM potential

# All samples
chromHMM_TE_state_dist = sample_distribution(chromHMM_TE_state,c(8:22),sample_counts["All","chromHMM"])
chromHMM_TE_state_cum = cumulative_distribution(chromHMM_TE_state,c(8:22),sample_counts["All","chromHMM"])
chromHMM_TE_state_cum$State = factor(rep(chromHMM_states,each=sample_counts["All","chromHMM"]),levels=chromHMM_states)
chromHMM_TE_state_dist_stats = potential_stats(chromHMM_TE_state_dist,15,sample_counts["All","chromHMM"])
chromHMM_TE_state_dist_stats$State = factor(chromHMM_states,levels=chromHMM_states)

# No cancer cell lines, IMR90
chromHMM_TE_state_dist_noCancer = sample_distribution(chromHMM_TE_state_noCancer,c(8:22),sample_counts["Include","chromHMM"])
chromHMM_TE_state_dist_noCancer_stats = potential_stats(chromHMM_TE_state_dist_noCancer,15,sample_counts["Include","chromHMM"])
chromHMM_TE_state_dist_noCancer_stats$State = factor(chromHMM_states,levels=chromHMM_states)

# Number of TEs in each chromHMM state by sample
state_sample_count = read.table("chromHMM/state_sample_counts_summit.txt",sep='\t')
colnames(state_sample_count) = c("Sample","State","Count")
state_sample_count[1905,] = c("E002","3_TxFlnk",0)
state_sample_count$Count = as.numeric(state_sample_count$Count)
state_sample_count$Proportion = ifelse(metadata[match(state_sample_count$Sample,metadata$Sample),]$chrY == "Yes",state_sample_count$Count/NUM_TE,state_sample_count$Count/NUM_TE_noY)
state_sample_count$State = factor(state_sample_count$State,levels=chromHMM_states)

# WGBS potential

# Cumulative distribution of methylation states and statistics
TE_meth_average_category = sample_distribution(TE_meth_average,c(46:49),sample_counts["All","WGBS"])
TE_meth_average_category_cum = cumulative_distribution(TE_meth_average,c(46:49),sample_counts["All","WGBS"])
TE_meth_average_category_stats = potential_stats(TE_meth_average_category,4,sample_counts["All","WGBS"])
TE_meth_average_category_stats$State = factor(rownames(TE_meth_average_category_stats),levels=meth_states)

# Cumulative distribution of methylation states and statistics, no IMR90
TE_meth_average_noIMR90_category = sample_distribution(TE_meth_average,c(50:53),sample_counts["Include","WGBS"])
colnames(TE_meth_average_noIMR90_category) = colnames(TE_meth_average_category) 
TE_meth_average_noIMR90_category_stats = potential_stats(TE_meth_average_noIMR90_category,4,sample_counts["Include","WGBS"])
TE_meth_average_noIMR90_category_stats$State = factor(rownames(TE_meth_average_noIMR90_category_stats),levels=meth_states)

# Number of TEs in each WGBS state by sample
WGBS_sample_state = melt(TE_meth_average[,8:44])
colnames(WGBS_sample_state) = c("Sample","Methylation")
WGBS_sample_state = ddply(WGBS_sample_state,.(Sample),summarise,
                                Hypomethylated=sum(na.omit(Methylation) < 0.3),Intermediate=sum(na.omit(Methylation) <= 0.7 & na.omit(Methylation) >= 0.3),
                                Hypermethylated=sum(na.omit(Methylation) > 0.7),Missing=sum(is.na(Methylation)))
WGBS_sample_state = melt(WGBS_sample_state,id.vars="Sample")
colnames(WGBS_sample_state) = c("Sample","State","Count")
WGBS_sample_state$Proportion = WGBS_sample_state$Count/NUM_TE_WGBS

# DNase potential

# Distribution of TEs overlapping DNase peaks
TE_DNase_potential = sample_distribution(TE_DNase_peaks,62,sample_counts["All","DNase"])
TE_DNase_potential[1,2] = NUM_TE-dim(TE_DNase_peaks)[1]
colnames(TE_DNase_potential)[2] = "DNase"
TE_DNase_potential_cum = cumulative_distribution(TE_DNase_peaks,62,sample_counts["All","DNase"])
TE_DNase_potential_cum$State = rep("DNase",sample_counts["All","DNase"])
TE_DNase_potential_stats = potential_stats(TE_DNase_potential,1,sample_counts["All","DNase"])
TE_DNase_potential_stats$State = "DNase"

# No cancer cell lines/IMR90
TE_DNase_potential_noCancer = sample_distribution(TE_DNase_peaks,63,sample_counts["Include","DNase"])
TE_DNase_potential_noCancer[1,2] = TE_DNase_potential_noCancer[1,2] + NUM_TE-dim(TE_DNase_peaks)[1]
colnames(TE_DNase_potential_noCancer)[2] = "DNase"
TE_DNase_potential_noCancer_stats = potential_stats(TE_DNase_potential_noCancer,1,sample_counts["Include","DNase"])
TE_DNase_potential_noCancer_stats$State = "DNase"

# Proportion of TEs overlapping DNase peaks, by sample
TE_DNase_peaks_sample = as.data.frame(apply(TE_DNase_peaks[,8:60],2,function(x) sum(x > 0)))
colnames(TE_DNase_peaks_sample) = "Count"
TE_DNase_peaks_sample$Proportion = ifelse(metadata[match(rownames(TE_DNase_peaks_sample),metadata$Sample),]$chrY == "Yes",
                                          TE_DNase_peaks_sample$Count/NUM_TE,
                                          TE_DNase_peaks_sample$Count/NUM_TE_noY)
TE_DNase_peaks_sample$Sample = rownames(TE_DNase_peaks_sample)
TE_DNase_peaks_sample$State = rep("DNase",sample_counts["All","DNase"])

# H3K27ac potential

# Distribution of TEs overlapping H3K27ac peaks
TE_H3K27ac_potential = sample_distribution(TE_H3K27ac_peaks,107,sample_counts["All","H3K27ac"])
TE_H3K27ac_potential[1,2] = NUM_TE-dim(TE_H3K27ac_peaks)[1]
colnames(TE_H3K27ac_potential)[2] = "H3K27ac"
TE_H3K27ac_potential_cum = cumulative_distribution(TE_H3K27ac_peaks,107,sample_counts["All","H3K27ac"])
TE_H3K27ac_potential_cum$State = rep("H3K27ac",sample_counts["All","H3K27ac"])
TE_H3K27ac_potential_stats = potential_stats(TE_H3K27ac_potential,1,sample_counts["All","H3K27ac"])
TE_H3K27ac_potential_stats$State = "H3K27ac"

# No cancer cell lines/IMR90
TE_H3K27ac_potential_noCancer = sample_distribution(TE_H3K27ac_peaks,108,sample_counts["Include","H3K27ac"])
TE_H3K27ac_potential_noCancer[1,2] = TE_H3K27ac_potential_noCancer[1,2] + NUM_TE-dim(TE_H3K27ac_peaks)[1]
colnames(TE_H3K27ac_potential_noCancer)[2] = "H3K27ac"
TE_H3K27ac_potential_noCancer_stats = potential_stats(TE_H3K27ac_potential_noCancer,1,sample_counts["Include","H3K27ac"]) #Check
TE_H3K27ac_potential_noCancer_stats$State = "H3K27ac"

# Proportion of TEs overlapping H3K27ac peaks, by sample
TE_H3K27ac_peaks_sample = as.data.frame(apply(TE_H3K27ac_peaks[,8:105],2,function(x) sum(x > 0)))
colnames(TE_H3K27ac_peaks_sample) = "Count"
TE_H3K27ac_peaks_sample$Proportion = ifelse(metadata[match(rownames(TE_H3K27ac_peaks_sample),metadata$Sample),]$chrY == "Yes",
                                            TE_H3K27ac_peaks_sample$Count/NUM_TE,
                                            TE_H3K27ac_peaks_sample$Count/NUM_TE_noY)
TE_H3K27ac_peaks_sample$Sample = rownames(TE_H3K27ac_peaks_sample)
TE_H3K27ac_peaks_sample$State = rep("H3K27ac",sample_counts["All","H3K27ac"])

# RNA-seq potential

# Distribution of TEs with RPKM >1
RNA_potential = sample_distribution(RNA_TE,64,sample_counts["All","RNA"])
colnames(RNA_potential)[2] = "Expressed_samples"
RNA_potential_cum = cumulative_distribution(RNA_TE,64,sample_counts["All","RNA"])
RNA_potential_cum$State = rep("Expressed_samples",sample_counts["All","RNA"])
RNA_potential_stats = potential_stats(RNA_potential,1,sample_counts["All","RNA"])
RNA_potential_stats$State = "Expressed_samples"

# Distribution of TEs with RPKM >1, no cancer cell lines/IMR90
RNA_potential_noCancer = sample_distribution(RNA_TE,65,sample_counts["Include","RNA"])
colnames(RNA_potential_noCancer)[2] = "Expressed_samples"
RNA_potential_noCancer_stats = potential_stats(RNA_potential_noCancer,1,sample_counts["Include","RNA"])
RNA_potential_noCancer_stats$State = "Expressed_samples"

# Proportion of TEs with RPKM >1, by sample
RNA_RPKM_sample = as.data.frame(apply(RNA_TE[,8:63],2,function(x) sum(x > 1)))
colnames(RNA_RPKM_sample) = "Count"
RNA_RPKM_sample$Proportion = ifelse(metadata[match(rownames(RNA_RPKM_sample),metadata$Sample),]$chrY == "Yes",
                                    RNA_RPKM_sample$Count/NUM_TE,
                                    RNA_RPKM_sample$Count/NUM_TE_noY)
RNA_RPKM_sample$Sample = rownames(RNA_RPKM_sample)
RNA_RPKM_sample$State = rep("Expressed_samples",sample_counts["All","RNA"])

# Combine marks
## TEs in state per sample
combine_boxplot = rbind(state_sample_count,WGBS_sample_state,TE_DNase_peaks_sample,TE_H3K27ac_peaks_sample,RNA_RPKM_sample)
combine_boxplot_noCancer_IMR90 = droplevels(combine_boxplot[which(metadata[match(combine_boxplot$Sample,metadata$Sample),]$Exclude == "Include"),])
rm(list=c("WGBS_sample_state","TE_DNase_peaks_sample","TE_H3K27ac_peaks_sample","RNA_RPKM_sample"))

## Combined samples in state
combine_potential = rbind(melt(chromHMM_TE_state_dist,id.vars="Samples"),melt(TE_meth_average_category,id.vars="Samples"),melt(TE_DNase_potential,id.var="Samples"),
                          melt(TE_H3K27ac_potential,id.var="Samples"),melt(RNA_potential,id.var="Samples"))
colnames(combine_potential) = c("Samples","State","Count")
combine_potential = ddply(combine_potential,.(State),transform,Sample.Proportion = Samples/(length(Samples)-1))

## Combined samples in state, no cancer/IMR90
combine_potential_noCancer = rbind(melt(chromHMM_TE_state_dist_noCancer,id.vars="Samples"),melt(TE_meth_average_noIMR90_category,id.vars="Samples"),melt(TE_DNase_potential_noCancer,id.var="Samples"),
                          melt(TE_H3K27ac_potential_noCancer,id.var="Samples"),melt(RNA_potential_noCancer,id.var="Samples"))
colnames(combine_potential_noCancer) = c("Samples","State","Count")
combine_potential_noCancer = ddply(combine_potential_noCancer,.(State),transform,Sample.Proportion = Samples/(length(Samples)-1))

rm(list=c("chromHMM_TE_state_dist_noCancer","TE_meth_average_noIMR90_category",
          "TE_DNase_potential","TE_DNase_potential_noCancer","TE_H3K27ac_potential","TE_H3K27ac_potential_noCancer","RNA_potential","RNA_potential_noCancer"))

## Stats
combine_stats = rbind(chromHMM_TE_state_dist_stats,TE_meth_average_category_stats,TE_DNase_potential_stats,TE_H3K27ac_potential_stats,RNA_potential_stats)
combine_stats$Group = factor(c(rep("chromHMM",15),rep("WGBS",4),"DNase","H3K27ac","Expressed_samples"),
                             levels=c("chromHMM","WGBS","DNase","H3K27ac","Expressed_samples"))

combine_stats_noCancer_IMR90 = rbind(chromHMM_TE_state_dist_noCancer_stats,TE_meth_average_noIMR90_category_stats,
                                     TE_DNase_potential_noCancer_stats,TE_H3K27ac_potential_noCancer_stats,RNA_potential_noCancer_stats)
combine_stats_noCancer_IMR90$Group = factor(c(rep("chromHMM",15),rep("WGBS",4),"DNase","H3K27ac","Expressed_samples"),
                                            levels=c("chromHMM","WGBS","DNase","H3K27ac","Expressed_samples"))

rm(list=c("chromHMM_TE_state_dist_noCancer_stats","TE_meth_average_category_stats","TE_meth_average_noIMR90_category_stats",
          "TE_DNase_potential_stats","TE_DNase_potential_noCancer_stats","TE_H3K27ac_potential_stats","TE_H3K27ac_potential_noCancer_stats",
          "RNA_potential_stats","RNA_potential_noCancer_stats"))

## Cumulative distribution
combine_cum = rbind(TE_meth_average_category_cum,TE_DNase_potential_cum,TE_H3K27ac_potential_cum,RNA_potential_cum)
rm(list=c("TE_meth_average_category_cum","TE_DNase_potential_cum","TE_H3K27ac_potential_cum","RNA_potential_cum"))