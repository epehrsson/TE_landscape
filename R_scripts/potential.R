# Creates dataframes of the number of TEs in each epigenetic state in each number of samples
# As well as the number of TEs in each state in at least one sample

## chromHMM_TE_state_dist: Number of TEs in each chromHMM state for each number of samples
## chromHMM_TE_state_dist_stats: Proportion of TEs ever in each chromHMM state and mean/SE proportion of samples in state
## state_sample_count: Number/proportion of TEs in each chromHMM state by sample
## TE_meth_average_category: Number of TEs in each methylation state for each number of samples
## combine_boxplot: TEs in each state per sample, all samples
## combine_potential: Number of TEs in each state in each number/proportion of samples, all samples
## combine_stats: Proportion of TEs ever in each state and mean/SE proportion of samples in state, all samples
## combine_boxplot_noCancer_IMR90: TEs in each state per sample, excluding cancer cell lines/IMR90
## combine_potential_noCancer: Number of TEs in each state in each number/proportion of samples, excluding cancer cell lines/IMR90
## combine_stats_noCancer_IMR90: Proportion of TEs ever in each state and mean/SE proportion of samples in state, excluding cancer cell lines/IMR90

# chromHMM potential

# All samples
## Number of TEs in each chromHMM state for each number of samples (0-127)
chromHMM_TE_state_dist = sample_distribution(chromHMM_TE_state,c(8:22),sample_counts["All","chromHMM"])

## Proportion of TEs ever in each state and mean/SE proportion of samples in state, for all TEs and those ever in the state
chromHMM_TE_state_dist_stats = potential_stats(chromHMM_TE_state_dist,15,sample_counts["All","chromHMM"])
chromHMM_TE_state_dist_stats$State = factor(chromHMM_states,levels=chromHMM_states)

# Excluding cancer cell lines/IMR90
## Number of TEs in each chromHMM state for each number of samples (0-121)
chromHMM_TE_state_dist_noCancer = sample_distribution(chromHMM_TE_state_noCancer,c(8:22),sample_counts["Include","chromHMM"])

## Proportion of TEs ever in each state and mean/SE proportion of samples in state, for all TEs and those ever in the state
chromHMM_TE_state_dist_noCancer_stats = potential_stats(chromHMM_TE_state_dist_noCancer,15,sample_counts["Include","chromHMM"])
chromHMM_TE_state_dist_noCancer_stats$State = factor(chromHMM_states,levels=chromHMM_states)

# Number/proportion of TEs in each chromHMM state by sample
state_sample_count = read.table("chromHMM/state_sample_counts_summit.txt",sep='\t',col.names=c("Sample","State","Count"))
state_sample_count[1905,] = c("E002","3_TxFlnk",0)
state_sample_count$Count = as.numeric(state_sample_count$Count)
state_sample_count$Proportion = ifelse(metadata[match(state_sample_count$Sample,metadata$Sample),]$chrY == "Yes",state_sample_count$Count/NUM_TE,state_sample_count$Count/NUM_TE_noY)
state_sample_count$State = factor(state_sample_count$State,levels=chromHMM_states)

# WGBS potential

# All samples
## Number of TEs in each methylation state for each number of samples (0-37)
TE_meth_average_category = sample_distribution(TE_meth_average,c(46:49),sample_counts["All","WGBS"])

## Proportion of TEs ever in each state and mean/SE proportion of samples in state, for all TEs and those ever in the state
TE_meth_average_category_stats = potential_stats(TE_meth_average_category,4,sample_counts["All","WGBS"])
TE_meth_average_category_stats$State = factor(rownames(TE_meth_average_category_stats),levels=meth_states)

# Excluding IMR90
## Number of TEs in each methylation state for each number of samples (0-36)
TE_meth_average_noIMR90_category = sample_distribution(TE_meth_average,c(50:53),sample_counts["Include","WGBS"])
colnames(TE_meth_average_noIMR90_category) = colnames(TE_meth_average_category) 

## Proportion of TEs ever in each state and mean/SE proportion of samples in state, for all TEs and those ever in the state
TE_meth_average_noIMR90_category_stats = potential_stats(TE_meth_average_noIMR90_category,4,sample_counts["Include","WGBS"])
TE_meth_average_noIMR90_category_stats$State = factor(rownames(TE_meth_average_noIMR90_category_stats),levels=meth_states)

# Number/proportion of TEs in each methylation state by sample
WGBS_sample_state = melt(TE_meth_average[,8:44],variable.name="Sample",value.name="Methylation")
WGBS_sample_state = ddply(WGBS_sample_state,.(Sample),summarise,
                                Hypomethylated=sum(na.omit(Methylation) < 0.3),Intermediate=sum(na.omit(Methylation) <= 0.7 & na.omit(Methylation) >= 0.3),
                                Hypermethylated=sum(na.omit(Methylation) > 0.7),Missing=sum(is.na(Methylation)))
WGBS_sample_state = melt(WGBS_sample_state,id.vars="Sample",variable.name="State",value.name="Count")
WGBS_sample_state$Proportion = WGBS_sample_state$Count/NUM_TE_WGBS

# DHS potential

# All samples
## Number of TEs overlapping a DHS peak summit for each number of samples (0-53)
TE_DNase_potential = sample_distribution(TE_DNase_peaks,62,sample_counts["All","DNase"])
TE_DNase_potential[1,2] = NUM_TE-dim(TE_DNase_peaks)[1]
colnames(TE_DNase_potential)[2] = "DNase"

## Proportion of TEs ever in each state and mean/SE proportion of samples in state, for all TEs and those ever in the state
TE_DNase_potential_stats = potential_stats(TE_DNase_potential,1,sample_counts["All","DNase"])
TE_DNase_potential_stats$State = "DNase"

# Excluding cancer cell lines/IMR90
## Number of TEs overlapping a DHS peak summit for each number of samples (0-48)
TE_DNase_potential_noCancer = sample_distribution(TE_DNase_peaks,63,sample_counts["Include","DNase"])
TE_DNase_potential_noCancer[1,2] = TE_DNase_potential_noCancer[1,2] + NUM_TE-dim(TE_DNase_peaks)[1]
colnames(TE_DNase_potential_noCancer)[2] = "DNase"

## Proportion of TEs ever in each state and mean/SE proportion of samples in state, for all TEs and those ever in the state
TE_DNase_potential_noCancer_stats = potential_stats(TE_DNase_potential_noCancer,1,sample_counts["Include","DNase"])
TE_DNase_potential_noCancer_stats$State = "DNase"

# Number/proportion of TEs overlapping a DHS peak summit by sample
TE_DNase_peaks_sample = as.data.frame(apply(TE_DNase_peaks[,8:60],2,function(x) sum(x > 0)))
colnames(TE_DNase_peaks_sample) = "Count"
TE_DNase_peaks_sample$Proportion = ifelse(metadata[match(rownames(TE_DNase_peaks_sample),metadata$Sample),]$chrY == "Yes",
                                          TE_DNase_peaks_sample$Count/NUM_TE,
                                          TE_DNase_peaks_sample$Count/NUM_TE_noY)
TE_DNase_peaks_sample$Sample = rownames(TE_DNase_peaks_sample)
TE_DNase_peaks_sample$State = rep("DNase",sample_counts["All","DNase"])

# H3K27ac potential

# All samples
## Number of TEs overlapping an H3K27ac peak summit for each number of samples (0-98)
TE_H3K27ac_potential = sample_distribution(TE_H3K27ac_peaks,107,sample_counts["All","H3K27ac"])
TE_H3K27ac_potential[1,2] = NUM_TE-dim(TE_H3K27ac_peaks)[1]
colnames(TE_H3K27ac_potential)[2] = "H3K27ac"

## Proportion of TEs ever in each state and mean/SE proportion of samples in state, for all TEs and those ever in the state
TE_H3K27ac_potential_stats = potential_stats(TE_H3K27ac_potential,1,sample_counts["All","H3K27ac"])
TE_H3K27ac_potential_stats$State = "H3K27ac"

# Excluding cancer cell lines/IMR90
## Number of TEs overlapping an H3K27ac peak summit for each number of samples (0-92)
TE_H3K27ac_potential_noCancer = sample_distribution(TE_H3K27ac_peaks,108,sample_counts["Include","H3K27ac"])
TE_H3K27ac_potential_noCancer[1,2] = TE_H3K27ac_potential_noCancer[1,2] + NUM_TE-dim(TE_H3K27ac_peaks)[1]
colnames(TE_H3K27ac_potential_noCancer)[2] = "H3K27ac"

## Proportion of TEs ever in each state and mean/SE proportion of samples in state, for all TEs and those ever in the state
TE_H3K27ac_potential_noCancer_stats = potential_stats(TE_H3K27ac_potential_noCancer,1,sample_counts["Include","H3K27ac"]) #Check
TE_H3K27ac_potential_noCancer_stats$State = "H3K27ac"

# Number/proportion of TEs overlapping an H3K27ac peak summit by sample
TE_H3K27ac_peaks_sample = as.data.frame(apply(TE_H3K27ac_peaks[,8:105],2,function(x) sum(x > 0)))
colnames(TE_H3K27ac_peaks_sample) = "Count"
TE_H3K27ac_peaks_sample$Proportion = ifelse(metadata[match(rownames(TE_H3K27ac_peaks_sample),metadata$Sample),]$chrY == "Yes",
                                            TE_H3K27ac_peaks_sample$Count/NUM_TE,
                                            TE_H3K27ac_peaks_sample$Count/NUM_TE_noY)
TE_H3K27ac_peaks_sample$Sample = rownames(TE_H3K27ac_peaks_sample)
TE_H3K27ac_peaks_sample$State = rep("H3K27ac",sample_counts["All","H3K27ac"])

# Expression potential

# All samples
## Number of TEs RPKM >1 for each number of samples (0-56)
RNA_potential = sample_distribution(RNA_TE,64,sample_counts["All","RNA"])
colnames(RNA_potential)[2] = "Expressed_samples"

## Proportion of TEs ever in each state and mean/SE proportion of samples in state, for all TEs and those ever in the state
RNA_potential_stats = potential_stats(RNA_potential,1,sample_counts["All","RNA"])
RNA_potential_stats$State = "Expressed_samples"

# Excluding cancer cell lines/IMR90
## Number of TEs RPKM >1 for each number of samples (0-52)
RNA_potential_noCancer = sample_distribution(RNA_TE,65,sample_counts["Include","RNA"])
colnames(RNA_potential_noCancer)[2] = "Expressed_samples"

## Proportion of TEs ever in each state and mean/SE proportion of samples in state, for all TEs and those ever in the state
RNA_potential_noCancer_stats = potential_stats(RNA_potential_noCancer,1,sample_counts["Include","RNA"])
RNA_potential_noCancer_stats$State = "Expressed_samples"

# Number/proportion of TEs RPKM >1 by sample
RNA_RPKM_sample = as.data.frame(apply(RNA_TE[,8:63],2,function(x) sum(x > 1)))
colnames(RNA_RPKM_sample) = "Count"
RNA_RPKM_sample$Proportion = ifelse(metadata[match(rownames(RNA_RPKM_sample),metadata$Sample),]$chrY == "Yes",
                                    RNA_RPKM_sample$Count/NUM_TE,
                                    RNA_RPKM_sample$Count/NUM_TE_noY)
RNA_RPKM_sample$Sample = rownames(RNA_RPKM_sample)
RNA_RPKM_sample$State = rep("Expressed_samples",sample_counts["All","RNA"])


# Combine dataframes for all techniques
## TEs in each state per sample, all samples
combine_boxplot = rbind(state_sample_count,WGBS_sample_state,TE_DNase_peaks_sample,TE_H3K27ac_peaks_sample,RNA_RPKM_sample)

## TEs in each state per sample, excluding cancer cell lines/IMR90
combine_boxplot_noCancer_IMR90 = droplevels(combine_boxplot[which(metadata[match(combine_boxplot$Sample,metadata$Sample),]$Exclude == "Include"),])
rm(list=c("WGBS_sample_state","TE_DNase_peaks_sample","TE_H3K27ac_peaks_sample","RNA_RPKM_sample"))

## Number of TEs in each state in each number/proportion of samples, all samples
combine_potential = rbind(melt(chromHMM_TE_state_dist,id.vars="Samples"),melt(TE_meth_average_category,id.vars="Samples"),melt(TE_DNase_potential,id.var="Samples"),
                          melt(TE_H3K27ac_potential,id.var="Samples"),melt(RNA_potential,id.var="Samples"))
colnames(combine_potential) = c("Samples","State","Count")
combine_potential = ddply(combine_potential,.(State),transform,Sample.Proportion = Samples/(length(Samples)-1))

## Number of TEs in each state in each number/proportion of samples, excluding cancer cell lines/IMR90
combine_potential_noCancer = rbind(melt(chromHMM_TE_state_dist_noCancer,id.vars="Samples"),melt(TE_meth_average_noIMR90_category,id.vars="Samples"),melt(TE_DNase_potential_noCancer,id.var="Samples"),
                          melt(TE_H3K27ac_potential_noCancer,id.var="Samples"),melt(RNA_potential_noCancer,id.var="Samples"))
colnames(combine_potential_noCancer) = c("Samples","State","Count")
combine_potential_noCancer = ddply(combine_potential_noCancer,.(State),transform,Sample.Proportion = Samples/(length(Samples)-1))

rm(list=c("chromHMM_TE_state_dist_noCancer","TE_meth_average_noIMR90_category",
          "TE_DNase_potential","TE_DNase_potential_noCancer","TE_H3K27ac_potential","TE_H3K27ac_potential_noCancer","RNA_potential","RNA_potential_noCancer"))

## Proportion of TEs ever in each state and mean/SE proportion of samples in state, all samples
combine_stats = rbind(chromHMM_TE_state_dist_stats,TE_meth_average_category_stats,TE_DNase_potential_stats,TE_H3K27ac_potential_stats,RNA_potential_stats)
combine_stats$Group = factor(c(rep("chromHMM",15),rep("WGBS",4),"DNase","H3K27ac","Expressed_samples"),
                             levels=c("chromHMM","WGBS","DNase","H3K27ac","Expressed_samples"))

## Proportion of TEs ever in each state and mean/SE proportion of samples in state, excluding cancer cell lines/IMR90 
combine_stats_noCancer_IMR90 = rbind(chromHMM_TE_state_dist_noCancer_stats,TE_meth_average_noIMR90_category_stats,
                                     TE_DNase_potential_noCancer_stats,TE_H3K27ac_potential_noCancer_stats,RNA_potential_noCancer_stats)
combine_stats_noCancer_IMR90$Group = factor(c(rep("chromHMM",15),rep("WGBS",4),"DNase","H3K27ac","Expressed_samples"),
                                            levels=c("chromHMM","WGBS","DNase","H3K27ac","Expressed_samples"))

rm(list=c("chromHMM_TE_state_dist_noCancer_stats","TE_meth_average_category_stats","TE_meth_average_noIMR90_category_stats",
          "TE_DNase_potential_stats","TE_DNase_potential_noCancer_stats","TE_H3K27ac_potential_stats","TE_H3K27ac_potential_noCancer_stats",
          "RNA_potential_stats","RNA_potential_noCancer_stats"))