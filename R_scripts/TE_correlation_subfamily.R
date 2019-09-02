# Creates a single dataframe with TE subfamily characteristics and number of enrichments in each epigenetic state, 
# then finds the Spearman correlation between subfamily characteristics and number of enrichments in each state

## rmsk_TE_subfamily_measure - For each TE subfamily, the proportion of samples enriched LOR > 1.5 in each state
## or representing >1% of the state, plus subfamily statistics

## correlate_subfamily - Spearman correlation between the proportion of samples a subfamily is enriched in each state and subfamily characteristics
## correlate_subfam_feature - Spearman correlation between subfamily characteristics
## correlate_subfam_state - Spearman correlation between the proportion of samples a subfamily is enriched in each state

# To the dataframe of TE subfamily statistics, adds the proportion of samples the subfamily is enriched LOR > 1.5 in each state
# including only those that pass subfamily member thresholds
rmsk_TE_subfamily_measure = merge(rmsk_TE_subfamily,
                                  dcast(subfamily_state_sample_counts,subfamily~State,
                                        value.var=c("Sample.Proportion"))[,c("subfamily",chromHMM_states,meth_states,"DNase","H3K27ac")],
                                  by="subfamily",all.x=TRUE)

# Adds the proportion of samples the subfamily represents >1% of each state
test = dcast(subfamily_state_sample_counts_pc,subfamily~State,
             value.var=c("Sample.Proportion"))[,c("subfamily",chromHMM_states,meth_states,"DNase","H3K27ac")]
colnames(test)[2:22] = lapply(colnames(test)[2:22],function(x) paste(x,"PC",sep="_"))
rmsk_TE_subfamily_measure = merge(rmsk_TE_subfamily_measure,test,by="subfamily",all.x=TRUE)
rm(test)

# Set all class levels as equivalent (none is baseline)
contrasts(rmsk_TE_subfamily_measure$class_update) <- contr.sum


all_metrics = c(cohorts,standard_chromosomes[1:24],measure_metrics_subfam)

# Finds the Spearman correlation across subfamilies between: 
# A) the proportion of samples each subfamily is enriched in each state or represents >1% of the state and
# B) the proportion of subfamily members overlapping each genic feature and on each chromosome, 
# and subfamily features (number of TEs, total length, median TE length, median mappability, median JC distance, 
# total number of CpGs, number of TEs with CpGs, mean CpGs per TE, and mean CpGs per kbp of subfamily),
# for all TEs and for each of the six TE classes
# FDR correction for p-values across all comparisons

correlate_subfamily = ldply(all_metrics,function(x) correlate_spearman(rmsk_TE_subfamily_measure,x,all_pc_states))
correlate_subfamily$Metric = rep(all_metrics,each=7*length(all_pc_states))
correlate_subfamily$Metric = factor(correlate_subfamily$Metric,levels=all_metrics)
correlate_subfamily$p.value = as.numeric(as.character(correlate_subfamily$p.value))
correlate_subfamily$p.adjust = p.adjust(correlate_subfamily$p.value,method="fdr")
correlate_subfamily$estimate.rho = as.numeric(as.character(correlate_subfamily$estimate.rho))
correlate_subfamily$class_update = factor(correlate_subfamily$class_update,levels=c("DNA","LINE","LTR","SINE","SVA","Other","All"))
correlate_subfamily$Category = ifelse(correlate_subfamily$Metric %in% standard_chromosomes,"Chromosome",
                                      ifelse(correlate_subfamily$Metric %in% cohorts,"Feature","Sequence"))

# Finds the Spearman correlation across subfamilies between: 
# the proportion of subfamily members overlapping each genic feature and on each chromosome, 
# and subfamily features (number of TEs, total length, median TE length, median mappability, median JC distance, 
# total number of CpGs, number of TEs with CpGs, mean CpGs per TE, and mean CpGs per kbp of subfamily),
# for all TEs and for each of the six TE classes
# FDR correction for p-values across all comparisons

correlate_subfam_feature = ldply(all_metrics,function(x) correlate_spearman(rmsk_TE_subfamily_measure,x,all_metrics))
correlate_subfam_feature$Metric = rep(all_metrics,each=7*length(all_metrics))
correlate_subfam_feature$p.value = as.numeric(as.character(correlate_subfam_feature$p.value))
correlate_subfam_feature$p.adjust = p.adjust(correlate_subfam_feature$p.value,method="fdr")
correlate_subfam_feature$estimate.rho = as.numeric(as.character(correlate_subfam_feature$estimate.rho))
correlate_subfam_feature$class_update = factor(correlate_subfam_feature$class_update,levels=c("DNA","LINE","LTR","SINE","SVA","Other","All"))

# Finds the Spearman correlation across subfamilies between: 
# the proportion of samples each subfamily is enriched in each state or represents >1% of the state,
# for all TEs and for each of the six TE classes
# FDR correction for p-values across all comparisons

correlate_subfam_state = ldply(all_pc_states,function(x) correlate_spearman(rmsk_TE_subfamily_measure,x,all_pc_states))
correlate_subfam_state$State2 = rep(all_pc_states,each=7*length(all_pc_states))
correlate_subfam_state$p.value = as.numeric(as.character(correlate_subfam_state$p.value))
correlate_subfam_state$p.adjust = p.adjust(correlate_subfam_state$p.value,method="fdr")
correlate_subfam_state$estimate.rho = as.numeric(as.character(correlate_subfam_state$estimate.rho))
correlate_subfam_state$class_update = factor(correlate_subfam_state$class_update,levels=c("DNA","LINE","LTR","SINE","SVA","Other","All"))