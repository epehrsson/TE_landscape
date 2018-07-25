all_metrics = c(cohorts,standard_chromosomes[1:24],measure_metrics_subfam)

# Correlate samples enriched/>1% in state with subfamily features
correlate_subfamily = ldply(all_metrics,function(x) correlate_spearman(rmsk_TE_subfamily_measure,x,all_pc_states))
correlate_subfamily$Metric = rep(all_metrics,each=7*length(all_pc_states))
correlate_subfamily$Metric = factor(correlate_subfamily$Metric,levels=all_metrics)
correlate_subfamily$p.value = as.numeric(as.character(correlate_subfamily$p.value))
correlate_subfamily$p.adjust = p.adjust(correlate_subfamily$p.value,method="fdr")
correlate_subfamily$estimate.rho = as.numeric(as.character(correlate_subfamily$estimate.rho))
correlate_subfamily$class_update = factor(correlate_subfamily$class_update,levels=c("DNA","LINE","LTR","SINE","SVA","Other","All"))

# Correlation between input variables (including feature ratio)
correlate_subfam_feature = ldply(all_metrics,function(x) correlate_spearman(rmsk_TE_subfamily_measure,x,all_metrics))
correlate_subfam_feature$Metric = rep(all_metrics,each=7*length(all_metrics))
correlate_subfam_feature$p.value = as.numeric(as.character(correlate_subfam_feature$p.value))
correlate_subfam_feature$p.adjust = p.adjust(correlate_subfam_feature$p.value,method="fdr")
correlate_subfam_feature$estimate.rho = as.numeric(as.character(correlate_subfam_feature$estimate.rho))
correlate_subfam_feature$class_update = factor(correlate_subfam_feature$class_update,levels=c("DNA","LINE","LTR","SINE","SVA","Other","All"))

# Compare states to each other
correlate_subfam_state = ldply(all_pc_states,function(x) correlate_spearman(rmsk_TE_subfamily_measure,x,all_pc_states))
correlate_subfam_state$State2 = rep(all_pc_states,each=7*length(all_pc_states))
correlate_subfam_state$p.value = as.numeric(as.character(correlate_subfam_state$p.value))
correlate_subfam_state$p.adjust = p.adjust(correlate_subfam_state$p.value,method="fdr")
correlate_subfam_state$estimate.rho = as.numeric(as.character(correlate_subfam_state$estimate.rho))
correlate_subfam_state$class_update = factor(correlate_subfam_state$class_update,levels=c("DNA","LINE","LTR","SINE","SVA","Other","All"))