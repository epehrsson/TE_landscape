# Spearman correlation between states and metrics
# Individual TE level

#source("R_scripts/TE_correlation.R")

a = length(c(states,measure_states_extra))

# Correlation between TE metrics and number of samples in state
correlate_TE_state = rbind(correlate_spearman(rmsk_TE_measure,"JC_distance",c(states,measure_states_extra)), 
                           correlate_spearman(rmsk_TE_measure,"Length",c(states,measure_states_extra)),  
                           correlate_spearman(rmsk_TE_measure,"mappability",c(states,measure_states_extra)), 
                           correlate_spearman(rmsk_TE_measure,"CpGs",c(states,measure_states_extra)), 
                           correlate_spearman(rmsk_TE_measure,"CpGs_per_length",c(states,measure_states_extra)))
correlate_TE_state$Metric = c(rep("JC_distance",a*7),rep("Length",a*7),rep("mappability",a*7),rep("CpGs",a*7),rep("CpGs_per_length",a*7))
correlate_TE_state$p.value = as.numeric(as.character(correlate_TE_state$p.value))
correlate_TE_state$estimate.rho = as.numeric(as.character(correlate_TE_state$estimate.rho))
correlate_TE_state$State = factor(correlate_TE_state$State,levels=c(states,measure_states_extra))
correlate_TE_state$class_update = factor(correlate_TE_state$class_update,levels=c("DNA","LINE","LTR","SINE","SVA","Other","All"))

# Correlation between TE metrics
correlate_TE_feature = rbind(correlate_spearman(rmsk_TE_measure,"Length",measure_metrics[2:5]),
                             correlate_spearman(rmsk_TE_measure,"mappability",measure_metrics[3:5]), 
                             correlate_spearman(rmsk_TE_measure,"JC_distance",measure_metrics[4:5]),
                             correlate_spearman(rmsk_TE_measure,"CpGs",measure_metrics[5]))
correlate_TE_feature$Metric = factor(c(rep("Length",28),rep("mappability",21),rep("JC_distance",14),rep("CpGs",14)),levels=c("Length","mappability","JC_distance","CpGs","CpGs_per_length"))
correlate_TE_feature$p.value = as.numeric(as.character(correlate_TE_feature$p.value))
correlate_TE_feature$estimate.rho = as.numeric(as.character(correlate_TE_feature$estimate.rho))
correlate_TE_feature$class_update = factor(correlate_TE_feature$class_update,levels=c("DNA","LINE","LTR","SINE","SVA","Other","All"))
correlate_TE_feature$State = factor(correlate_TE_feature$State,levels=measure_metrics)
correlate_TE_feature = correlate_TE_feature[which(correlate_TE_feature$State != correlate_TE_feature$Metric),]
