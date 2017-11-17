# Spearman correlation between states and metrics
# Individual TE level

#source("R_scripts/TE_correlation.R")

# Correlation between TE metrics and number of samples in state
correlate_TE_state = rbind(correlate_spearman(rmsk_TE_measure,"JC_distance",c(35:63)), 
                           correlate_spearman(rmsk_TE_measure,"Length",c(35:63)),  
                           correlate_spearman(rmsk_TE_measure,"mappability",c(35:63)), 
                           correlate_spearman(rmsk_TE_measure,"CpGs",c(35:63)), 
                           correlate_spearman(rmsk_TE_measure,"CpGs_per_length",c(35:63)))
correlate_TE_state$Metric = c(rep("JC_distance",203),rep("Length",203),rep("mappability",203),rep("CpGs",203),rep("CpGs_per_length",203))
correlate_TE_state$p.value = as.numeric(as.character(correlate_TE_state$p.value))
correlate_TE_state$estimate.rho = as.numeric(as.character(correlate_TE_state$estimate.rho))
correlate_TE_state$State = factor(correlate_TE_state$State,levels=unique(correlate_TE_state$State)[c(1:17,18,22,20,24,19,23,21,25:29)])
correlate_TE_state$class_update = factor(correlate_TE_state$class_update,levels=c("DNA","LINE","LTR","SINE","SVA","Other","All"))

# Correlation between TE metrics
correlate_TE_feature = rbind(correlate_spearman(rmsk_TE_measure,"Length",c(10,12:14)),
                             correlate_spearman(rmsk_TE_measure,"mappability",c(12:14)), 
                             correlate_spearman(rmsk_TE_measure,"JC_distance",c(13:14)),
                             correlate_spearman(rmsk_TE_measure,"CpGs",c(13:14)))
correlate_TE_feature$Metric = factor(c(rep("Length",28),rep("mappability",21),rep("JC_distance",14),rep("CpGs",14)),levels=c("Length","mappability","JC_distance","CpGs","CpGs_per_length"))
correlate_TE_feature$p.value = as.numeric(as.character(correlate_TE_feature$p.value))
correlate_TE_feature$estimate.rho = as.numeric(as.character(correlate_TE_feature$estimate.rho))
correlate_TE_feature$class_update = factor(correlate_TE_feature$class_update,levels=c("DNA","LINE","LTR","SINE","SVA","Other","All"))
correlate_TE_feature$State = factor(correlate_TE_feature$State,levels=c("Length","mappability","JC_distance","CpGs","CpGs_per_length"))
correlate_TE_feature = correlate_TE_feature[which(correlate_TE_feature$State != correlate_TE_feature$Metric),]
