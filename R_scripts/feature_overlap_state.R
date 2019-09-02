# Calculates the proportion of increase in the mean proportion of samples each TE is annotated with each epigenetic state,
# for TEs overlapping each genic feature versus those that do not

## feature_state_mean - proportion of increase in samples annotated with each state based on genic feature overlap, for all TEs
## feature_state_mean_class - proportion of increase in samples annotated with each state based on genic feature overlap, by class

# Proportion of increase in the mean proportion of samples each TE is annotated with each epigenetic state,
# For TEs overlapping each genic feature versus those that do not
feature_state_mean = apply(rmsk_TE_measure[,cohorts],2,function(y) 
  (colMeans(rmsk_TE_measure[which(y > 0),states],na.rm=TRUE)-colMeans(rmsk_TE_measure[which(y == 0),states],na.rm=TRUE))/colMeans(rmsk_TE_measure[which(y == 0),states],na.rm=TRUE))

rownames(feature_state_mean)[22] = "Expressed_samples"

feature_state_mean = melt(as.matrix(feature_state_mean))
colnames(feature_state_mean) = c("State","Cohort","Enrichment")
feature_state_mean$State = factor(feature_state_mean$State,levels=states)
## Splits genic feature and coding status
feature_state_mean = split_coding(feature_state_mean)

# Proportion of increase in the mean proportion of samples each TE is annotated with each epigenetic state,
# For TEs overlapping each genic feature versus those that do not, by class
feature_state_mean_class = ddply(rmsk_TE_measure,~class_update,function(z) 
  apply(z[,cohorts],2,function(y) 
    (colMeans(z[which(y > 0),states],na.rm=TRUE)-colMeans(z[which(y == 0),states],na.rm=TRUE))/colMeans(z[which(y == 0),states],na.rm=TRUE)))
feature_state_mean_class$State = factor(rep(c(chromHMM_states,meth_states[c(1,3,2,4)],"DNase","H3K27ac","Expressed_samples"),6),levels=states)

feature_state_mean_class = melt(feature_state_mean_class)
colnames(feature_state_mean_class) = c("Class","State","Cohort","Enrichment")
## Splits genic feature and coding status
feature_state_mean_class = split_coding(feature_state_mean_class)