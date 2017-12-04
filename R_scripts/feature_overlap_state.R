# Average number of samples in state by Refseq feature overlap
# Currently with IMR90 for methylation
# See 9/20/2016, 9/28/2016, 12/14/2016, 2/9/2017, 6/6/2017, 6/14/17, 7/24/17, 8/1/2017

#source("R_scripts/TE_correlation.R")
#source("R_scripts/potential.R")

# For TEs overlapping each feature, mean number of samples in state
feature_state_mean = apply(rmsk_TE_measure[,15:34],2,function(y) 
  (colMeans(rmsk_TE_measure[which(y == "yes"),c(35:49,52:55,60:62)],na.rm=TRUE)-colMeans(rmsk_TE_measure[which(y == "no"),c(35:49,52:55,60:62)],na.rm=TRUE))/colMeans(rmsk_TE_measure[which(y == "no"),c(35:49,52:55,60:62)],na.rm=TRUE))

rownames(feature_state_mean)[1:15] = chromHMM_states
rownames(feature_state_mean)[22] = "Expression"

feature_state_mean = melt(as.matrix(feature_state_mean))
colnames(feature_state_mean) = c("State","Cohort","Enrichment")
feature_state_mean$State = factor(feature_state_mean$State,levels=c(chromHMM_states,meth_states,"DNase","H3K27ac","Expression"))
feature_state_mean = split_coding(feature_state_mean,2)

# By class
feature_state_mean_class = ddply(rmsk_TE_measure,~class_update,function(z) 
  apply(z[,15:34],2,function(y) 
    (colMeans(z[which(y == "yes"),c(35:49,52:55,60:62)],na.rm=TRUE)-colMeans(z[which(y == "no"),c(35:49,52:55,60:62)],na.rm=TRUE))/colMeans(z[which(y == "no"),c(35:49,52:55,60:62)],na.rm=TRUE)))
feature_state_mean_class$State = factor(rep(c(chromHMM_states,meth_states[c(1,3,2,4)],"DNase","H3K27ac","Expression"),6),levels=c(chromHMM_states,meth_states,"DNase","H3K27ac","Expression"))

feature_state_mean_class = melt(feature_state_mean_class)
colnames(feature_state_mean_class) = c("Class","State","Cohort","Enrichment")
feature_state_mean_class = split_coding(feature_state_mean_class,3)

# Older code
#test = melt(as.matrix(apply(rmsk_TE_measure[,15:34],2,function(x) apply(rmsk_TE_measure[,35:63],2,function(y) {table = aggregate(data=rmsk_TE_measure,y~x,mean);(table$y[2]-table$y[1])/mean(y)}))))
#compare_TE_state = ddply(rmsk_TE_measure,~class_update,function(z) apply(z[,15:34],2,function(x) apply(z[,35:63],2,function(y) {table = aggregate(data=rmsk_TE_measure,y~x,mean);(table$y[2]-table$y[1])/mean(y)})))
