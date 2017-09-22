# Average number of samples in state by Refseq feature overlap
# Currently with IMR90 for methylation
# See 9/20/2016, 9/28/2016, 12/14/2016, 2/9/2017, 6/6/2017, 6/14/17, 7/24/17, 8/1/2017

source("R_scripts/TE_correlation.R")
source("R_scripts/potential.R")

# For TEs overlapping each feature, mean number of samples in state
feature_state_mean = apply(rmsk_TE_measure[,15:34],2,function(y) 
  (colMeans(rmsk_TE_measure[which(y == "yes"),c(35:49,52:55,60:62)],na.rm=TRUE)-colMeans(rmsk_TE_measure[which(y == "no"),c(35:49,52:55,60:62)],na.rm=TRUE))/colMeans(rmsk_TE_measure[,c(35:49,52:55,60:62)],na.rm=TRUE))

rownames(feature_state_mean)[1:15] = chromHMM_states
feature_state_mean[1:15,] = feature_state_mean[1:15,]/1.27
feature_state_mean[16:19,] = feature_state_mean[16:19,]/0.37
feature_state_mean[20,] = feature_state_mean[20,]/0.53
feature_state_mean[21,] = feature_state_mean[21,]/0.98
feature_state_mean[22,] = feature_state_mean[22,]/0.52

feature_state_mean = melt(as.matrix(feature_state_mean))
colnames(feature_state_mean) = c("State","Feature","Enrichment")
feature_state_mean$Feature = factor(feature_state_mean$Feature,levels=levels(feature_state_mean$Feature)[c(17,19,18,5,7,6,8:9,2,4,3,10,12,11,14,16,15,13,1,20)])

# By class
feature_state_mean_class = ddply(rmsk_TE_measure,~class_update,function(z) 
  apply(z[,15:34],2,function(y) 
    (colMeans(z[which(y == "yes"),c(35:49,52:55,60:62)],na.rm=TRUE)-colMeans(z[which(y == "no"),c(35:49,52:55,60:62)],na.rm=TRUE))/colMeans(z[,c(35:49,52:55,60:62)],na.rm=TRUE)))
feature_state_mean_class$State = factor(rep(c(chromHMM_states,colnames(rmsk_TE_measure)[c(52:55,60:62)]),6),levels=c(chromHMM_states,meth_states,"DNase","H3K27ac","Expressed_samples"))

feature_state_mean_class[which(feature_state_mean_class$State %in% chromHMM_states),2:21] = feature_state_mean_class[which(feature_state_mean_class$State %in% chromHMM_states),2:21]/1.27
feature_state_mean_class[which(feature_state_mean_class$State %in% meth_states),2:21] = feature_state_mean_class[which(feature_state_mean_class$State %in% meth_states),2:21]/0.37
feature_state_mean_class[which(feature_state_mean_class$State == "DNase"),2:21] = feature_state_mean_class[which(feature_state_mean_class$State == "DNase"),2:21]/0.53
feature_state_mean_class[which(feature_state_mean_class$State == "H3K27ac"),2:21] = feature_state_mean_class[which(feature_state_mean_class$State == "H3K27ac"),2:21]/0.98
feature_state_mean_class[which(feature_state_mean_class$State == "Expressed_samples"),2:21] = feature_state_mean_class[which(feature_state_mean_class$State == "Expressed_samples"),2:21]/0.52

feature_state_mean_class = melt(feature_state_mean_class)
colnames(feature_state_mean_class) = c("Class","State","Feature","Enrichment")
feature_state_mean_class$Feature = factor(feature_state_mean_class$Feature,levels=levels(feature_state_mean_class$Feature)[c(17,19,18,5,7,6,8:9,2,4,3,10,12,11,14,16,15,13,1,20)])

# Older code
#test = melt(as.matrix(apply(rmsk_TE_measure[,15:34],2,function(x) apply(rmsk_TE_measure[,35:63],2,function(y) {table = aggregate(data=rmsk_TE_measure,y~x,mean);(table$y[2]-table$y[1])/mean(y)}))))
#compare_TE_state = ddply(rmsk_TE_measure,~class_update,function(z) apply(z[,15:34],2,function(x) apply(z[,35:63],2,function(y) {table = aggregate(data=rmsk_TE_measure,y~x,mean);(table$y[2]-table$y[1])/mean(y)})))
