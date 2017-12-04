# Average TE metrics by Refseq feature overlap

#source("R_scripts/TE_correlation.R")

# For TEs overlapping each feature, metric average
feature_metric_mean = melt(as.matrix(apply(rmsk_TE_measure[,15:34],2,function(y) 
  (colMeans(rmsk_TE_measure[which(y == "yes"),c(8,10,12:14)],na.rm=TRUE)-colMeans(rmsk_TE_measure[which(y == "no"),c(8,10,12:14)],na.rm=TRUE))/colMeans(rmsk_TE_measure[which(y == "no"),c(8,10,12:14)],na.rm=TRUE))))

colnames(feature_metric_mean) = c("Metric","Cohort","Enrichment")
feature_metric_mean = split_coding(feature_metric_mean,2)

# By class
feature_metric_mean_class = ddply(rmsk_TE_measure,~class_update,function(z) 
  apply(z[,15:34],2,function(y) 
    (colMeans(z[which(y == "yes"),c(8,10,12:14)],na.rm=TRUE)-colMeans(z[which(y == "no"),c(8,10,12:14)],na.rm=TRUE))/colMeans(z[which(y == "no"),c(8,10,12:14)],na.rm=TRUE)))
feature_metric_mean_class$Metric = rep(colnames(rmsk_TE_measure)[c(8,10,12:14)],6)
feature_metric_mean_class = melt(feature_metric_mean_class,id.vars=c("class_update","Metric"))
colnames(feature_metric_mean_class) = c("Class","Metric","Cohort","Enrichment")
feature_metric_mean_class = split_coding(feature_metric_mean_class,3)

# Older code
#test = melt(as.matrix(apply(rmsk_TE_measure[,15:34],2,function(x) apply(rmsk_TE_measure[,c(8,10,12:14)],2,function(y) {table = aggregate(data=rmsk_TE_measure,y~x,mean);(table$y[2]-table$y[1])/mean(y)}))))

#compare_TE_feature = ddply(rmsk_TE_measure,~class_update,function(z) apply(z[,15:34],2,function(x) apply(z[,c(8,10,12:14)],2,function(y) {table = aggregate(data=z,y~x,mean);(table$y[2]-table$y[1])/mean(y)})))
#compare_TE_feature$Metric = rep(colnames(rmsk_TE_measure)[c(8,10,12:14)],6)