# Potential for WGBS 
# See 6/2/2016, 7/11/2016, 8/25/2016, 8/26/2016, 9/16/2016, 9/20/2016, 9/28/2016, 12/15/2016, 2/6/2017, 5/10/2017, 7/22/2017, 8/1/2017

# Cumulative distribution of methylation states and statistics
TE_meth_average_category = sample_distribution(TE_meth_average,c(46:49),37)
TE_meth_average_category_cum = cumulative_distribution(TE_meth_average,c(46:49),37)
TE_meth_average_category_cum_long = melt(as.matrix(TE_meth_average_category_cum))
colnames(TE_meth_average_category_cum_long) = c("Samples","Methylation","Proportion")
TE_meth_average_category_stats = as.data.frame(t(rbind(apply(TE_meth_average_category[,2:5],2,function(x) sum(x[2:38])/3200428),apply(TE_meth_average_category[,2:5],2,function(x) sum(as.numeric(x)*seq(0,37))/sum(x))/37,apply(TE_meth_average_category[2:38,2:5],2,function(x) sum(as.numeric(x)*seq(1,37))/sum(x))/37)))
TE_meth_average_category_stats$State = rownames(TE_meth_average_category_stats)
colnames(TE_meth_average_category_stats)[1:3] = c("Proportion_ever","Samples_avg_all","Samples_avg_ever")
TE_meth_average_category_stats$State = factor(TE_meth_average_category_stats$State,levels=as.vector(TE_meth_average_category_stats$State)[c(4,2,3,1)])

# Cumulative distribution of methylation states and statistics, no IMR90
TE_meth_average_noIMR90_category = sample_distribution(TE_meth_average,c(50:53),36)
TE_meth_average_noIMR90_category_cum = cumulative_distribution(TE_meth_average,c(50:53),36)
TE_meth_average_noIMR90_category_cum_long = melt(as.matrix(TE_meth_average_noIMR90_category_cum))
colnames(TE_meth_average_noIMR90_category_cum_long) = c("Samples","Methylation","Proportion")
TE_meth_average_noIMR90_category_stats = as.data.frame(t(rbind(apply(TE_meth_average_noIMR90_category[,2:5],2,function(x) sum(x[2:37])/3200428),apply(TE_meth_average_noIMR90_category[,2:5],2,function(x) sum(as.numeric(x)*seq(0,36))/sum(x))/36,apply(TE_meth_average_noIMR90_category[2:37,2:5],2,function(x) sum(as.numeric(x)*seq(1,36))/sum(x))/36)))
TE_meth_average_noIMR90_category_stats$State = factor(c("Hypomethylated","Hypermethylated","Intermediate","Missing"),levels=c("Missing","Hypermethylated","Intermediate","Hypomethylated"))
colnames(TE_meth_average_noIMR90_category_stats)[1:3] = c("Proportion_ever","Samples_avg_all","Samples_avg_ever")