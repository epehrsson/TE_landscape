# Potential for WGBS 
# See 7/22/2017, 8/1/2017

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

# Cumulative distribution of methylation states and statistics, by class
TE_meth_average_class = by(TE_meth_average,TE_meth_average$class_update,function(x) sample_distribution(x,c(46:49),37))
TE_meth_average_class_cum = by(TE_meth_average,TE_meth_average$class_update,function(x) cumulative_distribution(x,c(46:49),37))
TE_meth_average_class_cum_long = rbind(melt(as.matrix(TE_meth_average_class_cum$DNA)),melt(as.matrix(TE_meth_average_class_cum$LINE)),melt(as.matrix(TE_meth_average_class_cum$LTR)),melt(as.matrix(TE_meth_average_class_cum$SINE)),melt(as.matrix(TE_meth_average_class_cum$Other)),melt(as.matrix(TE_meth_average_class_cum$RC)),melt(as.matrix(TE_meth_average_class_cum$Unconfident)))
TE_meth_average_class_cum_long$Class = c(rep("DNA",148),rep("LINE",148),rep("LTR",148),rep("SINE",148),rep("Other",148),rep("RC",148),rep("Unconfident",148))
colnames(TE_meth_average_class_cum_long)[1:3] = c("Samples","Methylation","Proportion")
TE_meth_average_class_cum_long$Group = paste(TE_meth_average_class_cum_long$Methylation,TE_meth_average_class_cum_long$Class,sep="_")
TE_meth_average_class_stats = ldply(TE_meth_average_class,function(y) as.data.frame(t(rbind(apply(y[,2:5],2,function(x) sum(x[2:38])/(sum(x)/100)),apply(y[,2:5],2,function(x) sum(as.numeric(x)*seq(0,37))/sum(x))/37,apply(y[2:38,2:5],2,function(x) sum(as.numeric(x)*seq(1,37))/sum(x))/37))))
colnames(TE_meth_average_class_stats) = c("Class","Proportion_ever","Samples_avg_all","Samples_avg_ever")
TE_meth_average_class_stats$Class = factor(TE_meth_average_class_stats$Class,levels=c("DNA","LINE","LTR","SINE","Other","RC","Unconfident"))
TE_meth_average_class_stats$State = factor(rep(meth_states[c(1,3,2,4)],7),levels=meth_states)
TE_meth_average_class_stats[,2:4] = apply(TE_meth_average_class_stats[,2:4],2,function(x) as.numeric(x))

# Cumulative distribution of methylation states and statistics, no IMR90, by class
TE_meth_average_class_noIMR90 = by(TE_meth_average,TE_meth_average$class_update,function(x) sample_distribution(x,c(50:53),36))
TE_meth_average_class_noIMR90_cum = by(TE_meth_average,TE_meth_average$class_update,function(x) cumulative_distribution(x,c(50:53),36))
TE_meth_average_class_noIMR90_cum_long = rbind(melt(as.matrix(TE_meth_average_class_noIMR90_cum$DNA)),melt(as.matrix(TE_meth_average_class_noIMR90_cum$LINE)),melt(as.matrix(TE_meth_average_class_noIMR90_cum$LTR)),melt(as.matrix(TE_meth_average_class_noIMR90_cum$SINE)),melt(as.matrix(TE_meth_average_class_noIMR90_cum$Other)),melt(as.matrix(TE_meth_average_class_noIMR90_cum$RC)),melt(as.matrix(TE_meth_average_class_noIMR90_cum$Unconfident)))
TE_meth_average_class_noIMR90_cum_long$Class = c(rep("DNA",144),rep("LINE",144),rep("LTR",144),rep("SINE",144),rep("Other",144),rep("RC",144),rep("Unconfident",144))
colnames(TE_meth_average_class_noIMR90_cum_long)[1:3] = c("Samples","Methylation","Proportion")
TE_meth_average_class_noIMR90_cum_long$Group = paste(TE_meth_average_class_noIMR90_cum_long$Methylation,TE_meth_average_class_noIMR90_cum_long$Class,sep="_")
TE_meth_average_class_noIMR90_stats = ldply(TE_meth_average_class_noIMR90,function(y) as.data.frame(t(rbind(apply(y[,2:5],2,function(x) sum(x[2:37])/(sum(x)/100)),apply(y[,2:5],2,function(x) sum(as.numeric(x)*seq(0,36))/sum(x))/36,apply(y[2:37,2:5],2,function(x) sum(as.numeric(x)*seq(1,36))/sum(x))/36))))
colnames(TE_meth_average_class_noIMR90_stats) = c("Class","Proportion_ever","Samples_avg_all","Samples_avg_ever")
TE_meth_average_class_noIMR90_stats$Class = factor(TE_meth_average_class_noIMR90_stats$Class,levels=c("DNA","LINE","LTR","SINE","Other","RC","Unconfident"))
TE_meth_average_class_noIMR90_stats$State = factor(rep(meth_states[c(1,3,2,4)],7),levels=meth_states)
TE_meth_average_class_noIMR90_stats[,2:4] = apply(TE_meth_average_class_noIMR90_stats[,2:4],2,function(x) as.numeric(x))
