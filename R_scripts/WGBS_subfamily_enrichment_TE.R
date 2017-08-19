# Enrichment of hypomethylated TEs by subfamily and sample
# See 5/9/2016, 5/25/2016, 5/26/2016, 7/10/2016, 7/11/2016, 8/29/2016, 8/30/2016, 9/6/2016, 9/7/2016, 9/8/2016, 9/19/2016, 9/27/2016, 11/5/2016, 11/7/2016, 11/9/2016, 12/16/2016, 12/24/2016, 12/27/2016, 1/13/2017, 2/6/2017, 2/9/2017, 2/10/2017, 2/15/2017, 2/16/2017, 2/23/2017, 3/5/2017, 3/16/2017, 5/14/2017, 5/15/2017, 5/16/2017, 5/17/2017, 5/18/2017, 6/14/2017, 6/15/2017, 7/23/2017, 7/24/2017

load("R_scripts/TE_meth_average.RData")

# Proportion of subfamily members hypomethylated
TE_meth_subfamily_hypo = aggregate(data=TE_meth_average[,c(4:6,8:44,54)],.~subfamily+family+class_update,function(x) sum(na.omit(x) < 0.3)/length(na.omit(x)),na.action=na.pass)

# Average across all samples
TE_meth_subfamily_hypo$Mean = apply(TE_meth_subfamily_hypo,1,function(x) mean(as.numeric(x[4:40]),na.rm=TRUE))
TE_meth_subfamily_hypo$Mean_noIMR90 = apply(TE_meth_subfamily_hypo,1,function(x) mean(as.numeric(x[c(4:13,15:40)]),na.rm=TRUE))

# Range across all samples
TE_meth_subfamily_hypo$Range = apply(TE_meth_subfamily_hypo,1,function(x) max(as.numeric(na.omit(x[4:40])))-min(as.numeric(na.omit(x[4:40]))))
TE_meth_subfamily_hypo$Range_noIMR90 = apply(TE_meth_subfamily_hypo,1,function(x) max(as.numeric(na.omit(x[c(4:13,15:40)])))-min(as.numeric(na.omit(x[c(4:13,15:40)]))))

# Maximum across all samples 
TE_meth_subfamily_hypo$Max = apply(TE_meth_subfamily_hypo,1,function(x) max(as.numeric(x[4:40])))
TE_meth_subfamily_hypo$Max_noIMR90 = apply(TE_meth_subfamily_hypo,1,function(x) max(as.numeric(x[c(4:13,15:40)])))