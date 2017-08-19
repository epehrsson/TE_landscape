# Proportion of subfamily members hypermethylated
# See 9/7/2016, 9/8/2017, 12/16/2016, 2/6/2017, 3/16/2017, 7/23/2017

load("R_scripts/TE_meth_average.RData")

TE_meth_subfamily_hyper = aggregate(data=TE_meth_average[,c(4:6,8:44,54)],.~subfamily+family+class_update,function(x) sum(na.omit(x) > 0.7)/length(na.omit(x)),na.action=na.pass)