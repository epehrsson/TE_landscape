# Proportion of subfamily members hypermethylated
# See 9/7/2016, 9/8/2017, 12/16/2016, 2/6/2017, 3/16/2017, 7/23/2017

TE_meth_subfamily_hyper = aggregate(data=TE_meth_average[,c(4:6,8:44)],.~subfamily+class+family,function(x) sum(na.omit(x) > 0.7)/length(na.omit(x)),na.action=na.pass)

# Subfamilies hypermethylated in all samples
sort(apply(TE_meth_subfamily_hyper[,c(4:13,15:40)],1,function(x) sum(na.omit(x) == 1)))

# Average proportion hypermethylated, with/without IMR90
mean(unlist(TE_meth_subfamily_hyper[,c(4:13,15:40)]),na.rm=TRUE)
mean(unlist(TE_meth_subfamily_hyper[,c(4:40)]),na.rm=TRUE)
