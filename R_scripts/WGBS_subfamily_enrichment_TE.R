# Enrichment of methylation state TEs by subfamily and sample
# See 5/9/2016, 5/25/2016, 5/26/2016, 7/10/2016, 7/11/2016, 8/29/2016, 8/30/2016, 9/6/2016, 9/7/2016, 9/8/2016, 9/19/2016, 9/27/2016, 11/5/2016, 11/7/2016, 11/9/2016, 12/16/2016, 12/24/2016, 12/27/2016, 1/13/2017, 2/6/2017, 2/9/2017, 2/10/2017, 2/15/2017, 2/16/2017, 2/23/2017, 3/5/2017, 3/16/2017, 5/14/2017, 5/15/2017, 5/16/2017, 5/17/2017, 5/18/2017, 6/14/2017, 6/15/2017, 7/23/2017, 7/24/2017

# Proportion of subfamily members in each methylation state per sample
TE_meth_subfamily = list(aggregate(data=TE_meth_average[,c(4,6,8:44,54)],.~subfamily+family+class_update,function(x) sum(na.omit(x) < 0.3),na.action=na.pass),
                              aggregate(data=TE_meth_average[,c(4,6,8:44,54)],.~subfamily+family+class_update,function(x) sum(na.omit(x) >= 0.3 & na.omit(x) <= 0.7),na.action=na.pass),
                              aggregate(data=TE_meth_average[,c(4,6,8:44,54)],.~subfamily+family+class_update,function(x) sum(na.omit(x) > 0.7),na.action=na.pass),
                              aggregate(data=TE_meth_average[,c(4,6,8:44,54)],.~subfamily+family+class_update,function(x) sum(is.na(x)),na.action=na.pass))
names(TE_meth_subfamily) = meth_states