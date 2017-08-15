# Enrichment of hypomethylated TEs by subfamily and sample
# See 5/9/2016, 5/25/2016, 5/26/2016, 7/10/2016, 7/11/2016, 8/29/2016, 8/30/2016, 9/6/2016, 9/7/2016, 9/8/2016, 9/19/2016, 9/27/2016, 11/5/2016, 11/7/2016, 11/9/2016, 12/16/2016, 12/24/2016, 12/27/2016, 1/13/2017, 2/6/2017, 2/9/2017, 2/10/2017, 2/15/2017, 2/16/2017, 2/23/2017, 3/5/2017, 3/16/2017, 5/14/2017, 5/15/2017, 5/16/2017, 5/17/2017, 5/18/2017, 6/14/2017, 6/15/2017, 7/23/2017, 7/24/2017

# Proportion of subfamily members hypomethylated
TE_meth_subfamily_hypo = aggregate(data=TE_meth_average[,c(4:6,8:44)],.~subfamily+class+family,function(x) sum(na.omit(x) < 0.3)/length(na.omit(x)),na.action=na.pass)

# Average across all samples
TE_meth_subfamily_hypo$Mean = apply(TE_meth_subfamily_hypo,1,function(x) mean(as.numeric(x[4:40]),na.rm=TRUE))
TE_meth_subfamily_hypo$Mean_noIMR90 = apply(TE_meth_subfamily_hypo,1,function(x) mean(as.numeric(x[c(4:13,15:40)]),na.rm=TRUE))

# Range across all samples
TE_meth_subfamily_hypo$Range = apply(TE_meth_subfamily_hypo,1,function(x) max(as.numeric(na.omit(x[4:40])))-min(as.numeric(na.omit(x[4:40]))))
TE_meth_subfamily_hypo$Range_noIMR90 = apply(TE_meth_subfamily_hypo,1,function(x) max(as.numeric(na.omit(x[c(4:13,15:40)])))-min(as.numeric(na.omit(x[c(4:13,15:40)]))))
mean(TE_meth_subfamily_hypo$Range)
mean(TE_meth_subfamily_hypo$Range_noIMR90)
mean(TE_meth_subfamily_hypo$Mean)
mean(TE_meth_subfamily_hypo$Mean_noIMR90)

# Maximum across all samples 
TE_meth_subfamily_hypo$Max = apply(TE_meth_subfamily_hypo,1,function(x) max(as.numeric(x[4:40])))
TE_meth_subfamily_hypo$Max_noIMR90 = apply(TE_meth_subfamily_hypo,1,function(x) max(as.numeric(x[c(4:13,15:40)])))

# Proportion of all hypomethylated TEs per subfamily x sample
TE_meth_subfamily_hypo_num = aggregate(data=TE_meth_average[,c(4:6,8:44)],.~subfamily+class+family,function(x) sum(na.omit(x) < 0.3),na.action=na.pass)
TE_meth_subfamily_hypo_num[,4:40] = apply(TE_meth_subfamily_hypo_num[,4:40],2,function(x) x/sum(x))

# Subfamilies with >1% hypomethylated members
TE_meth_subfamily_hypo_num$Samples_1per = apply(TE_meth_subfamily_hypo_num,1,function(x) sum(as.numeric(na.omit(x[c(4:40)])) > 0.01))

# Average across all samples, with/without IMR90
mean(unlist(TE_meth_subfamily_hypo_num[,c(4:13,15:40)]),na.rm=TRUE)
mean(unlist(TE_meth_subfamily_hypo_num[,4:40]),na.rm=TRUE)

TE_meth_subfamily_hypo_num_long = melt(TE_meth_subfamily_hypo_num[,1:40],id.vars=c("class","family","subfamily"))
colnames(TE_meth_subfamily_hypo_num_long)[4:5] = c("Sample","Proportion")
TE_meth_subfamily_hypo_num_long = merge(TE_meth_subfamily_hypo_num_long,EID_metadata[,c(1,5,7:11)],by=c("Sample"))

# Enrichment (LOR) of hypomethylated TEs in a subfamily
TE_meth_subfamily_hypo_lor = TE_meth_subfamily_hypo[,1:40]
TE_meth_subfamily_hypo_lor[,4:40] = t(apply(TE_meth_subfamily_hypo_lor,1,function(x) log2(as.numeric(x[4:40])/TE_meth_average_state[colnames(TE_meth_subfamily_hypo_lor[,4:40]),]$Hypomethylated)))
TE_meth_subfamily_hypo_lor$Samples_1.5 = apply(TE_meth_subfamily_hypo_lor,1,function(x) sum(as.numeric(na.omit(x[4:40])) > 1.5)) 
dim(TE_meth_subfamily_hypo_lor[which(!(TE_meth_subfamily_hypo_lor$subfamily %in% TEother_meth_exclude) & TE_meth_subfamily_hypo_lor$Samples_1.5 > 0),])
sum(na.omit(unlist(TE_meth_subfamily_hypo_lor[which(!(TE_meth_subfamily_hypo_lor$subfamily %in% TEother_meth_exclude)),4:40])) > 1.5)

# Adding Unconfident class
TE_meth_subfamily_hypo_lor$class_update = TE_meth_subfamily_hypo_lor$class
TE_meth_subfamily_hypo_lor$class_update = factor(TE_meth_subfamily_hypo_lor$class_update,levels=c("DNA","LINE","LTR","SINE","RC","Other","Unconfident"))
TE_meth_subfamily_hypo_lor[which(TE_meth_subfamily_hypo_lor$class %in% c("DNA?","LINE?","LTR?","SINE?","Unknown","Unknown?")),]$class_update = "Unconfident"

TE_meth_subfamily_hypo_lor_long = melt(TE_meth_subfamily_hypo_lor[,c(1:40,42)],id.vars=c("subfamily","class","family","class_update"))
colnames(TE_meth_subfamily_hypo_lor_long)[5:6] = c("Sample","Enrichment")
TE_meth_subfamily_hypo_lor_long = merge(TE_meth_subfamily_hypo_lor_long,EID_metadata[,c(1,5,7:11)],by=c("Sample"),all.x=TRUE)
