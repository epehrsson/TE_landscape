# Methylation state of TEs by sample
# See 6/1/2016, 7/10/2016, 9/16/2016, 9/17/2016, 12/15/2016, 12/16/2016, 1/24/2017, 2/6/2017, 2/9/2017, 3/8/2017, 5/11/2017, 5/12/2017, 6/14/2017, 7/22/2017, 8/1/2017

# Average TE methylation level by sample
TE_meth_average_mean = apply(TE_meth_average[,c(8:44)],2,function(x) mean(na.omit(x)))

# Summaries with/without IMR90
mean(TE_meth_average_mean)
sd(TE_meth_average_mean)
mean(TE_meth_average_mean[c(1:10,12:37)])
sd(TE_meth_average_mean[c(1:10,12:37)])

# TEs per methylation state per sample
TE_meth_average_state = as.data.frame(cbind(apply(TE_meth_average[,8:44],2,function(x) sum(na.omit(x) < 0.3)/length(x)),apply(TE_meth_average[,8:44],2,function(x) sum(na.omit(x) > 0.7)/length(x)),apply(TE_meth_average[,8:44],2,function(x)  sum(na.omit(x) <= 0.7 & na.omit(x) >= 0.3)/length(x)),apply(TE_meth_average[,8:44],2,function(x) sum(is.na(x))/length(x))))
colnames(TE_meth_average_state) = c("Hypomethylated","Hypermethylated","Intermediate","Missing")
TE_meth_average_state = TE_meth_average_state[order(TE_meth_average_state$Hypermethylated + TE_meth_average_state$Missing),]

# Average proportion of TEs hypomethylated, excluding IMR90
mean(TE_meth_average_state$Hypomethylated[2:37])
sd(TE_meth_average_state$Hypomethylated[2:37])

# Proportion of TEs hypomethylated by sample metadata
TE_meth_average_state_long = melt(as.matrix(TE_meth_average_state))
colnames(TE_meth_average_state_long) = c("Sample","State","Proportion")
TE_meth_average_state_long = merge(TE_meth_average_state_long,EID_metadata[,c(1,5,7:11)],by=c("Sample"))
TE_meth_average_state_long$State = factor(TE_meth_average_state_long$State,levels=levels(TE_meth_average_state_long$State)[c(4,2,3,1)])

# Hypomethylated TEs by sample grouping
kruskal.test(TE_meth_average_state_long[which(TE_meth_average_state_long$State == "Hypomethylated"),]$Proportion,TE_meth_average_state_long[which(TE_meth_average_state_long$State == "Hypomethylated"),]$Type)
kruskal.test(TE_meth_average_state_long[which(TE_meth_average_state_long$State == "Hypomethylated"),]$Proportion,TE_meth_average_state_long[which(TE_meth_average_state_long$State == "Hypomethylated"),]$Group)
kruskal.test(TE_meth_average_state_long[which(TE_meth_average_state_long$State == "Hypomethylated"),]$Proportion,TE_meth_average_state_long[which(TE_meth_average_state_long$State == "Hypomethylated"),]$Anatomy)
kruskal.test(TE_meth_average_state_long[which(TE_meth_average_state_long$State == "Hypomethylated"),]$Proportion,TE_meth_average_state_long[which(TE_meth_average_state_long$State == "Hypomethylated"),]$Germline)
kruskal.test(TE_meth_average_state_long[which(TE_meth_average_state_long$State == "Hypomethylated"),]$Proportion,TE_meth_average_state_long[which(TE_meth_average_state_long$State == "Hypomethylated"),]$Age)
apply(TE_meth_average_state,2,function(x) wilcox_to_all(x,droplevels(EID_metadata[match(rownames(TE_meth_average_state),EID_metadata$Sample),]$Group)))

# Hypomethylated TEs by sample grouping, no IMR90
kruskal.test(TE_meth_average_state_long[which(TE_meth_average_state_long$State == "Hypomethylated" & TE_meth_average_state_long$Sample != "E017"),]$Proportion,TE_meth_average_state_long[which(TE_meth_average_state_long$State == "Hypomethylated" & TE_meth_average_state_long$Sample != "E017"),]$Type)
kruskal.test(TE_meth_average_state_long[which(TE_meth_average_state_long$State == "Hypomethylated" & TE_meth_average_state_long$Sample != "E017"),]$Proportion,TE_meth_average_state_long[which(TE_meth_average_state_long$State == "Hypomethylated" & TE_meth_average_state_long$Sample != "E017"),]$Group)
kruskal.test(TE_meth_average_state_long[which(TE_meth_average_state_long$State == "Hypomethylated" & TE_meth_average_state_long$Sample != "E017"),]$Proportion,TE_meth_average_state_long[which(TE_meth_average_state_long$State == "Hypomethylated" & TE_meth_average_state_long$Sample != "E017"),]$Anatomy)
kruskal.test(TE_meth_average_state_long[which(TE_meth_average_state_long$State == "Hypomethylated" & TE_meth_average_state_long$Sample != "E017"),]$Proportion,TE_meth_average_state_long[which(TE_meth_average_state_long$State == "Hypomethylated" & TE_meth_average_state_long$Sample != "E017"),]$Germline)
kruskal.test(TE_meth_average_state_long[which(TE_meth_average_state_long$State == "Hypomethylated" & TE_meth_average_state_long$Sample != "E017"),]$Proportion,TE_meth_average_state_long[which(TE_meth_average_state_long$State == "Hypomethylated" & TE_meth_average_state_long$Sample != "E017"),]$Age)
apply(TE_meth_average_state[which(rownames(TE_meth_average_state) != "E017"),],2,function(x) wilcox_to_all(x,droplevels(EID_metadata[match(rownames(TE_meth_average_state[which(rownames(TE_meth_average_state) != "E017"),]),EID_metadata$Sample),]$Group)))
apply(TE_meth_average_state[which(rownames(TE_meth_average_state) != "E017"),],2,function(x) wilcox_to_all(x,droplevels(EID_metadata[match(rownames(TE_meth_average_state[which(rownames(TE_meth_average_state) != "E017"),]),EID_metadata$Sample),]$Age)))
apply(TE_meth_average_state[which(rownames(TE_meth_average_state) != "E017"),],2,function(x) wilcox_to_all(x,droplevels(EID_metadata[match(rownames(TE_meth_average_state[which(rownames(TE_meth_average_state) != "E017"),]),EID_metadata$Sample),]$Germline)))
apply(TE_meth_average_state[which(rownames(TE_meth_average_state) != "E017"),],2,function(x) wilcox_to_all(x,droplevels(EID_metadata[match(rownames(TE_meth_average_state[which(rownames(TE_meth_average_state) != "E017"),]),EID_metadata$Sample),]$Type)))
