# Proportion of feature in each chromHMM state
# See 4/18/2016, 4/20/2016, 4/25/2016, 4/26/2016, 4/27/2016, 5/10/2016, 6/27/2016, 9/8/2016, 9/25/2016, 2/3/2017, 2/10/2017, 3/8/2017, 5/8/2017, 8/7/2017

# Genome
# Number of bases in each state in each sample
t(read.table("mnemonics_state.txt",sep='\t',header=TRUE,row.names=1))

# Number of bp in each state and overall
colSums(mnemonics_states)
mnemonics_states_long = melt(mnemonics_states)
colnames(mnemonics_states_long) = c("Sample","State","Length")

# Number of bases in each state in each sample, normalized by total annotated
as.data.frame(mnemonics_states)[,2:16]/as.data.frame(mnemonics_states)$Total
test = melt(mnemonics_states_normalized)
colnames(test) = c("State","Proportion")
test$Cohort = rep("Genome",1905)

# Mean, median, sd of proportion of bases in each state
cbind(colMeans(mnemonics_states_normalized),apply(mnemonics_states_normalized,2,median),apply(mnemonics_states_normalized,2,sd))
colnames(test) = c("Mean","Median","SD")

# TEs
# Number of bases in each state in each sample in merged TEs
test = as.data.frame(t(read.table("mnemonics_TEmerge_states.txt",sep='\t',header=TRUE,row.names=1)))
test[is.na(test)] = 0

# Number of bp in each state and overall in merged TEs
colSums(mnemonics_states_TEmerge)

# Number of bases in each state in each sample in merged TEs, normalized by total annotated
mnemonics_states_TEmerge_normalized = mnemonics_states_TEmerge[2:16]/mnemonics_states_TEmerge$Total
mnemonics_states_TEmerge_normalized_long = melt(as.matrix(mnemonics_states_TEmerge_normalized))
mnemonics_states_TEmerge_normalized_long = cbind(mnemonics_states_TEmerge_normalized_long[,2:3],rep("TEmerge",1905))
colnames(mnemonics_states_TEmerge_normalized_long) = c("State","Proportion","Cohort")

# Mean, median, sd of proportion of bases in each state in merged TEs
mnemonics_states_TEmerge_normalized_average = cbind(colMeans(mnemonics_states_TEmerge_normalized),apply(mnemonics_states_TEmerge_normalized,2,median),apply(mnemonics_states_TEmerge_normalized,2,sd))
colnames(mnemonics_states_TEmerge_normalized_average) = c("Mean","Median","SD")

# Number of bases in each state in each sample in merged all TEs
mnemonics_states_TEother = as.data.frame(t(read.table("other/mnemonics_TEother_merge_states.txt",sep='\t',header=TRUE,row.names=1)))
mnemonics_states_TEother[is.na(mnemonics_states_TEother)] = 0

# Number of bases in each state in each sample in merged all TEs, normalized by total annotated
mnemonics_states_TEother_normalized = mnemonics_states_TEother[2:16]/mnemonics_states_TEother$Total

# Mean, median, sd of proportion of bases in each state in merged all TEs
mnemonics_states_TEother_normalized_average = cbind(colMeans(mnemonics_states_TEother_normalized),apply(mnemonics_states_TEother_normalized,2,median),apply(mnemonics_states_TEother_normalized,2,sd))
colnames(mnemonics_states_TEother_normalized_average) = c("Mean","Median","SD")

# Number of bases in each state and composite states overall, in merged all TEs, and in promoters, normalized by total
mnemonics_states_quad_normalized_long = rbind(melt(as.matrix(mnemonics_states_normalized)),melt(as.matrix(mnemonics_states_TEother_normalized)))
colnames(mnemonics_states_quad_normalized_long) = c("Sample","State","Proportion")
mnemonics_states_quad_normalized_long$Cohort = c(rep("Genome",1905),rep("TEs",1905))
mnemonics_states_quad_normalized_long = rbind(mnemonics_states_quad_normalized_long,mnemonics_states_promoter_normalized[,c(2,1,3:4)])
mnemonics_states_quad_normalized_long$Cohort = factor(mnemonics_states_quad_normalized_long$Cohort,levels=c("Genome","TEs","Promoter"))

# Adding composite states
test = rbind(aggregate(data=mnemonics_states_quad_normalized_long[which(mnemonics_states_quad_normalized_long$State %in% chromHMM_states[1:8]),],Proportion~Sample+Cohort,sum),aggregate(data=mnemonics_states_quad_normalized_long[which(mnemonics_states_quad_normalized_long$State %in% chromHMM_states[9:15]),],Proportion~Sample+Cohort,sum),aggregate(data=mnemonics_states_quad_normalized_long[which(mnemonics_states_quad_normalized_long$State %in% chromHMM_states[c(1:3,6:7)]),],Proportion~Sample+Cohort,sum),aggregate(data=mnemonics_states_quad_normalized_long[which(mnemonics_states_quad_normalized_long$State %in% chromHMM_states[4:5]),],Proportion~Sample+Cohort,sum),aggregate(data=mnemonics_states_quad_normalized_long[which(mnemonics_states_quad_normalized_long$State %in% chromHMM_states[10:12]),],Proportion~Sample+Cohort,sum),aggregate(data=mnemonics_states_quad_normalized_long[which(mnemonics_states_quad_normalized_long$State %in% chromHMM_states[c(9,13:14)]),],Proportion~Sample+Cohort,sum))
test$State = c(rep("Active",381),rep("Inactive",381),rep("Active regulatory",381),rep("Transcribed",381),rep("Poised regulatory",381),rep("Repressed",381))
mnemonics_states_quad_normalized_long = rbind(mnemonics_states_quad_normalized_long,test)

# Wilcox test with Bonferroni correction on state proportions in genome, merged TEs
p.adjust(by(mnemonics_states_quad_normalized_long[which(mnemonics_states_quad_normalized_long$Cohort %in% c("Genome","TEs")),],mnemonics_states_quad_normalized_long[which(mnemonics_states_quad_normalized_long$Cohort %in% c("Genome","TEs")),]$State,function(x) unlist(wilcox.test(x$Proportion~x$Cohort))["p.value"]),method="bonf")

# Wilcox test with Bonferroni correction on state proportions in genome, promoters
p.adjust(by(mnemonics_states_quad_normalized_long[which(mnemonics_states_quad_normalized_long$Cohort %in% c("Genome","Promoter")),],mnemonics_states_quad_normalized_long[which(mnemonics_states_quad_normalized_long$Cohort %in% c("Genome","Promoter")),]$State,function(x) unlist(wilcox.test(x$Proportion~x$Cohort))["p.value"]),method="bonf")

# Refseq features
# Number of bases in each state in each sample in merged genic features
mnemonics_states_features = read.table("chromHMM_features/chromHMM_feature_states.txt",sep='\t')
colnames(mnemonics_states_features) = c("State","Sample","Bases","Cohort")
test = t(matrix(c("3_TxFlnk","E002",0,"3UTR","3_TxFlnk","E002",0,"intergenic","3_TxFlnk","E002",0,"repeats","10_TssBiv","E098",0,"3UTR","13_ReprPC","E098",0,"3UTR"),nrow=4)) #Adding those with no bases
colnames(test) = c("State","Sample","Bases","Cohort")
test = as.data.frame(test)
test$Bases = as.numeric(test$Bases)
mnemonics_states_features = rbind(mnemonics_states_features,test)
mnemonics_states_features[15236:15240,]$Bases = 0
test = aggregate(data=mnemonics_states_features,Bases~Sample+Cohort,sum)
mnemonics_states_features$Proportion = apply(mnemonics_states_features,1,function(x) as.numeric(x[3])/test[which(test$Sample == x[2] & test$Cohort == x[4]),]$Bases)

# Adding TEs and Genome
test = mnemonics_states_quad_normalized_long[which(mnemonics_states_quad_normalized_long$State %in% chromHMM_states & mnemonics_states_quad_normalized_long$Cohort %in% c("Genome","TEs")),]
test$Bases = rep("NA",3810)
mnemonics_states_features = rbind(mnemonics_states_features,test)
mnemonics_states_features$Cohort = factor(mnemonics_states_features$Cohort,levels=levels(mnemonics_states_features$Cohort)[c(9,10,7,2,3,1,4,6,5,8)])

# Number of bases in each state in each sample in merged genic features, no TEs
mnemonics_states_features_noTE = read.table("chromHMM/Refseq_features/chromHMM_feature_noTE_states.txt",sep='\t')
colnames(mnemonics_states_features_noTE) = c("State","Sample","Bases","Cohort")
mnemonics_expand = expand.grid(Sample = EID_metadata$Sample,State = chromHMM_states, Cohort = unique(mnemonics_states_features_noTE$Cohort))
mnemonics_states_features_noTE = merge(mnemonics_states_features_noTE,mnemonics_expand,all.y=TRUE,all.x=TRUE,by=c("Sample","State","Cohort"))
mnemonics_states_features_noTE[which(is.na(mnemonics_states_features_noTE$Bases)),]$Bases = 0
test = aggregate(data=mnemonics_states_features_noTE,Bases~Sample+Cohort,sum)
mnemonics_states_features_noTE$Proportion = apply(mnemonics_states_features_noTE,1,function(x) as.numeric(x[4])/test[which(test$Sample == x[1] & test$Cohort == x[3]),]$Bases)
