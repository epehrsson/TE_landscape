# Contribution (TEs, Refseq)
# See 4/26/2016, 4/27/2016, 5/2/2016, 5/18/2016, 5/19/2016, 6/27/2016, 2/3/2017, 2/9/2017, 3/2/2017, 5/12/2017, 8/7/2017

# ChromHMM contribution 
contribution = read.table(file="chromHMM/TEother_contribution.txt",header=TRUE,sep='\t',row.names=1)
colnames(contribution) = c("Total",chromHMM_states)

# Refseq genic features
# Contribution to each chromHMM state
# feature_contribution = dcast(aggregate(data=mnemonics_states_features_noTE,Bases~Cohort+State,sum),Cohort~State)[,c(1:2,9:16,3:8)]
# feature_contribution$Total = apply(feature_contribution,1,function(x) sum(as.numeric(x[2:16])))
# feature_contribution = as.data.frame(cbind(feature_contribution,rowSums(feature_contribution[,2:9]),rowSums(feature_contribution[,10:16]),rowSums(feature_contribution[,c(2:4,7:8)]),rowSums(feature_contribution[,5:6]),rowSums(feature_contribution[,c(10,14:15)]),rowSums(feature_contribution[,11:13])))
# colnames(feature_contribution) = c("Cohort",chromHMM_states,"Total","Active","Inactive","Active Regulatory","Transcribed","Repressed","Poised Regulatory")

# feature_contribution_long = melt(feature_contribution,id.vars = "Cohort")
# colnames(feature_contribution_long)[2:3] = c("State","Bases")
# feature_contribution_long$Proportion = apply(feature_contribution_long,1,function(x) as.numeric(x[3])/contribution[1,x[2]])
# feature_contribution_long = feature_contribution_long[order(feature_contribution_long$Cohort,feature_contribution_long$Proportion),]

# Contribution to CpGs in each state, across all samples
# feature_CpG_meth_contribution = aggregate(data=feature_CpG_meth[,1:5],.~Cohort,sum)
# rownames(feature_CpG_meth_contribution) =  feature_CpG_meth_contribution$Cohort
# feature_CpG_meth_contribution =  feature_CpG_meth_contribution[,2:5]
# feature_CpG_meth_contribution = t(apply(feature_CpG_meth_contribution,1,function(x) x/colSums(all_CpG_meth)))
