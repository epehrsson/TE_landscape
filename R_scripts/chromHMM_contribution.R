# Contribution to chromHMM states (TEs, Refseq)
# See 4/26/2016, 4/27/2016, 5/2/2016, 5/18/2016, 5/19/2016, 6/27/2016, 2/3/2017, 2/9/2017, 3/2/2017, 5/12/2017, 8/7/2017

# Contribution of TEs/TE classes to each chromHMM state (TE vs. non-TE, total is all bases across all tissues; classes, total is all unmerged TE bases across all tissues)
contribution = read.table(file="TEother_contribution.txt",header=TRUE,sep='\t',row.names=1)
contribution = as.data.frame(cbind(contribution,rowSums(contribution[,2:9]),rowSums(contribution[,10:16]),rowSums(contribution[,c(2:4,7:8)]),rowSums(contribution[,5:6]),rowSums(contribution[,c(10,14:15)]),rowSums(contribution[,11:13])))
colnames(contribution) = c("Total",chromHMM_states,"Active","Inactive","Active Regulatory","Transcribed","Repressed","Poised Regulatory")

# Contribution of TEs/non-TEs to each chromHMM state
contribution_TE = as.data.frame(t(contribution[1:3,]))
contribution_TE = (contribution_TE/contribution_TE$Total)[,2:3]

contribution_TE_long = melt(as.matrix(contribution_TE))
colnames(contribution_TE_long) = c("State","Class","Proportion")
contribution_TE_long = contribution_TE_long[order(contribution_TE_long$Class,contribution_TE_long$Proportion),]

# Contribution of TEs vs. non-TEs (Reordered)
order_states = contribution_TE_long$State[1:22]

# Contribution of Refseq genic features to each chromHMM state
feature_contribution = dcast(aggregate(data=mnemonics_states_features_noTE,Bases~Cohort+State,sum),Cohort~State)[,c(1:2,9:16,3:8)]
feature_contribution$Total = apply(feature_contribution,1,function(x) sum(as.numeric(x[2:16])))
feature_contribution = as.data.frame(cbind(feature_contribution,rowSums(feature_contribution[,2:9]),rowSums(feature_contribution[,10:16]),rowSums(feature_contribution[,c(2:4,7:8)]),rowSums(feature_contribution[,5:6]),rowSums(feature_contribution[,c(10,14:15)]),rowSums(feature_contribution[,11:13])))
colnames(feature_contribution) = c("Cohort",chromHMM_states,"Total","Active","Inactive","Active Regulatory","Transcribed","Repressed","Poised Regulatory")

feature_contribution_long = melt(feature_contribution,id.vars = "Cohort")
colnames(feature_contribution_long)[2:3] = c("State","Bases")
feature_contribution_long$Proportion = apply(feature_contribution_long,1,function(x) as.numeric(x[3])/contribution[1,x[2]])
feature_contribution_long = feature_contribution_long[order(feature_contribution_long$Cohort,feature_contribution_long$Proportion),]
