# Enrichment core vs. periphery 
# See 8/29/2016, 9/8/2016, 1/20/2017, 2/10/2017, 3/16/2017

# Average proportion of a subfamily in each state across all samples
subfamily_state_members =merge(aggregate(data=subfamily_state_sample_members,Percent~Subfamily+State,FUN=sum),aggregate(data=subfamily_state_sample_members,Percent~Subfamily+State,FUN=max),by=c("Subfamily","State"))
colnames(subfamily_state_members)[3:4] = c("Total_Percent","Max_Percent")
subfamily_state_members$Avg_Percent = subfamily_state_members$Total_Percent/127
subfamily_state_members[which(subfamily_state_members$Total_Percent > 1),]$Total_Percent = 1
subfamily_state_members = merge(subfamily_state_members,potential_subfamily_state,by=c("Subfamily","State"))
mean(na.omit(subfamily_state_members$Members_in_state/subfamily_state_members$Max_Percent))
mean(subfamily_state_members$Members_in_state/subfamily_state_members$Total_Percent,na.rm=TRUE)

# Proportion of subfamily ever in state; on average in state
potential_subfamily_state = merge(melt(aggregate(data=potential_TE_state[,c(4,8:22)],.~subfamily,function(x) sum(x > 0)/length(x)),id.vars=c("subfamily")),melt(aggregate(data=potential_TE_state[,c(4,8:22)],.~subfamily,function(x) mean(x)),id.vars=c("subfamily")),by=c("subfamily","variable"))
colnames(potential_subfamily_state) = c("Subfamily","State","Members_in_state","TE_average_samples")
potential_subfamily_state$State = mapvalues(potential_subfamily_state$State,from=names(chromHMM_colors_X),to=chromHMM_states)

# Shannon entropy on subfamily proportion in state
library(entropy)
potential_subfamily_state_entropy = ddply(potential_TE_state[,c(4,8:22)],.(subfamily),function(x) apply(x[,2:16],2,function(y) entropy.empirical(y+0.000001)))
apply(potential_subfamily_state_entropy[,2:16],2,function(x) mean(na.omit(x)))

# Average proportion of a subfamily in each state across all samples, other TEs
subfamily_state_members_other =merge(aggregate(data=subfamily_state_sample_members_other,Percent~Subfamily+State,FUN=sum),aggregate(data=subfamily_state_sample_members_other,Percent~Subfamily+State,FUN=max),by=c("Subfamily","State"))
colnames(subfamily_state_members_other)[3:4] = c("Total_Percent","Max_Percent")
subfamily_state_members_other$Avg_Percent = subfamily_state_members_other$Total_Percent/127
subfamily_state_members_other[which(subfamily_state_members_other$Total_Percent > 1),]$Total_Percent = 1

# Proportion of subfamily ever in state; on average in state, other TEs
potential_subfamily_state_other = merge(melt(aggregate(data=potential_other_state[,c(4,8:22)],.~subfamily,function(x) sum(x > 0)/length(x)),id.vars=c("subfamily")),melt(aggregate(data=potential_other_state[,c(4,8:22)],.~subfamily,function(x) mean(x)),id.vars=c("subfamily")),by=c("subfamily","variable"))
colnames(potential_subfamily_state_other) = c("Subfamily","State","Members_in_state","TE_average_samples")
potential_subfamily_state_other$State = mapvalues(potential_subfamily_state_other$State,from=names(chromHMM_colors_X),to=chromHMM_states)
subfamily_state_members_other = merge(subfamily_state_members_other,potential_subfamily_state_other,by=c("Subfamily","State"))
mean(na.omit(rbind(subfamily_state_members,subfamily_state_members_other)$Members_in_state/rbind(subfamily_state_members,subfamily_state_members_other)$Max_Percent))
mean(rbind(subfamily_state_members,subfamily_state_members_other)$Members_in_state/rbind(subfamily_state_members,subfamily_state_members_other)$Total_Percent,na.rm=TRUE)

# Proportion of active TEs in subfamily by state
sort(by(rbind(subfamily_state_members,subfamily_state_members_other),rbind(subfamily_state_members,subfamily_state_members_other)$State,function(x) mean(x$Members_in_state/x$Max_Percent,na.rm=TRUE)))

# Shannon entropy on subfamily proportion in state, other TEs
potential_subfamily_state_entropy_other = ddply(potential_other_state[,c(4,8:22)],.(subfamily),function(x) apply(x[,2:16],2,function(y) entropy.empirical(y+0.000001)))
