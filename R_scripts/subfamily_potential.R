# Potential at the subfamily level

subfamily_state_potential = ddply(subfamily_state_sample_combined,.(subfamily,State),summarise,Samples=sum(Length_ijk > 0))
subfamily_state_potential$State = factor(subfamily_state_potential$State,levels=states[1:21])

subfam_states = merge(ddply(subfamily_state_potential,.(subfamily),summarise,States=length(unique(State[which(Samples > 0)]))),
                      ddply(ddply(subfamily_state_sample_combined,.(subfamily,Sample),summarise,Intra=length(unique(State[which(Length_ijk > 0)]))),
                            .(subfamily),function(x) x[which.max(x$Intra),c("subfamily","Intra")]),
                      by="subfamily",all=TRUE)
subfam_states = melt(subfam_states,id.var="subfamily")