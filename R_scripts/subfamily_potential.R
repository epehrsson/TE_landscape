# Potential at the subfamily level

subfamily_state_potential = ddply(subfamily_state_sample_combined,.(subfamily,State),summarise,Samples=sum(Length_ijk > 0))
subfamily_state_potential$State = factor(subfamily_state_potential$State,levels=states[1:21])
subfamily_state_potential$Metric = ifelse(subfamily_state_potential$State %in% chromHMM_states,"chromHMM",
                                          ifelse(subfamily_state_potential$State %in% meth_states,"WGBS",as.character(subfamily_state_potential$State)))
subfamily_state_potential$Sample.Proportion = as.numeric(subfamily_state_potential$Samples/sample_counts["All",subfamily_state_potential$Metric])

subfam_intra = list(ddply(ddply(subfamily_state_sample_combined[which(subfamily_state_sample_combined$State %in% chromHMM_states),],.(subfamily,Sample),
                           summarise,Intra=length(State[which(Length_ijk > 0)])),.(subfamily),function(x) x[which.max(x$Intra),c("subfamily","Intra")]),
                    ddply(droplevels(ddply(subfamily_state_sample_combined[which(subfamily_state_sample_combined$State %in% meth_states),],.(subfamily,Sample),
                              summarise,Intra=length(State[which(Length_ijk > 0)]))),.(subfamily),function(x) x[which.max(x$Intra),c("subfamily","Intra")]),
                    ddply(ddply(subfamily_state_sample_combined[which(subfamily_state_sample_combined$State == "DNase"),],.(subfamily,Sample),
                                 summarise,Intra=length(State[which(Length_ijk > 0)])),.(subfamily),function(x) x[which.max(x$Intra),c("subfamily","Intra")]),
                    ddply(ddply(subfamily_state_sample_combined[which(subfamily_state_sample_combined$State == "H3K27ac"),],.(subfamily,Sample),
                                 summarise,Intra=length(State[which(Length_ijk > 0)])),.(subfamily),function(x) x[which.max(x$Intra),c("subfamily","Intra")]))
names(subfam_intra) = c("chromHMM","WGBS","DNase","H3K27ac")
subfam_intra = ldply(subfam_intra)
colnames(subfam_intra)[1] = "Metric"

subfam_states = merge(ddply(subfamily_state_potential,.(subfamily,Metric),summarise,States=length(State[which(Samples > 0)])),subfam_intra,by=c("subfamily","Metric"),all=TRUE)
subfam_states = melt(subfam_states,id.vars=c("subfamily","Metric"))
subfam_states = merge(subfam_states,rmsk_TE_subfamily[,c("subfamily","family","class_update","Total_length","Count","CpGs")],by="subfamily")
subfam_states$Metric = factor(subfam_states$Metric,levels=c("chromHMM","WGBS","DNase","H3K27ac"))
