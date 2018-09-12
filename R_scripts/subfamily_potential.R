# Potential at the subfamily level

subfamily_state_potential = ddply(subfamily_state_sample_combined,.(subfamily,State),summarise,Samples=sum(Length_ijk > 0))
subfamily_state_potential$State = factor(subfamily_state_potential$State,levels=states[1:21])
subfamily_state_potential$Sample.Proportion = ifelse(subfamily_state_potential$State %in% chromHMM_states,subfamily_state_potential$Samples/sample_counts["All","chromHMM"],
       ifelse(subfamily_state_potential$State %in% meth_states,subfamily_state_potential$Samples/sample_counts["All","WGBS"],
              ifelse(subfamily_state_potential$State == "DNase",subfamily_state_potential$Samples/sample_counts["All","DNase"],
                     ifelse(subfamily_state_potential$State == "H3K27ac",subfamily_state_potential$Samples/sample_counts["All","H3K27ac"],"NA"))))
subfamily_state_potential$Sample.Proportion = as.numeric(subfamily_state_potential$Sample.Proportion)

subfam_states = merge(ddply(subfamily_state_potential,.(subfamily),summarise,States=length(unique(State[which(Samples > 0)]))),
                      ddply(ddply(subfamily_state_sample_combined,.(subfamily,Sample),summarise,Intra=length(unique(State[which(Length_ijk > 0)]))),
                            .(subfamily),function(x) x[which.max(x$Intra),c("subfamily","Intra")]),
                      by="subfamily",all=TRUE)
subfam_states = melt(subfam_states,id.var="subfamily")
subfam_states = merge(subfam_states,rmsk_TE_subfamily[,c("subfamily","family","class_update","Total_length","Count","CpGs")],by="subfamily")
