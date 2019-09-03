# Calculates the number of samples each TE subfamily is annotated with each epigenetic state,
# as well as the number of unique states each TE subfamily is annotated with across all samples. 
# Similar to potential for individual TEs to be annotated with each epigenetic state (Figure 2)

## subfamily_state_potential - number/proportion of samples in which each TE subfamily is annotated with each epigenetic state
## subfam_states - number of unique states each TE subfamily overlaps across all samples

# Number/proportion of samples in which each TE subfamily is annotated with each epigenetic state
## chromHMM: length of overlap >= 1bp, WGBS: number of CpGs >= 1, DHS/H3K27ac: number of peak summits >= 1
subfamily_state_potential = ddply(subfamily_state_sample_combined,.(subfamily,State),summarise,Samples=sum(Length_ijk > 0))
subfamily_state_potential$State = factor(subfamily_state_potential$State,levels=states[1:21])
subfamily_state_potential$Metric = ifelse(subfamily_state_potential$State %in% chromHMM_states,"chromHMM",
                                          ifelse(subfamily_state_potential$State %in% meth_states,"WGBS",as.character(subfamily_state_potential$State)))
subfamily_state_potential$Sample.Proportion = as.numeric(subfamily_state_potential$Samples/sample_counts["All",subfamily_state_potential$Metric])

# Number of unique states each TE subfamily overlaps across all samples, by epigenetic mark
subfam_states = ddply(subfamily_state_potential,.(subfamily,Metric),summarise,States=length(State[which(Samples > 0)]))
## Add family and class information, total length of subfamily, number of members, and number of CpGs overlapping the subfamily
subfam_states = merge(subfam_states,rmsk_TE_subfamily[,c("subfamily","family","class_update","Total_length","Count","CpGs")],by="subfamily")
subfam_states$Metric = factor(subfam_states$Metric,levels=c("chromHMM","WGBS","DNase","H3K27ac"))