# 50-state chromHMM models

# Number of TEs in each chromHMM state by sample
state50_files = list.files(path="/scratch/ecp/TE_landscape/state50",pattern="state_sample_counts_rmsk_TEother_50state_",full.names = TRUE)
state_sample_50state = lapply(state50_files,function(x) read.table(x,sep='\t',quote=""))
names(state_sample_50state) = gsub(".txt","",gsub("/scratch/ecp/TE_landscape/state50/state_sample_counts_rmsk_TEother_50state_","",state50_files))
state_sample_50state = ldply(state_sample_50state)
colnames(state_sample_50state) = c("Sample","State50","Count")

# Corresponding 18-state model state
state50_state18 = read.table("sample_lists/chromHMM_50state.txt",sep='\t',quote="",header=TRUE)
state50_state18$State50 = paste("E",state50_state18$State50,sep="")

# Combine
state_sample_50state = merge(state_sample_50state,state50_state18,by=c("Sample","State50"))
state_sample_50state$Proportion = state_sample_50state$Count/NUM_TE
state_sample_50state$State50 = factor(state_sample_50state$State50,levels=paste("E",seq(1,50,1),sep=""))

# Plot
ggplot(state_sample_50state,aes(x=State50,y=Proportion,fill=State18)) + geom_bar(stat="identity",color="black") + scale_fill_manual(values=chromHMM_states_18,guide=FALSE) + coord_flip() +
  facet_wrap(~Sample) + labs(x="State with 50-state model",y="Proportion of TEs in state") + scale_x_discrete(limits=rev(levels(state_sample_50state$State50)))