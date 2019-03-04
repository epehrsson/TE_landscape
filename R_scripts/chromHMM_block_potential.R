# Block
## Per sample
state_sample_count_blocks = read.table("chromHMM/block/rmsk_TE_block_summit_counts.txt",sep='\t')
colnames(state_sample_count_blocks) = c("State","Sample","Count")
state_sample_count_blocks[1905,] = c("3_TxFlnk","E002",0)
state_sample_count_blocks$Count = as.numeric(state_sample_count_blocks$Count)

## Ever
block_potential = read.table("chromHMM/block/rmsk_TE_block_potential.txt",sep = '\t')
colnames(block_potential) = c(TE_coordinates[c(1:4,6,5,7)],"State","Samples")
block_potential_ever = ddply(block_potential,.(State),summarise,Ever=length(State[which(Samples > 0)])/NUM_TE)

### Any active regulatory state
dim(unique(block_potential[which(block_potential$State %in% chromHMM_states[c(1:3,6:7)]),TE_coordinates]))[1]/NUM_TE

# Filtered to summit only
## By-sample
state_sample_count_summit = read.table("chromHMM/state_sample_counts_summit_only.txt",sep='\t')
colnames(state_sample_count_summit) = c("Sample","State","Count")
state_sample_count_summit[1905,] = c("E002","3_TxFlnk",0)
state_sample_count_summit$Count = as.numeric(state_sample_count_summit$Count)

## Ever
summit_potential = chromHMM_TE_state[which(chromHMM_TE_state$Category == "summit"),]
summit_potential_dist = sample_distribution(summit_potential,c(8:22),sample_counts["All","chromHMM"])
summit_potential_dist = melt(summit_potential_dist,id.var="Samples")
colnames(summit_potential_dist)[2:3] = c("State","Count")
summit_potential_ever = ddply(summit_potential_dist,.(State),summarise,Ever=sum(Count[which(Samples > 0)])/NUM_TE)

### Any active regulatory state
sum(apply(summit_potential[c(8:10,13:14)],1,function(x) sum(x > 0)) > 0)/dim(summit_potential)[1]

# Combine
## By sample
test = merge(merge(state_sample_count[,c("State","Sample","Count")],
                   state_sample_count_blocks,by=c("State","Sample")),
             state_sample_count_summit,by=c("State","Sample"))
colnames(test)[3:5] = c("Count","Block","Summit") 
test = melt(test[,c("State","Sample","Count","Block","Summit")],id.vars=c("State","Sample","Count"))
colnames(test)[4:5] = c("Category","Count.Filter")
test$Ratio = test$Count.Filter/test$Count

test2 = ddply(test,.(State,Category),summarise,Median=median(na.omit(Ratio)))
test2[order(test2$Median),]

## Ever
potential_ever = merge(merge(chromHMM_TE_state_dist_stats[,c("Proportion_ever","State")],block_potential_ever,by="State"),
                             summit_potential_ever,by="State")
colnames(potential_ever)[3:4] = c("Block","Summit")
potential_ever = melt(potential_ever,id.vars=c("State","Proportion_ever"))
colnames(potential_ever)[3:4] = c("Category","Proportion_ever.Filter")
potential_ever$Ratio = potential_ever$Proportion_ever.Filter/potential_ever$Proportion_ever

# Plot (x-axis in increasing order of median block size)
ggplot(test,aes(x=State,y=Ratio,fill=Category)) + geom_boxplot() + ylim(0,1) + 
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) + 
  ylab("Proportion of paper TEs") + scale_x_discrete(limits=chromHMM_states[c(2,10:12,3,6:8,1,13,9,4:5,14:15)]) 

ggplot(potential_ever,aes(x=State,y=Ratio,color=Category)) + geom_point(size=4) + ylim(0,1) +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) + 
  ylab("Proportion of paper TEs") + scale_x_discrete(limits=chromHMM_states[c(2,10:12,3,6:8,1,13,9,4:5,14:15)])
