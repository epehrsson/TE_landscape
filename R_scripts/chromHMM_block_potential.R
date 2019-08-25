# Creates a combined dataframe with the potential for TEs to be in each chromHMM state
# Using the standard rules, requiring overlap with the center of a 200bp bin, 
# Or requiring overlap with the center of a chromHMM block

## block_potential: Number of samples each TE is in each chromHMM state, overlapping chromHMM block center
## summit_potential: Number of samples each TE is in each chromHMM state, overlapping chromHMM 200bp bin center
## state_sample_count_rules: Number of TEs in state per sample, all 3 rules, plus ratio of new to original results
## potential_ever: Proportion of TEs ever in each state, all 3 rules, plus ratio of new to original results

# Requiring overlap with chromHMM block center

## Number of TEs in state per sample
print("Load block per sample")
state_sample_count_blocks = read.table("chromHMM/block/rmsk_TE_block_summit_counts.txt",sep='\t')
colnames(state_sample_count_blocks) = c("State","Sample","Count")
state_sample_count_blocks[1905,] = c("3_TxFlnk","E002",0)
state_sample_count_blocks$Count = as.numeric(state_sample_count_blocks$Count)

## Proportion of TEs ever in each state
print("Block ever")

### Number of samples each TE is in each chromHMM state
block_potential = read.table("chromHMM/block/rmsk_TE_block_potential.txt",sep = '\t')
colnames(block_potential) = c(TE_coordinates[c(1:4,6,5,7)],"State","Samples")
block_potential_ever = ddply(block_potential,.(State),summarise,Ever=length(State[which(Samples > 0)])/NUM_TE)

# Requiring overlap with chromHMM 200bp bin center

## Number of TEs in state per sample
print("Load summit per sample")
state_sample_count_summit = read.table("chromHMM/state_sample_counts_summit_only.txt",sep='\t')
colnames(state_sample_count_summit) = c("Sample","State","Count")
state_sample_count_summit[1905,] = c("E002","3_TxFlnk",0)
state_sample_count_summit$Count = as.numeric(state_sample_count_summit$Count)

## Proportion of TEs ever in each state
print("Summit ever")

### Number of samples each TE is in each chromHMM state,
### Filtered to only TEs that overlap the center of a 200bp bin center
summit_potential = chromHMM_TE_state[which(chromHMM_TE_state$Category == "summit"),]

### Number of TEs in each state at each number of samples
summit_potential_dist = sample_distribution(summit_potential,c(8:22),sample_counts["All","chromHMM"])
summit_potential_dist = melt(summit_potential_dist,id.var="Samples")
colnames(summit_potential_dist)[2:3] = c("State","Count")
summit_potential_ever = ddply(summit_potential_dist,.(State),summarise,Ever=sum(Count[which(Samples > 0)])/NUM_TE)

# Combine dataframes for all 3 rules

## Number of TEs in state per sample
print("Combine per sample")
state_sample_count_rules = merge(merge(state_sample_count[,c("State","Sample","Count")],
                   state_sample_count_blocks,by=c("State","Sample")),
             state_sample_count_summit,by=c("State","Sample"))
colnames(state_sample_count_rules)[3:5] = c("Count","Block","Summit") 
state_sample_count_rules = melt(state_sample_count_rules[,c("State","Sample","Count","Block","Summit")],id.vars=c("State","Sample","Count"))
colnames(state_sample_count_rules)[4:5] = c("Category","Count.Filter")

### Divide number of TEs in state per sample using new rule by number using standard rule
state_sample_count_rules$Ratio = state_sample_count_rules$Count.Filter/state_sample_count_rules$Count

rm(list=c("state_sample_count_blocks","state_sample_count_summit"))

## Proportion of TEs ever in each state
print("Combine ever")
potential_ever = merge(merge(chromHMM_TE_state_dist_stats[,c("Proportion_ever","State")],block_potential_ever,by="State"),
                             summit_potential_ever,by="State")
colnames(potential_ever)[3:4] = c("Block","Summit")
potential_ever = melt(potential_ever,id.vars=c("State","Proportion_ever"))
colnames(potential_ever)[3:4] = c("Category","Proportion_ever.Filter")

### Divide proportion of TEs ever in state using new rule by number using standard rule
potential_ever$Ratio = potential_ever$Proportion_ever.Filter/potential_ever$Proportion_ever

rm(list=c("block_potential_ever","summit_potential_dist","summit_potential_ever"))