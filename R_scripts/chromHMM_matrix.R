# Creates matrices of hg19 TEs (row) x number of samples the TE is in each chromHMM state (column)
## chromHMM_TE_state: All samples 
## chromHMM_TE_state_noCancer: Excluding cancer cell lines/IMR90

# Load dataframe with number of samples each TE is in each chromHMM state, all samples
print("Load chromHMM")

chromHMM_TE_state = read.table("chromHMM/potential/rmsk_TEother_chromHMM_summit_potential.txt",sep='\t',header=TRUE)
colnames(chromHMM_TE_state)[8:22] = chromHMM_states

# Number of unique states across all epigenomes per TE
chromHMM_TE_state$States = apply(chromHMM_TE_state[,8:22],1,function(x) sum(x > 0))

# Total number of states across all epigenomes per TE
chromHMM_TE_state$Tissues = apply(chromHMM_TE_state[,8:22],1,sum)

# Add maximum number of states annotating each TE in a single sample
print("Add max intra")

chromHMM_max_intra = read.table("chromHMM/rmsk_TEother_chromHMM_summit_max.txt",sep='\t')
colnames(chromHMM_max_intra) = c(TE_coordinates[c(1:4,6,5,7)],"Max_states_intra")
chromHMM_TE_state = merge(chromHMM_TE_state,chromHMM_max_intra,by=TE_coordinates,all.x=TRUE)
chromHMM_TE_state[is.na(chromHMM_TE_state)] = 0
rm(chromHMM_max_intra)

# Update class assignments
print("Update class")

chromHMM_TE_state$class_update = convert_class(chromHMM_TE_state$class)

# Add column indicating whether the TE overlaps the center of a 200bp window
print("Summit assignments")

summit = read.table("features/TEs/rmsk_TEother_summit.txt",sep='\t')
colnames(summit) = TE_coordinates[c(1:4,6,5,7)]
summit$Category = rep("summit",dim(summit)[1])
chromHMM_TE_state = merge(chromHMM_TE_state,summit,by=TE_coordinates,all.x=TRUE)
chromHMM_TE_state[which(is.na(chromHMM_TE_state$Category)),]$Category = "majority"
rm(summit)
  
# Load dataframe with number of samples each TE is in each chromHMM state, excluding cancer cell lines and IMR90
print("Load chromHMM noCancer")

chromHMM_TE_state_noCancer = read.table("chromHMM/potential/rmsk_TEother_chromHMM_summit_potential_noCancer.txt",sep='\t',header=TRUE)
colnames(chromHMM_TE_state_noCancer)[8:22] = chromHMM_states

# Save dataframes
save(list=c("chromHMM_TE_state","chromHMM_TE_state_noCancer"),file="R_datasets/chromHMM_TE_state.RData")