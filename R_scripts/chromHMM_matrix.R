# Load matrices of TE x number of samples in each state (all samples and no cancer cell lines/IMR90)
# See 4/25/2016, 4/26/2016, 4/27/2016, 5/3/2016, 5/11/2016, 5/12/2016, 6/27/2016, 8/30/2016, 9/20/2016, 9/28/2016, 2/3/2017, 
# 2/6/2017, 5/10/2017, 6/14/2017, 8/1/2017, 8/2/2017
# Switched to TEs in summit only 5/1/18

# Number of samples each TE is in each chromHMM state
print("Load chromHMM")

chromHMM_TE_state = read.table("chromHMM/potential/rmsk_TEother_chromHMM_summit_potential.txt",sep='\t',header=TRUE)
colnames(chromHMM_TE_state)[8:22] = chromHMM_states
chromHMM_TE_state$States = apply(chromHMM_TE_state[,8:22],1,function(x) sum(x > 0))
chromHMM_TE_state$Tissues = apply(chromHMM_TE_state[,8:22],1,sum)

# Maximum number of states in a single sample per TE
print("Add max intra")

chromHMM_max_intra = read.table("chromHMM/rmsk_TEother_chromHMM_summit_max.txt",sep='\t')
colnames(chromHMM_max_intra) = c(TE_coordinates[c(1:4,6,5,7)],"Max_states_intra")
chromHMM_TE_state = merge(chromHMM_TE_state,chromHMM_max_intra,by=TE_coordinates,all.x=TRUE)
chromHMM_TE_state[is.na(chromHMM_TE_state)] = 0
rm(chromHMM_max_intra)

# Updating classes
print("Update class")

chromHMM_TE_state$class_update = convert_class(chromHMM_TE_state$class)

# Summit or majority assignment
print("Summit assignments")

summit = read.table("/scratch/ecp/pandas/E001",sep='\t',header=TRUE)
summit = unique(summit[,c(TE_coordinates,"Category")])
chromHMM_TE_state = merge(chromHMM_TE_state,summit,by=TE_coordinates,all.x=TRUE)
rm(summit)
  
# Number of samples each TE is in each chromHMM state, no cancer cell lines or IMR90
print("Load chromHMM noCancer")

chromHMM_TE_state_noCancer = read.table("chromHMM/potential/rmsk_TEother_chromHMM_summit_potential_noCancer.txt",sep='\t',header=TRUE)
colnames(chromHMM_TE_state_noCancer)[8:22] = chromHMM_states

save(list=c("chromHMM_TE_state","chromHMM_TE_state_noCancer"),file="R_datasets/chromHMM_TE_state.RData")
