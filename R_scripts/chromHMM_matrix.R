# Loads matrix of TE x number of samples in each state (all samples, no cancer cell lines/IMR90)
# See 4/25/2016, 4/26/2016, 4/27/2016, 5/3/2016, 5/11/2016, 5/12/2016, 6/27/2016, 8/30/2016, 9/20/2016, 9/28/2016, 2/3/2017, 2/6/2017, 5/10/2017, 6/14/2017, 8/1/2017, 8/2/2017

# Number of samples each TE is in each chromHMM state
chromHMM_TE_state = rbind(read.table("chromHMM/potential/all_chromHMM_TE_potential_0_nodiv.txt",sep='\t',header=TRUE),read.table("chromHMM/potential/all_chromHMM_other_potential.txt",sep='\t',header=TRUE))
colnames(chromHMM_TE_state)[8:22] = chromHMM_states
chromHMM_TE_state$States = apply(chromHMM_TE_state[,8:22],1,function(x) sum(x > 0))
chromHMM_TE_state$Tissues = apply(chromHMM_TE_state[,8:22],1,sum)

# Maximum number of states in a single sample per TE
chromHMM_max_intra = rbind(read.table("chromHMM/all_chromHMM_TE_max.txt",sep='\t'),read.table("chromHMM/all_chromHMM_other_max.txt",sep='\t'))
colnames(chromHMM_max_intra) = c(TE_coordinates[c(1:4,6,5,7)],"Max_states_intra")
chromHMM_TE_state = merge(chromHMM_TE_state,chromHMM_max_intra,by=TE_coordinates)
rm(chromHMM_max_intra)

# Updating classes
chromHMM_TE_state$class_update = convert_class(chromHMM_TE_state$class)

# Number of samples each TE is in each chromHMM state, no cancer cell lines or IMR90
chromHMM_TE_state_noCancer = rbind(read.table("chromHMM/potential/all_chromHMM_TE_noCancer_potential_0.txt",sep='\t',header=TRUE),read.table("chromHMM/potential/all_chromHMM_other_potential_noCancer.txt",sep='\t',header=TRUE))
colnames(chromHMM_TE_state_noCancer)[8:22] = chromHMM_states

save(list=c("chromHMM_TE_state","chromHMM_TE_state_noCancer"),file="R_datasets/chromHMM_TE_state.RData")
