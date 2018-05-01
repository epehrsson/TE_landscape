# Load matrices of TE x sample x chromHMM block overlap (summit only)

## Number of samples each TE is in each chromHMM state
TE_chromHMM_potential = read.table("chromHMM/potential/rmsk_TEother_chromHMM_summit_potential.txt",sep='\t',header=TRUE)
colnames(TE_chromHMM_potential)[8:22] = chromHMM_states
TE_chromHMM_potential$States = apply(TE_chromHMM_potential[,8:22],1,function(x) sum(x > 0))
TE_chromHMM_potential$Tissues = apply(TE_chromHMM_potential[,8:22],1,sum)

print("Load chromHMM")

# Maximum number of states in a single sample per TE
TE_chromHMM_max_intra = read.table("chromHMM/rmsk_TEother_chromHMM_summit_max.txt",sep='\t')
colnames(TE_chromHMM_max_intra) = c(TE_coordinates[c(1:4,6,5,7)],"Max_states_intra")
TE_chromHMM_potential = merge(TE_chromHMM_potential,TE_chromHMM_max_intra,by=TE_coordinates,all.x=TRUE)
TE_chromHMM_potential[is.na(TE_chromHMM_potential)] = 0
rm(TE_chromHMM_max_intra)

print("Add max intra")

# Updating classes
TE_chromHMM_potential$class_update = convert_class(TE_chromHMM_potential$class)

print("Update class")

# Number of samples each TE is in each chromHMM state, no cancer cell lines or IMR90
TE_chromHMM_potential_noCancer = read.table("chromHMM/potential/rmsk_TEother_chromHMM_summit_potential_noCancer.txt",sep='\t',header=TRUE)
colnames(TE_chromHMM_potential_noCancer)[8:22] = chromHMM_states

print("Load chromHMM noCancer")

save(list=c("TE_chromHMM_potential","TE_chromHMM_potential_noCancer"),file="R_datasets/TE_chromHMM_summit.RData")