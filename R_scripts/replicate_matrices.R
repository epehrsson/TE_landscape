# Load TE x sample x state
# Includes all TEs in the state in each sample, plus total samples in state per TE
print("Load state matrices")

## chromHMM
chromHMM_active_matrix = lapply(states[c(1:3,6:7)],function(x) load_state(x))
names(chromHMM_active_matrix) = states[c(1:3,6:7)]
chromHMM_active_matrix = ldply(chromHMM_active_matrix)
colnames(chromHMM_active_matrix)[1] = "State"

## H3K27ac
H3K27ac_pairs = melt(TE_H3K27ac_peaks[,c(TE_coordinates,as.vector(metadata[which(!is.na(metadata$H3K27ac)),]$Sample),"Samples")],
                     id.vars=c(TE_coordinates,"Samples"))
colnames(H3K27ac_pairs)[8:10] = c("Total","Sample","Peaks")
H3K27ac_pairs = H3K27ac_pairs[which(H3K27ac_pairs$Peaks > 0),c(TE_coordinates,"Total","Sample")]
H3K27ac_pairs$State = rep("H3K27ac",dim(H3K27ac_pairs)[1])

chromHMM_active_matrix = rbind(chromHMM_active_matrix,H3K27ac_pairs)

# Generate matrices of TE x samples in Group by state
print("Add metadata")
chromHMM_active_group = merge(chromHMM_active_matrix,metadata[,c("Sample","Group")],by="Sample")

print("Count samples in group")
chromHMM_active_group = ddply(chromHMM_active_group,.(State,chromosome,start,stop,subfamily,family,class,strand,Group),summarise,Group.Count=length(Group))