# Load TE x sample x state
print("Load state matrices")

## chromHMM
chromHMM_active_matrix = lapply(states[c(1:3,6:7)],function(x) read.table(paste("chromHMM/chromHMM_",x,".txt",sep=""),sep='\t',
                                                                          col.names=c(TE_coordinates[c(1:4,6,5,7)],"Sample","Overlap","State","Category"))[,c(1:8)])
names(chromHMM_active_matrix) = states[c(1:3,6:7)]
chromHMM_active_matrix = ldply(chromHMM_active_matrix)
colnames(chromHMM_active_matrix)[1] = "State"

## H3K27ac
H3K27ac_pairs = melt(TE_H3K27ac_peaks[,c(TE_coordinates,as.vector(metadata[which(!is.na(metadata$H3K27ac)),]$Sample),"Samples")],
                     id.vars=c(TE_coordinates,"Samples"))
colnames(H3K27ac_pairs)[8:10] = c("Total","Sample","Peaks")
H3K27ac_pairs = H3K27ac_pairs[which(H3K27ac_pairs$Peaks > 0),c(TE_coordinates,"Total","Sample")]
H3K27ac_pairs$State = rep("H3K27ac",dim(H3K27ac_pairs)[1])

# Load whole genome
windows_active = read.table("chromHMM/genome/windows/windows_active_reg.bed",sep='\t')
colnames(windows_active) = c("chr","start","stop","State","Sample")