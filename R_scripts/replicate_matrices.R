# Load TE x sample x state
print("Load state matrices")

## chromHMM
chromHMM_active_matrix = lapply(states[c(1:3,6:7)],function(x) read.table(paste("chromHMM/chromHMM_",x,".txt",sep=""),sep='\t',
                                                                          col.names=c(TE_coordinates[c(1:4,6,5,7)],"Sample","Overlap","State","Category"))[,c(1:8)])
names(chromHMM_active_matrix) = states[c(1:3,6:7)]
chromHMM_active_matrix = ldply(chromHMM_active_matrix)
colnames(chromHMM_active_matrix)[1] = "State"

# Load whole genome
windows_active = read.table("chromHMM/genome/windows/windows_active_reg.bed",sep='\t')
colnames(windows_active) = c("chr","start","stop","State","Sample")