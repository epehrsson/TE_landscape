# Loads data frames of TEs and 200bp windows annotated with a chromHMM active regulatory state
# For all samples (n=127)

## chromHMM_active_matrix: TEs in an active regulatory state in a sample
## windows_active: 200bp windows in an active regulatory state in a sample

# Load data frames of TEs in each active regulatory state
chromHMM_active_matrix = lapply(states[c(1:3,6:7)],function(x) read.table(paste("chromHMM/chromHMM_",x,".txt",sep=""),sep='\t',
                                                                          col.names=c(TE_coordinates[c(1:4,6,5,7)],"Sample","Overlap","State","Category"))[,c(1:8)])
names(chromHMM_active_matrix) = states[c(1:3,6:7)]

## Combine into a single dataframe
chromHMM_active_matrix = ldply(chromHMM_active_matrix)
colnames(chromHMM_active_matrix)[1] = "State"

# Load data frame of 200bp windows in an active regulatory state
windows_active = read.table("chromHMM/genome/windows/windows_active_reg.bed",sep='\t')
colnames(windows_active) = c("chr","start","stop","State","Sample")