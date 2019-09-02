# Creates a dataframe ("all_lengths") of the width of all chromHMM annotation blocks, 
# DHS peaks, and H3K27ac peaks, by sample and state

# Load lengths of all chromHMM blocks by sample and state
chromHMM_lengths = read.table("chromHMM/chromHMM_blocks.txt",sep='\t',
                              col.names=c("Length","State","Sample"))
chromHMM_lengths$State = factor(chromHMM_lengths$State,levels=chromHMM_states)

# Load lengths of all DHS peaks by sample
DNase_peaks = read.table(file="DNase/peak_widths.txt",sep='\t',
                         col.names=c("Length","Sample"))
DNase_peaks$State = rep("DNase",dim(DNase_peaks)[1])

# Load lengths of all H3K27ac peaks by sample
H3K27ac_peaks = read.table(file="H3K27ac/peak_widths.txt",sep='\t',
                           col.names=c("Length","Sample"))
H3K27ac_peaks$State = rep("H3K27ac",dim(H3K27ac_peaks)[1])

# Combine into one dataframe
all_lengths = rbind(chromHMM_lengths,DNase_peaks,H3K27ac_peaks)

rm(list=c("chromHMM_lengths","DNase_peaks","H3K27ac_peaks"))