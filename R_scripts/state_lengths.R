# State lengths/peak widths
# See 5/16/2016, 7/31/2017, 8/3/2017
# Update 5/2/18 by sample

# Lengths of all chromHMM blocks by sample and state
chromHMM_lengths = read.table("chromHMM/chromHMM_blocks.txt",sep='\t')
colnames(chromHMM_lengths) = c("Length","State","Sample")
chromHMM_lengths$State = factor(chromHMM_lengths$State,levels=chromHMM_states)

# Lengths of DNase peaks by sample
DNase_peaks = read.table(file="DNase/peak_widths.txt",sep='\t')
colnames(DNase_peaks) = c("Length","Sample")
DNase_peaks$State = rep("DNase",dim(DNase_peaks)[1])

# Lengths of H3K27ac peaks by sample
H3K27ac_peaks = read.table(file="H3K27ac/peak_widths.txt",sep='\t')
colnames(H3K27ac_peaks) = c("Length","Sample")
H3K27ac_peaks$State = rep("H3K27ac",dim(H3K27ac_peaks)[1])

all_lengths = rbind(chromHMM_lengths,DNase_peaks,H3K27ac_peaks)
rm(chromHMM_lengths)
rm(DNase_peaks)
rm(H3K27ac_peaks)