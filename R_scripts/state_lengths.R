# State lengths/peak widths
# See 5/16/2016, 7/31/2017, 8/3/2017

# Lengths of all chromHMM blocks in all samples
chromHMM_lengths = read.table("chromHMM/chromHMM_blocks.txt",sep=' ')
colnames(chromHMM_lengths) = c("Length","State")
chromHMM_lengths$State = factor(chromHMM_lengths$State,levels=chromHMM_states)

# Lengths of DNase peaks by sample
DNase_peaks = read.table(file="DNase/peak_widths.txt",sep='\t')
colnames(DNase_peaks) = c("Length","Sample")
DNase_peaks$State = rep("DNase",dim(DNase_peaks)[1])

# Lengths of H3K27ac peaks by sample
H3K27ac_peaks = read.table(file="raw_data/H3K27ac/H3K27ac_narrow_peaks/peak_widths.txt",sep='\t')
colnames(H3K27ac_peaks) = c("Length","Sample")
H3K27ac_peaks$State = rep("H3K27ac",dim(H3K27ac_peaks)[1])

all_lengths = rbind(chromHMM_lengths,DNase_peaks[,c(1,3)],H3K27ac_peaks[,c(1,3)])
rm(chromHMM_lengths)
rm(DNase_peaks)
rm(H3K27ac_peaks)