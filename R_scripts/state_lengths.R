# State lengths/peak widths
# See 5/16/2016, 7/31/2017, 8/3/2017

# Lengths of all chromHMM blocks in all samples
chromHMM_lengths = read.table("chromHMM_blocks.txt",sep=' ')
colnames(chromHMM_lengths) = c("Length","State")
chromHMM_lengths$State = factor(chromHMM_lengths$State,levels=chromHMM_states)

# Lengths of DNase peaks by sample
DNase_peaks = read.table(file="DNase/peak_widths.txt",sep='\t')

# Lengths of H3K27ac peaks by sample
H3K27ac_peaks = read.table(file="raw_data/H3K27ac/H3K27ac_narrow_peaks/peak_widths.txt",sep='\t')
mean(H3K27ac_peaks$V1)
sd(H3K27ac_peaks$V1)
mean(aggregate(data=H3K27ac_peaks,V1~V2,sum)$V1)
sd(aggregate(data=H3K27ac_peaks,V1~V2,sum)$V1)