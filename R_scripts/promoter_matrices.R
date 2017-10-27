# Creates matrices of promoter x number of samples in each state for all 4 measures

library(reshape2)

# Number of samples each promoter is in each chromHMM state
chromHMM_promoter_state = read.table("chromHMM/Refseq_promoters/potential_refseq_promoters_unique.txt",sep='\t',header=TRUE)

# Methylation level for all promoters
promoter_meth_average = read.table("WGBS/Refseq_promoters/refseq_promoter_unique_CpG_Meth_average.txt",sep='\t',header=TRUE)

# Number of samples promoter is in each methylation state (IMR90, no IMR90)
promoter_meth_average$Hypomethylated = apply(promoter_meth_average,1,function(x) sum(x[5:41] < 0.3,na.rm=TRUE))
promoter_meth_average$Hypermethylated = apply(promoter_meth_average,1,function(x) sum(x[5:41] > 0.7,na.rm=TRUE))
promoter_meth_average$Intermediate = apply(promoter_meth_average,1,function(x) sum(x[5:41] <= 0.7 & x[5:41] >= 0.3,na.rm=TRUE))
promoter_meth_average$Missing = apply(promoter_meth_average,1,function(x) sum(is.na(x[5:41])))
promoter_meth_average$Hypomethylated_noIMR90 = apply(promoter_meth_average,1,function(x) sum(x[c(5:14,16:41)] < 0.3,na.rm=TRUE))
promoter_meth_average$Hypermethylated_noIMR90 = apply(promoter_meth_average,1,function(x) sum(x[c(5:14,16:41)] > 0.7,na.rm=TRUE))
promoter_meth_average$Intermediate_noIMR90 = apply(promoter_meth_average,1,function(x) sum(x[c(5:14,16:41)] <= 0.7 & x[c(5:14,16:41)] >= 0.3,na.rm=TRUE))
promoter_meth_average$Missing_noIMR90 = apply(promoter_meth_average,1,function(x) sum(is.na(x[c(5:14,16:41)])))

# Number of overlaps between promoters, DNase peaks
promoter_DNase_peaks = read.table(file="DNase/Refseq_promoters/refseq_promoters_unique_DNase_peaks.txt",sep='\t')
colnames(promoter_DNase_peaks) = c("chromosome","start","stop","strand","Sample","Peaks","Overlap")
promoter_DNase_peaks = dcast(promoter_DNase_peaks,chromosome+start+stop+strand~Sample,value.var="Peaks")
promoter_DNase_peaks[is.na(promoter_DNase_peaks)] = 0

# Number of samples a promoter overlaps a DNase peak
promoter_DNase_peaks$Samples = apply(promoter_DNase_peaks,1,function(x) sum(as.numeric(x[5:57]) > 0))

# Number of overlaps between promoters, H3K27ac peaks
promoter_H3K27ac_peaks = read.table(file="H3K27ac/Refseq_promoters/refseq_promoters_unique_H3K27ac_peaks.txt",sep='\t')
colnames(promoter_H3K27ac_peaks) = c("chromosome","start","stop","strand","Sample","Peaks","Overlap")
promoter_H3K27ac_peaks = dcast(promoter_H3K27ac_peaks,chromosome+start+stop+strand~Sample,value.var="Peaks")
promoter_H3K27ac_peaks[is.na(promoter_H3K27ac_peaks)] = 0

# Number of samples a promoter overlaps a H3K27ac peak
promoter_H3K27ac_peaks$Samples = apply(promoter_H3K27ac_peaks,1,function(x) sum(as.numeric(x[5:102]) > 0))

save(list=c("chromHMM_promoter_state","promoter_meth_average","promoter_DNase_peaks","promoter_H3K27ac_peaks"),file="R_datasets/promoter_matrices.RData")