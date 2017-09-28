# H3K27ac subfamily enrichment
# See 7/4/2017, 7/24/2017

source("R_scripts/TE_subfamily_stats.R")
source("R_scripts/H3K27ac_overlap.R")
load("R_datasets/TE_H3K27ac_peaks.RData")

# Length of H3K27ac peak overlapping subfamily per sample
subfamily_H3K27ac_sample = read.table("H3K27ac/subfamily_H3K27ac_sample.txt",sep='\t')
colnames(subfamily_H3K27ac_sample) = c("subfamily","Sample","Length_ijk")

# Length of subfamily
subfamily_H3K27ac_sample_expand = expand.grid(subfamily = levels(rmsk_TE_subfamily$subfamily),Sample = levels(subfamily_H3K27ac_sample$Sample))
subfamily_H3K27ac_sample_expand = join(subfamily_H3K27ac_sample_expand,rmsk_TE_subfamily[,c(1:3,33)],by=c("subfamily"))
subfamily_H3K27ac_sample = join(subfamily_H3K27ac_sample_expand,subfamily_H3K27ac_sample,by=c("subfamily","Sample"),type="left")
colnames(subfamily_H3K27ac_sample)[5] = "Length_ik"
subfamily_H3K27ac_sample[is.na(subfamily_H3K27ac_sample$Length_ijk),]$Length_ijk = 0
rm(subfamily_H3K27ac_sample_expand)

# Total length of H3K27ac peaks per sample
subfamily_H3K27ac_sample = merge(subfamily_H3K27ac_sample,H3K27ac_stats[,c(1,4)],by=c("Sample"),all.x=TRUE)
colnames(subfamily_H3K27ac_sample)[7] = "Length_jk"

# Total length of sample
mnemonics_states_genome = melt(as.matrix(t(read.table("chromHMM/genome/mnemonics_state.txt",sep='\t',header=TRUE,row.names=1)[1,])))[,c(1,3)]
colnames(mnemonics_states_genome) = c("Sample","Length_k")
subfamily_H3K27ac_sample = join(subfamily_H3K27ac_sample,mnemonics_states_genome,by=c("Sample"),type="left")
rm(mnemonics_states_genome)

# LOR enrichment
subfamily_H3K27ac_sample$Enrichment = log2((subfamily_H3K27ac_sample$Length_ijk/subfamily_H3K27ac_sample$Length_ik)/(subfamily_H3K27ac_sample$Length_jk/subfamily_H3K27ac_sample$Length_k))

# Total proportion of H3K27ac peaks in subfamily
subfamily_H3K27ac_sample$Length_percent_jk = subfamily_H3K27ac_sample$Length_ijk/subfamily_H3K27ac_sample$Length_jk

# Add metadata
subfamily_H3K27ac_sample = merge(subfamily_H3K27ac_sample,metadata[,c(1,4:9)],by=c("Sample"),all.x=TRUE)

# Number of members overlapping H3K27ac peak (needs matrix)
TE_H3K27ac_peaks_members = merge(melt(aggregate(data=TE_H3K27ac_peaks[,c(4,8:105)],.~subfamily,function(x) sum(as.numeric(x))),id.vars=c("subfamily")),melt(aggregate(data=TE_H3K27ac_peaks[,c(4,8:105)],.~subfamily,function(x) sum(as.numeric(x) > 0)),id.vars=c("subfamily")),by=c("subfamily","variable"))
colnames(TE_H3K27ac_peaks_members)[2:4] = c("Sample","H3K27ac_peaks","TEs_overlapping")
TE_H3K27ac_peaks_members = merge(TE_H3K27ac_peaks_members,rmsk_TE_subfamily[,c(1,4)],by="subfamily")
TE_H3K27ac_peaks_members[which(TE_H3K27ac_peaks_members$Sample %in% c("E116","E117","E123","E124","E126","E127")),]$Count = rmsk_TE_subfamily[match(TE_H3K27ac_peaks_members[which(TE_H3K27ac_peaks_members$Sample %in% c("E116","E117","E123","E124","E126","E127")),]$subfamily,rmsk_TE_subfamily$subfamily),]$Count_noY

# Proportion of members overlapping H3K27ac peak
TE_H3K27ac_peaks_members$Percent = TE_H3K27ac_peaks_members$TEs_overlapping/TE_H3K27ac_peaks_members$Count

# Added to enrichment matrix
subfamily_H3K27ac_sample = merge(subfamily_H3K27ac_sample,TE_H3K27ac_peaks_members,by=c("subfamily","Sample"),all.x=TRUE)
subfamily_H3K27ac_sample[which(is.na(subfamily_H3K27ac_sample$Percent)),17:20] = 0
subfamily_H3K27ac_sample[which(subfamily_H3K27ac_sample$subfamily == "Tigger2a_Car"),]$Count = 2
subfamily_H3K27ac_sample$State = rep("H3K27ac",dim(subfamily_H3K27ac_sample)[1])
subfamily_H3K27ac_sample = subfamily_H3K27ac_sample[,c(1,3:4,21,2,6,5,7:20)]

# Number of samples subfamily is enriched
subfamily_H3K27ac_sample_counts = ddply(subfamily_H3K27ac_sample,.(class_update,family,subfamily,State),function(x) sum(x$Enrichment > 1.5 & x$Length_ijk >= 600 & x$Length_ik > 5000))

# Number of >1% per subfamily x state
subfamily_H3K27ac_sample_counts_pc = ddply(subfamily_H3K27ac_sample,.(class_update,family,subfamily,State),function(x) sum(x$Length_percent_jk > 0.01 & x$Length_ijk >= 600))
