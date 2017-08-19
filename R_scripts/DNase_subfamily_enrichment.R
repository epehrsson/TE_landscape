# DNase subfamily enrichment
# See 6/5/17, 6/6/17, 6/7/17, 6/9/17

source("R_scripts/TE_subfamily_stats.R")
source("R_scripts/DNase_overlap.R")
load("R_scripts/TE_DNase_peaks.RData")

# Length of DNase peak overlapping subfamily per sample
subfamily_DNase_sample = read.table("DNase/subfamily_DNase_sample.txt",sep='\t')
colnames(subfamily_DNase_sample) = c("subfamily","Sample","Length_ijk")

# Length of subfamily
subfamily_DNase_sample_expand = expand.grid(subfamily = levels(rmsk_TE_subfamily$subfamily),Sample = levels(subfamily_DNase_sample$Sample))
subfamily_DNase_sample_expand = join(subfamily_DNase_sample_expand,rmsk_TE_subfamily[,c(1:3,21)],by=c("subfamily"))
subfamily_DNase_sample = join(subfamily_DNase_sample_expand,subfamily_DNase_sample,by=c("subfamily","Sample"),type="left")
colnames(subfamily_DNase_sample)[5] = "Length_ik"
subfamily_DNase_sample[is.na(subfamily_DNase_sample$Length_ijk),]$Length_ijk = 0
rm(subfamily_DNase_sample_expand)

# Total length of DNase peaks per sample
subfamily_DNase_sample = merge(subfamily_DNase_sample,DNase_stats[,c(1,4)],by=c("Sample"),all.x=TRUE)
colnames(subfamily_DNase_sample)[7] = "Length_jk"

# Total length of sample
mnemonics_states_genome = melt(as.matrix(t(read.table("chromHMM/genome/mnemonics_state.txt",sep='\t',header=TRUE,row.names=1)[1,])))[,c(1,3)]
colnames(mnemonics_states_genome) = c("Sample","Length_k")
subfamily_DNase_sample = join(subfamily_DNase_sample,mnemonics_states_genome,by=c("Sample"),type="left")
rm(mnemonics_states_genome)

# LOR enrichment
subfamily_DNase_sample$Enrichment = log2((subfamily_DNase_sample$Length_ijk/subfamily_DNase_sample$Length_ik)/(subfamily_DNase_sample$Length_jk/subfamily_DNase_sample$Length_k))

# Total proportion of DNase peaks in subfamily
subfamily_DNase_sample$Length_percent_jk = subfamily_DNase_sample$Length_ijk/subfamily_DNase_sample$Length_jk

# Add metadata
subfamily_DNase_sample = merge(subfamily_DNase_sample,metadata[,c(1,4:9)],by=c("Sample"),all.x=TRUE)

# Number of members overlapping DNase peak (needs matrix)
TE_DNase_peaks_members = merge(melt(aggregate(data=TE_DNase_peaks[,c(4,8:60)],.~subfamily,function(x) sum(as.numeric(x))),id.vars=c("subfamily")),melt(aggregate(data=TE_DNase_peaks[,c(4,8:60)],.~subfamily,function(x) sum(as.numeric(x) > 0)),id.vars=c("subfamily")),by=c("subfamily","variable"))
colnames(TE_DNase_peaks_members)[2:4] = c("Sample","DNase_peaks","TEs_overlapping")
TE_DNase_peaks_members = merge(TE_DNase_peaks_members,rmsk_TE_subfamily[,c(1,4)],by="subfamily")
TE_DNase_peaks_members[which(TE_DNase_peaks_members$Sample %in% c("E116","E117","E123","E124","E126","E127")),]$Count = rmsk_TE_subfamily[match(TE_DNase_peaks_members[which(TE_DNase_peaks_members$Sample %in% c("E116","E117","E123","E124","E126","E127")),]$subfamily,rmsk_TE_subfamily$subfamily),]$Count_noY

# Proportion of members overlapping DNase peak
TE_DNase_peaks_members$Percent = TE_DNase_peaks_members$TEs_overlapping/TE_DNase_peaks_members$Count

# Added to enrichment matrix
subfamily_DNase_sample = merge(subfamily_DNase_sample,TE_DNase_peaks_members,by=c("subfamily","Sample"),all.x=TRUE)
subfamily_DNase_sample[which(is.na(subfamily_DNase_sample$Percent)),17:20] = 0
subfamily_DNase_sample[which(subfamily_DNase_sample$subfamily == "Tigger2a_Car"),]$Count = 2
subfamily_DNase_sample = subfamily_DNase_sample[,c(1,3:4,2,6,5,7:20)]

# Number of samples subfamily is enriched
subfamily_DNase_sample_counts = ddply(subfamily_DNase_sample,.(class_update,family,subfamily),function(x) sum(x$Enrichment > 1.5 & x$Length_ijk >= 600 & x$Length_ik > 5000))
subfamily_DNase_sample_counts$State = rep("DNase",968)