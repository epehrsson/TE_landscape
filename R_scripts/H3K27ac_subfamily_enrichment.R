# H3K27ac subfamily enrichment
# See 7/4/2017, 7/24/2017

# Enrichment of H3K27ac peaks (width) in TE subfamilies
subfamily_H3K27ac_sample = read.table("H3K27ac/subfamily_H3K27ac_sample.txt",sep='\t')
colnames(subfamily_H3K27ac_sample) = c("Subfamily","Sample","Length_ijk")
subfamily_H3K27ac_sample_expand = expand.grid(Subfamily = levels(subfamily_state_sample$Subfamily),Sample = levels(subfamily_H3K27ac_sample$Sample))
subfamily_H3K27ac_sample_expand = join(subfamily_H3K27ac_sample_expand,rbind(rmsk_TE_stats_subfamily,rmsk_other_stats_subfamily)[,1:3],by=c("Subfamily"))
subfamily_H3K27ac_sample = join(subfamily_H3K27ac_sample_expand,subfamily_H3K27ac_sample,by=c("Subfamily","Sample"),type="left")
subfamily_H3K27ac_sample[is.na(subfamily_H3K27ac_sample$Length_ijk),]$Length_ijk = 0
subfamily_H3K27ac_sample = merge(subfamily_H3K27ac_sample,H3K27ac_stats[,c(1,4)],by=c("Sample"))
colnames(subfamily_H3K27ac_sample)[6] = "Length_jk"
subfamily_H3K27ac_sample$Length_k = apply(subfamily_H3K27ac_sample,1,function(x) mnemonics_states_TEother[x[1],]$Total)
subfamily_H3K27ac_sample = merge(subfamily_H3K27ac_sample,unique(subfamily_state_sample[,c(1:3,5,7:12,14)]),by=c("Subfamily","Family","Class","Sample"))
subfamily_H3K27ac_sample = subfamily_H3K27ac_sample[,c(1:4,8:13,5:6,14,7)]
subfamily_H3K27ac_sample$Enrichment = log2((subfamily_H3K27ac_sample$Length_ijk/subfamily_H3K27ac_sample$Length_ik)/(subfamily_H3K27ac_sample$Length_jk/subfamily_H3K27ac_sample$Length_k))

# Total proportion of H3K27ac peaks in subfamily
subfamily_H3K27ac_sample$Length_percent_jk = subfamily_H3K27ac_sample$Length_ijk/subfamily_H3K27ac_sample$Length_jk

# Filtered with thresholds
subfamily_H3K27ac_sample_filter = subfamily_H3K27ac_sample[which(subfamily_H3K27ac_sample$Length_ijk >= 600 & subfamily_H3K27ac_sample$Length_ik > 5000),]

# Number of members overlapping H3K27ac peak (needs matrix)
TE_H3K27ac_peaks_members = merge(melt(aggregate(data=TE_H3K27ac_peaks[,c(4,8:105)],.~subfamily,function(x) sum(as.numeric(x))),id.vars=c("subfamily")),melt(aggregate(data=TE_H3K27ac_peaks[,c(4,8:105)],.~subfamily,function(x) sum(as.numeric(x) > 0)),id.vars=c("subfamily")),by=c("subfamily","variable"))
colnames(TE_H3K27ac_peaks_members)[2:4] = c("Sample","H3K27ac_peaks","TEs_overlapping")
TE_H3K27ac_peaks_members = merge(TE_H3K27ac_peaks_members,rmsk_TEother_stats_subfamily[,3:4],by.x=c("subfamily"),by.y=c("Subfamily"))
TE_H3K27ac_peaks_members[which(TE_H3K27ac_peaks_members$Sample %in% c("E116","E117","E123","E124","E126","E127")),]$Count = rmsk_TEother_stats_subfamily[match(TE_H3K27ac_peaks_members[which(TE_H3K27ac_peaks_members$Sample %in% c("E116","E117","E123","E124","E126","E127")),]$subfamily,rmsk_TEother_stats_subfamily$Subfamily),]$Count_noY

# Proportion of members overlapping H3K27ac peak
TE_H3K27ac_peaks_members$Percent = TE_H3K27ac_peaks_members$TEs_overlapping/TE_H3K27ac_peaks_members$Count

# Added to enrichment matrix
subfamily_H3K27ac_sample = merge(subfamily_H3K27ac_sample,TE_H3K27ac_peaks_members,by.x=c("Subfamily","Sample"),by.y=c("subfamily","Sample"),all.x=TRUE)
subfamily_H3K27ac_sample[which(subfamily_H3K27ac_sample$Subfamily == "Tigger2a_Car"),c(17:18,20)] = 0
subfamily_H3K27ac_sample[which(subfamily_H3K27ac_sample$Subfamily == "Tigger2a_Car"),]$Count = 2

# Number of samples subfamily is enriched
subfamily_H3K27ac_sample_counts = ddply(subfamily_H3K27ac_sample,.(Class,Family,Subfamily),function(x) sum(x$Enrichment > 1.5 & x$Length_ijk >= 600 & x$Length_ik > 5000))
