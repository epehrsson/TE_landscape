# DNase contribution, proportion overlapping peaks, proportion of peaks overlapping TEs by sample grouping
# Proportion of Refseq feature overlapping peaks
# See 5/24/2017, 5/30/2017, 6/5/2017, 8/3/2017, 8/7/2017

# Number of peaks per sample, overall and in TEs 
DNase_stats = read.table("DNase_peaks/DNase_stats.txt",sep='\t')
colnames(DNase_stats) = c("Peaks","Sample","Total_width","Peaks_in_TE")
DNase_stats = DNase_stats[,c(2,1,4,3)]

# Contribution
# Proportion of all peaks in TEs
sum(DNase_stats$Peaks_in_TE)/sum(DNase_stats$Peaks)

# Total DNase width in TEs, all samples
test = read.table("DNase_peaks/rmsk_TEother_merge_DNase_contribution.txt",sep='\t')
colnames(test) = c("Sample","Total_width_in_TE")
DNase_stats = merge(DNase_stats,test,by=c("Sample"))
sum(as.numeric(DNase_stats$Total_width_in_TE))/sum(as.numeric(DNase_stats$Total_width))

# Proportion of Dnase peaks overlapping TEs by sample classification
mean(DNase_stats$Total_width_in_TE/DNase_stats$Total_width)
kruskal.test(DNase_stats$Total_width_in_TE/DNase_stats$Total_width,EID_metadata[match(DNase_stats$Sample,EID_metadata$Sample),]$Type)
wilcox_to_all(DNase_stats$Total_width_in_TE/DNase_stats$Total_width,droplevels(EID_metadata[match(DNase_stats$Sample,EID_metadata$Sample),]$Group))

# Proportion of TEs overlapping DNase peaks, by sample
mean(DNase_stats$Total_width_in_TE/mnemonics_states_TEother[as.vector(DNase_stats$Sample),]$Total)

# Number of bases of Dnase peak in Refseq genic features, by sample
DNase_features = read.table("DNase/Refseq_features/refseq_features_DNase.txt",sep='\t')
colnames(DNase_features) = c("Sample","Bases","Feature")