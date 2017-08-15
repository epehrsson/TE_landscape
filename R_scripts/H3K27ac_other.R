# H3K27ac contribution, proportion of peaks overlapping TEs by sample grouping
# Proportion of Refseq feature overlapping peaks
# See 7/4/2017, 8/3/2017, 8/7/2017

# Number and width of peaks per sample, overall and in TEs 
H3K27ac_stats = read.table("H3K27ac/H3K27ac_stats.txt",sep='\t',header=TRUE)
test = read.table("H3K27ac/rmsk_TEother_merge_H3K27ac_contribution.txt",sep='\t')
colnames(test) = c("Sample","Total_width_in_TE")
H3K27ac_stats = merge(H3K27ac_stats,test,by=c("Sample"))

# Contribution
# Proportion of all peaks in TEs
sum(as.numeric(H3K27ac_stats$Total_width_in_TE))/sum(as.numeric(H3K27ac_stats$Total_width))

# By-sample
# Proportion of H3K27ac peaks overlapping TEs by sample classification
mean(H3K27ac_stats$Total_width_in_TE/H3K27ac_stats$Total_width) 
kruskal.test(H3K27ac_stats$Total_width_in_TE/H3K27ac_stats$Total_width,EID_metadata[match(H3K27ac_stats$Sample,EID_metadata$Sample),]$Group)
wilcox_to_all(H3K27ac_stats$Total_width_in_TE/H3K27ac_stats$Total_width,droplevels(EID_metadata[match(H3K27ac_stats$Sample,EID_metadata$Sample),]$Group))

# Number of bases of H3K27ac peak in Refseq genic features, by sample
H3K27ac_features = read.table("H3K27ac/Refseq_features/refseq_features_H3K27ac.txt",sep='\t')
colnames(H3K27ac_features) = c("Sample","Bases","Feature")
