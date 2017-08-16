# H3K27ac contribution, proportion of peaks overlapping TEs by sample grouping
# Proportion of Refseq feature overlapping peaks
# See 7/4/2017, 8/3/2017, 8/7/2017

# Number and width of peaks per sample, overall and in TEs 
H3K27ac_stats = read.table("H3K27ac/H3K27ac_stats.txt",sep='\t',header=TRUE)
test = read.table("H3K27ac/rmsk_TEother_merge_H3K27ac_contribution.txt",sep='\t')
colnames(test) = c("Sample","Total_width_in_TE")
H3K27ac_stats = merge(H3K27ac_stats,test,by=c("Sample"))
rm(test)

# Number of bases of H3K27ac peak in Refseq genic features, by sample
H3K27ac_features = read.table("H3K27ac/Refseq_features/refseq_features_H3K27ac.txt",sep='\t')
colnames(H3K27ac_features) = c("Sample","Bases","Cohort")
H3K27ac_features$Cohort = gsub("refseq_", "", H3K27ac_features$Cohort)
H3K27ac_features$Cohort = gsub("_merge_noTE_H3K27ac", "", H3K27ac_features$Cohort)
H3K27ac_features$Cohort = gsub("_noTE_H3K27ac", "", H3K27ac_features$Cohort)