# DNase contribution, proportion overlapping peaks, proportion of peaks overlapping TEs by sample grouping
# Proportion of Refseq feature overlapping peaks
# See 5/24/2017, 5/30/2017, 6/5/2017, 8/3/2017, 8/7/2017

# Number and width of peaks per sample, overall and in TEs 
DNase_stats = read.table("DNase/DNase_stats.txt",sep='\t')
colnames(DNase_stats) = c("Peaks","Sample","Total_width","Peaks_in_TE")
test = read.table("DNase/rmsk_TEother_merge_DNase_contribution.txt",sep='\t')
colnames(test) = c("Sample","Total_width_in_TE")
DNase_stats = merge(DNase_stats,test,by=c("Sample"))[,c(1:2,4,3,5)]
rm(test)

# Number of bases of DNase peak in Refseq genic features, by sample
DNase_features = read.table("DNase/Refseq_features/refseq_features_DNase.txt",sep='\t')
colnames(DNase_features) = c("Sample","Bases","Cohort")
DNase_features$Cohort = gsub("refseq_", "", DNase_features$Cohort)
DNase_features$Cohort = gsub("_merge_noTE_DNase", "", DNase_features$Cohort)
DNase_features$Cohort = gsub("_noTE_DNase", "", DNase_features$Cohort)
DNase_features = DNase_features[which(!(DNase_features$Cohort %in% c("genes","genes_nc","genes_pc"))),]