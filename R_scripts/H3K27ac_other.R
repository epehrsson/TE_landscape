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

# Proportion
# Proportion of RefSeq features overlapping H3K27ac peaks, averaged
H3K27ac_features = merge(H3K27ac_features,aggregate(data=mnemonics_states_features_noTE,Bases~Cohort+Sample,sum),by=c("Sample","Cohort"))
colnames(H3K27ac_features)[3:4] = c("Total_width_in_feature","Feature_width")
H3K27ac_features$Proportion = H3K27ac_features$Total_width_in_feature/H3K27ac_features$Feature_width
H3K27ac_proportion = aggregate(data=H3K27ac_features,Proportion~Cohort,mean)

# Proportion of TEs overlapping H3K27ac peaks, averaged
H3K27ac_proportion[20,] = c("Genome",mean(H3K27ac_stats$Total_width/as.numeric(mnemonics_states_genome[1,as.vector(H3K27ac_stats$Sample)])))
H3K27ac_proportion[21,] = c("TE",mean(H3K27ac_stats$Total_width_in_TE/as.numeric(mnemonics_states_TE[1,as.vector(H3K27ac_stats$Sample)])))
H3K27ac_proportion$Cohort = factor(H3K27ac_proportion$Cohort,levels=c("TE","Genome","genome_noTE","promoters","promoters_pc","promoters_nc","5UTR","5UTR_pc","5UTR_nc","coding_exon","coding_exon_pc","3UTR","3UTR_pc","3UTR_nc","exons","exons_pc","exons_nc","introns","introns_pc","introns_nc","intergenic"))
H3K27ac_proportion$Proportion = as.numeric(H3K27ac_proportion$Proportion)

# Contribution
# Proportion of all peaks in TEs
sum(as.numeric(H3K27ac_stats$Total_width_in_TE))/sum(as.numeric(H3K27ac_stats$Total_width))

# By-sample
# Proportion of H3K27ac peaks overlapping TEs by sample classification
mean(H3K27ac_stats$Total_width_in_TE/H3K27ac_stats$Total_width) 
kruskal.test(H3K27ac_stats$Total_width_in_TE/H3K27ac_stats$Total_width,EID_metadata[match(H3K27ac_stats$Sample,EID_metadata$Sample),]$Group)
wilcox_to_all(H3K27ac_stats$Total_width_in_TE/H3K27ac_stats$Total_width,droplevels(EID_metadata[match(H3K27ac_stats$Sample,EID_metadata$Sample),]$Group))