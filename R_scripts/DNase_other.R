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

# Proportion
# Proportion of RefSeq features overlapping DNase peaks, averaged
DNase_features = merge(DNase_features,aggregate(data=mnemonics_states_features_noTE,Bases~Cohort+Sample,sum),by=c("Sample","Cohort"))
colnames(DNase_features)[3:4] = c("Total_width_in_feature","Feature_width")
DNase_features$Proportion = DNase_features$Total_width_in_feature/DNase_features$Feature_width
DNase_proportion = aggregate(data=DNase_features,Proportion~Cohort,mean)

# Proportion of TEs overlapping DNase peaks, averaged
DNase_proportion[20,] = c("Genome",mean(DNase_stats$Total_width/as.numeric(mnemonics_states_genome[1,as.vector(DNase_stats$Sample)])))
DNase_proportion[21,] = c("TE",mean(DNase_stats$Total_width_in_TE/as.numeric(mnemonics_states_TE[1,as.vector(DNase_stats$Sample)])))
DNase_proportion$Cohort = factor(DNase_proportion$Cohort,levels=c("TE","Genome","genome_noTE","promoters","promoters_pc","promoters_nc","5UTR","5UTR_pc","5UTR_nc","coding_exon","coding_exon_pc","3UTR","3UTR_pc","3UTR_nc","exons","exons_pc","exons_nc","introns","introns_pc","introns_nc","intergenic"))
DNase_proportion$Proportion = as.numeric(DNase_proportion$Proportion)

# Contribution
# Proportion of all peaks in TEs
sum(DNase_stats$Peaks_in_TE)/sum(DNase_stats$Peaks)

# Total DNase width in TEs, all samples
sum(as.numeric(DNase_stats$Total_width_in_TE))/sum(as.numeric(DNase_stats$Total_width))

# Proportion of Dnase peaks overlapping TEs by sample classification
mean(DNase_stats$Total_width_in_TE/DNase_stats$Total_width)
kruskal.test(DNase_stats$Total_width_in_TE/DNase_stats$Total_width,EID_metadata[match(DNase_stats$Sample,EID_metadata$Sample),]$Type)
wilcox_to_all(DNase_stats$Total_width_in_TE/DNase_stats$Total_width,droplevels(EID_metadata[match(DNase_stats$Sample,EID_metadata$Sample),]$Group))