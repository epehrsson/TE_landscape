# Overlap with features (bp)
# See 8/24/2016, 8/25/2016, 9/28/2016, 2/6/2017, 2/9/2017, 6/15/2017

# Proportion of genome, TEs overlapping each RefSeq feature and vice-versa
feature_overlap$Genome_percent_genome = feature_overlap$Genome/GENOME_WIDTH
feature_overlap$TEs_percent_TE = feature_overlap$TEs/MERGED_TE_WIDTH
feature_overlap$TEs_strand_percent_TE = feature_overlap$TEs_strand/MERGED_TE_WIDTH
feature_overlap$Genome_percent_TE = feature_overlap$TEs/feature_overlap$Genome

feature_overlap_long = melt(feature_overlap[,c(2:3,8:11)],id.vars = c("Feature","Coding"))
colnames(feature_overlap_long)[3:4] = c("Measure","Percent")
