# From the matrix of RefSeq feature lengths and lengths of overlap with TEs,
# calculates the proportion of the genome within each feature, 
# the proportion of TE length overlapping each feature (same and either strand), and
# the proportion of each feature overlapping TEs

feature_overlap$Genome_percent_genome = feature_overlap$Genome/GENOME_WIDTH
feature_overlap$TEs_percent_TE = feature_overlap$TEs/MERGED_TE_WIDTH
feature_overlap$TEs_strand_percent_TE = feature_overlap$TEs_strand/MERGED_TE_WIDTH
feature_overlap$Genome_percent_TE = feature_overlap$TEs/feature_overlap$Genome

feature_overlap_long = melt(feature_overlap[,c(2:3,8:11)],id.vars = c("Feature","Coding"))
colnames(feature_overlap_long)[3:4] = c("Measure","Percent")