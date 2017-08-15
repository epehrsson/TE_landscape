# Overlap with features (bp)
# See 8/24/2016, 8/25/2016, 9/28/2016, 2/6/2017, 2/9/2017, 6/15/2017

# Proportion of genome, TEs overlapping each RefSeq feature and vice-versa
feature_overlap_cast = read.table("feature_overlap.txt",sep='\t',header=TRUE)
feature_overlap_cast$Genome_percent_genome = feature_overlap_cast$Genome/3095693983
feature_overlap_cast$TEs_percent_TE = feature_overlap_cast$TEs/1389947349
feature_overlap_cast$TEs_strand_percent_TE = feature_overlap_cast$TEs_strand/1389947349
feature_overlap_cast$Genome_percent_TE = feature_overlap_cast$TEs/feature_overlap_cast$Genome
feature_overlap_long = melt(feature_overlap_cast,id.vars = c("Feature","Cohort"))
colnames(feature_overlap_long)[3:4] = c("Measure","Percent")
feature_overlap_long = feature_overlap_long[which(feature_overlap_long$Measure %in% colnames(feature_overlap_cast)[6:9]),]
feature_overlap_long$Cohort = factor(feature_overlap_long$Cohort,levels=levels(feature_overlap_long$Cohort)[c(1,3,2)])
feature_overlap_long$Feature = factor(feature_overlap_long$Feature,levels=levels(feature_overlap_long$Feature)[c(7,2,3,1,4,6,5)])
