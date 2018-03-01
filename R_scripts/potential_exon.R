# Expression potential for exons

NUM_EXON = dim(RNA_exon_unique)[1]
NUM_EXON_noY = dim(RNA_exon_unique[which(RNA_exon_unique$chromosome != "chrY"),])[1]

# Distribution of exons with RPKM >1
RNA_potential_exon = sample_distribution(RNA_exon_unique,58,sample_counts["All","RNA"])
colnames(RNA_potential_exon)[2] = "Expression"
RNA_potential_exon_cum = cumulative_distribution(RNA_exon_unique,58,sample_counts["All","RNA"])
RNA_potential_exon_cum$State = rep("Expression",sample_counts["All","RNA"])
RNA_potential_exon_stats = potential_stats(RNA_potential_exon,1,sample_counts["All","RNA"])
RNA_potential_exon_stats$State = "Expression"

# Proportion of exons with RPKM >1, by sample
RNA_exon_sample = as.data.frame(apply(RNA_exon_unique[,6:57],2,function(x) sum(x > 1)))
colnames(RNA_exon_sample) = "Count"
RNA_exon_sample$Proportion = RNA_exon_sample$Count/ifelse(metadata[match(RNA_exon_sample$Sample,metadata$Sample),]$chrY == "Yes",NUM_EXON,NUM_EXON_noY)
RNA_exon_sample$Sample = rownames(RNA_exon_sample)
RNA_exon_sample$State = rep("Expression",sample_counts["All","RNA"])