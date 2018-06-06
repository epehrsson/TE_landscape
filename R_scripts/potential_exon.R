# Expression potential for exons

NUM_EXON = dim(RNA_refseq_exon)[1]
NUM_EXON_noY = dim(RNA_refseq_exon[which(RNA_refseq_exon$chromosome != "chrY"),])[1]

# Distribution of exons with RPKM >1
RNA_potential_exon = sample_distribution(RNA_refseq_exon,61,sample_counts["All","RNA"])
colnames(RNA_potential_exon)[2] = "Expressed_samples"
RNA_potential_exon_cum = cumulative_distribution(RNA_refseq_exon,61,sample_counts["All","RNA"])
RNA_potential_exon_cum$State = rep("Expressed_samples",sample_counts["All","RNA"])
RNA_potential_exon_stats = potential_stats(RNA_potential_exon,1,sample_counts["All","RNA"])
RNA_potential_exon_stats$State = "Expressed_samples"

# Proportion of exons with RPKM >1, by sample
RNA_exon_sample = as.data.frame(apply(RNA_refseq_exon[,5:60],2,function(x) sum(x > 1)))
colnames(RNA_exon_sample) = "Count"
RNA_exon_sample$Proportion = RNA_exon_sample$Count/ifelse(metadata[match(rownames(RNA_exon_sample),metadata$Sample),]$chrY == "Yes",NUM_EXON,NUM_EXON_noY)
RNA_exon_sample$Sample = rownames(RNA_exon_sample)
RNA_exon_sample$State = rep("Expressed_samples",sample_counts["All","RNA"])

RNA_potential_exon_long = melt(RNA_potential_exon,id.var="Samples")
colnames(RNA_potential_exon_long) = c("Samples","State","Count")
RNA_potential_exon_long$Sample.Proportion = RNA_potential_exon_long$Samples/(length(RNA_potential_exon_long$Samples)-1)
