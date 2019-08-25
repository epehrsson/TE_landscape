# Creates dataframes of the number of exons expressed RPKM > 1 in each number of samples
# As well as the number of exons expressed RPKM > 1 in at least one sample

## RNA_potential_exon: Number of exons expressed RPKM >1 for each number/proportion of samples
## RNA_potential_exon_stats: Proportion of exons ever expressed RPKM > 1 and mean/SE proportion of samples in state
## RNA_exon_sample: Number/proportion of exons expressed RPKM >1 by sample

# Number of unique exons, with and without chrY
NUM_EXON = dim(RNA_refseq_exon)[1]
NUM_EXON_noY = dim(RNA_refseq_exon[which(RNA_refseq_exon$chromosome != "chrY"),])[1]

# Number of exons expressed RPKM >1 for each number of samples (0-56)
RNA_potential_exon = sample_distribution(RNA_refseq_exon,61,sample_counts["All","RNA"])
colnames(RNA_potential_exon)[2] = "Expressed_samples"

# Proportion of exons ever expressed RPKM > 1 and mean/SE proportion of samples in state, for all exons and those ever in the state
RNA_potential_exon_stats = potential_stats(RNA_potential_exon,1,sample_counts["All","RNA"])
RNA_potential_exon_stats$State = "Expressed_samples"

# Number/proportion of exons expressed RPKM >1 by sample
RNA_exon_sample = as.data.frame(apply(RNA_refseq_exon[,5:60],2,function(x) sum(x > 1)))
colnames(RNA_exon_sample) = "Count"
RNA_exon_sample$Proportion = RNA_exon_sample$Count/ifelse(metadata[match(rownames(RNA_exon_sample),metadata$Sample),]$chrY == "Yes",NUM_EXON,NUM_EXON_noY)
RNA_exon_sample$Sample = rownames(RNA_exon_sample)
RNA_exon_sample$State = rep("Expressed_samples",sample_counts["All","RNA"])

# Number of exons expressed RPKM > 1 in each number/proportion of samples
RNA_potential_exon = melt(RNA_potential_exon,id.var="Samples")
colnames(RNA_potential_exon) = c("Samples","State","Count")
RNA_potential_exon$Sample.Proportion = RNA_potential_exon$Samples/(length(RNA_potential_exon$Samples)-1)