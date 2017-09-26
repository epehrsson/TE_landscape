# Expression potential for exons

library(reshape2)

load("R_datasets/rna.RData")

# Reduce matrix to unique exons
RNA_exon_unique = unique(RNA_refseq_exon_agnostic[,c(1:3,7:61)])

# Distribution of exons with RPKM >1
RNA_potential_exon = sample_distribution(RNA_exon_unique,57,52)
colnames(RNA_potential_exon)[2] = "Expression"
RNA_potential_exon_cum = cumulative_distribution(RNA_exon_unique,57,52)
RNA_potential_exon_cum$State = rep("Expression",52)
RNA_potential_exon_stats = potential_stats(RNA_potential_exon,1,52)
RNA_potential_exon_stats$State = "Expression"

# Proportion of exons with RPKM >1, by sample
RNA_exon_sample = as.data.frame(apply(RNA_exon_unique[,5:56],2,function(x) sum(x > 1)))
colnames(RNA_exon_sample) = "Count"
RNA_exon_sample$Proportion = RNA_exon_sample$Count/263413
RNA_exon_sample[which(metadata[match(RNA_exon_sample$Sample,metadata$Sample),]$chrY == "No"),]$Proportion = RNA_exon_sample[which(metadata[match(RNA_exon_sample$Sample,metadata$Sample),]$chrY == "No"),]$Count/262307
RNA_exon_sample$Sample = rownames(RNA_exon_sample)
RNA_exon_sample$State = rep("Expression",52)