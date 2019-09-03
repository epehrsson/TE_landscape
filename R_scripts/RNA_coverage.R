# RNA-seq read coverage (total and average) over the entire genome, TEs, and each TE class and subfamily, by sample

## RNA_coverage - RNA-seq read coverage (total and average) over the entire genome and TEs, by sample
## RNA_proportion - Ratio of average read coverage over TEs compared to the entire genome, by sample
## class_RNA - RNA-seq read coverage (total and average) over each TE class, by sample, 
## and the ratio of average read coverage over each class compared to the entire genome
## subfamily_RNA - RNA-seq read coverage (total and average) over each TE subfamily, by sample

# Load dataframe of total RNA-seq read coverage and length of coverage for the entire genome, by sample
genome_RNA = read.table("RNAseq/features/Genome_average.txt",sep='\t',col.names=c("Sample","Coverage","Cov_length"))
## Total RNA-seq read coverage and length of coverage on chrY for four samples without chrY
genome_RNA_chrY = read.table("RNAseq/features/Genome_average_chrY.txt",sep='\t',col.names=c("Sample","Coverage","Cov_length"))
## Remove chrY coverage for four samples without chrY
genome_RNA[match(genome_RNA_chrY$Sample,genome_RNA$Sample),2:3] = genome_RNA[match(genome_RNA_chrY$Sample,genome_RNA$Sample),2:3]-genome_RNA_chrY[,2:3]
genome_RNA$Feature = rep("Genome",dim(genome_RNA)[1])
## Add length of genome
genome_RNA$Length = ifelse(metadata[match(genome_RNA$Sample,metadata$Sample),]$chrY == "Yes",GENOME_WIDTH,GENOME_WIDTH_noY)

# Load dataframe of total RNA-seq read coverage and length of coverage over TEs, by sample
TEs_RNA = read.table("RNAseq/features/TE_average.txt",sep='\t',col.names=c("Sample","Coverage","Cov_length"))
## Total RNA-seq read coverage and length of coverage on chrY for four samples without chrY
TEs_RNA_chrY = read.table("RNAseq/features/TE_average_chrY.txt",sep='\t',col.names=c("Sample","Coverage","Cov_length"))
## Remove chrY coverage for four samples without chrY
TEs_RNA[match(TEs_RNA_chrY$Sample,TEs_RNA$Sample),2:3] = TEs_RNA[match(TEs_RNA_chrY$Sample,TEs_RNA$Sample),2:3]-TEs_RNA_chrY[,2:3]
TEs_RNA$Feature = rep("TE",dim(TEs_RNA)[1])
## Add length of TEs
TEs_RNA$Length = ifelse(metadata[match(TEs_RNA$Sample,metadata$Sample),]$chrY == "Yes",MERGED_TE_WIDTH,MERGED_TE_WIDTH_noY)

# Combine into a single data frame
RNA_coverage = rbind(genome_RNA,TEs_RNA)
rm(list=c("genome_RNA","genome_RNA_chrY","TEs_RNA","TEs_RNA_chrY"))

# Average read coverage over the entire genome and TEs, by sample
RNA_coverage$Average = RNA_coverage$Coverage/RNA_coverage$Length

# Ratio of average read coverage over TEs compared to the entire genome
RNA_proportion = dcast(RNA_coverage,Sample~Feature,value.var="Average")
RNA_proportion$Ratio = RNA_proportion$TE/RNA_proportion$Genome
## Add sample metadata
RNA_proportion = merge(RNA_proportion,metadata[,c("Sample",sample_categories)])


# Load dataframe of total RNA-seq read coverage and length of coverage over each TE class, by sample
class_RNA = read.table("RNAseq/features/TEother_class_merge_average.txt",sep='\t',col.names=c("Sample","class","Coverage","Cov_length"))
## Total RNA-seq read coverage and length of coverage on chrY for four samples without chrY
class_RNA_chrY = read.table("RNAseq/features/TEother_class_merge_average_chrY.txt",sep='\t',col.names=c("Sample","class","Coverage","Cov_length"))
## Remove chrY coverage for four samples without chrY
class_RNA = merge(class_RNA,class_RNA_chrY,by=c("Sample","class"),all.x=TRUE)
class_RNA$Coverage = ifelse(is.na(class_RNA$Coverage.y),class_RNA$Coverage.x,class_RNA$Coverage.x-class_RNA$Coverage.y)
class_RNA$Cov_length = ifelse(is.na(class_RNA$Cov_length.y),class_RNA$Cov_length.x,class_RNA$Cov_length.x-class_RNA$Cov_length.y)
class_RNA = class_RNA[,c("Sample","class","Coverage","Cov_length")]
## Remove deprecated classes
class_RNA = class_RNA[which(!(class_RNA$class %in% c("Unconfident","RC"))),]
## Update class assignments
class_RNA$class_update = convert_class(class_RNA$class)
class_RNA$class_update = factor(class_RNA$class_update,levels=c("LINE","SINE","LTR","DNA","SVA","Other"))
## Add length of class
class_RNA$Length = ifelse(metadata[match(class_RNA$Sample,metadata$Sample),]$chrY == "Yes",
                          rmsk_TE_class[match(class_RNA$class_update,rmsk_TE_class$class_update),]$Total_length,
                          rmsk_TE_class[match(class_RNA$class_update,rmsk_TE_class$class_update),]$Total_length_noY)

# Average RNA-seq read coverage over each TE class, by sample
class_RNA$Average = class_RNA$Coverage/class_RNA$Length

# Combine average RNA-seq read coverage by class and for all TEs/entire genome
class_RNA = merge(class_RNA,RNA_proportion[,c("Sample","Genome","TE",sample_categories)],by="Sample",all.x=TRUE)

# Ratio of average read coverage over each TE class to coverage over the entire genome
class_RNA$Ratio = class_RNA$Average/class_RNA$Genome


# Load dataframe of total RNA-seq read coverage and length of coverage over each TE subfamily, by sample
subfamily_RNA = read.table("RNAseq/features/TEother_subfamily_merge_average.txt",sep='\t',col.names=c("Sample","subfamily","Coverage","Cov_length"))
## Total RNA-seq read coverage and length of coverage on chrY for four samples without chrY
subfamily_RNA_chrY = read.table("RNAseq/features/TEother_subfamily_merge_average_chrY.txt",sep='\t',col.names=c("Sample","subfamily","Coverage","Cov_length"))
## Remove chrY coverage for four samples without chrY
subfamily_RNA = merge(subfamily_RNA,subfamily_RNA_chrY,by=c("Sample","subfamily"),all.x=TRUE)
subfamily_RNA$Coverage = ifelse(is.na(subfamily_RNA$Coverage.y),subfamily_RNA$Coverage.x,subfamily_RNA$Coverage.x-subfamily_RNA$Coverage.y)
subfamily_RNA$Cov_length = ifelse(is.na(subfamily_RNA$Cov_length.y),subfamily_RNA$Cov_length.x,subfamily_RNA$Cov_length.x-subfamily_RNA$Cov_length.y)
subfamily_RNA = subfamily_RNA[,c("Sample","subfamily","Coverage","Cov_length")]
## Add subfamily x sample combinations with no RNA-seq coverage
subfamily_RNA = merge(subfamily_RNA,expand.grid(subfamily=rmsk_TE_subfamily$subfamily,Sample=metadata[which(!is.na(metadata$RNA)),]$Sample),
                      by=c("subfamily","Sample"),all=TRUE)
subfamily_RNA[is.na(subfamily_RNA)] = 0
## Add length of subfamily
subfamily_RNA$Length = ifelse(metadata[match(subfamily_RNA$Sample,metadata$Sample),]$chrY == "Yes",
                          rmsk_TE_subfamily[match(subfamily_RNA$subfamily,rmsk_TE_subfamily$subfamily),]$Total_length,
                          rmsk_TE_subfamily[match(subfamily_RNA$subfamily,rmsk_TE_subfamily$subfamily),]$Total_length_noY)

# Average RNA-seq read coverage over each TE subfamily, by sample
subfamily_RNA$Average = subfamily_RNA$Coverage/subfamily_RNA$Length

rm(list=c("class_RNA_chrY","subfamily_RNA_chrY"))