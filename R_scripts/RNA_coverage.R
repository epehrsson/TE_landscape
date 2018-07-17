# Genome
genome_RNA = read.table("RNAseq/features/Genome_average.txt",sep='\t')
genome_RNA_chrY = read.table("RNAseq/features/Genome_average_chrY.txt",sep='\t')
genome_RNA[match(genome_RNA_chrY$V1,genome_RNA$V1),2:3] = genome_RNA[match(genome_RNA_chrY$V1,genome_RNA$V1),2:3]-genome_RNA_chrY[,2:3]
colnames(genome_RNA) = c("Sample","Coverage","Cov_length")
genome_RNA$Feature = rep("Genome",dim(genome_RNA)[1])
genome_RNA$Length = ifelse(metadata[match(genome_RNA$Sample,metadata$Sample),]$chrY == "Yes",GENOME_WIDTH,GENOME_WIDTH_noY)

# TE
TEs_RNA = read.table("RNAseq/features/TE_average.txt",sep='\t')
TEs_RNA_chrY = read.table("RNAseq/features/TE_average_chrY.txt",sep='\t')
TEs_RNA[match(TEs_RNA_chrY$V1,TEs_RNA$V1),2:3] = TEs_RNA[match(TEs_RNA_chrY$V1,TEs_RNA$V1),2:3]-TEs_RNA_chrY[,2:3]
colnames(TEs_RNA) = c("Sample","Coverage","Cov_length")
TEs_RNA$Feature = rep("TE",dim(TEs_RNA)[1])
TEs_RNA$Length = ifelse(metadata[match(TEs_RNA$Sample,metadata$Sample),]$chrY == "Yes",MERGED_TE_WIDTH,MERGED_TE_WIDTH_noY)

## Combine
RNA_coverage = rbind(genome_RNA,TEs_RNA)
rm(list=c("genome_RNA","genome_RNA_chrY","TEs_RNA","TEs_RNA_chrY"))

## Average expression
RNA_coverage$Average = RNA_coverage$Coverage/RNA_coverage$Length

## TE relative to genome
RNA_proportion = dcast(RNA_coverage,Sample~Feature,value.var="Average")
RNA_proportion$Ratio = RNA_proportion$TE/RNA_proportion$Genome
RNA_proportion = merge(RNA_proportion,metadata[,c("Sample",sample_categories)])

# Class
class_RNA = read.table("RNAseq/features/TEother_class_merge_average.txt",sep='\t')
class_RNA_chrY = read.table("RNAseq/features/TEother_class_merge_average_chrY.txt",sep='\t')
class_RNA = merge(class_RNA,class_RNA_chrY,by=c("V1","V2"),all.x=TRUE)
class_RNA$Coverage = ifelse(is.na(class_RNA$V3.y),class_RNA$V3.x,class_RNA$V3.x-class_RNA$V3.y)
class_RNA$Cov_length = ifelse(is.na(class_RNA$V4.y),class_RNA$V4.x,class_RNA$V4.x-class_RNA$V4.y)
class_RNA = class_RNA[,c(1:2,7:8)]
colnames(class_RNA) = c("Sample","class","Coverage","Cov_length")
class_RNA = class_RNA[which(!(class_RNA$class %in% c("Unconfident","RC"))),]
class_RNA$class_update = convert_class(class_RNA$class)
class_RNA$class_update = factor(class_RNA$class_update,levels=c("LINE","SINE","LTR","DNA","SVA","Other"))
class_RNA$Length = ifelse(metadata[match(class_RNA$Sample,metadata$Sample),]$chrY == "Yes",
                          rmsk_TE_class[match(class_RNA$class_update,rmsk_TE_class$class_update),]$Total_length,
                          rmsk_TE_class[match(class_RNA$class_update,rmsk_TE_class$class_update),]$Total_length_noY)
class_RNA$Average = class_RNA$Coverage/class_RNA$Length

## Combine with genome
class_RNA = merge(class_RNA,RNA_proportion[,c("Sample","Genome","TE",sample_categories)],by="Sample",all.x=TRUE)
class_RNA$Ratio = class_RNA$Average/class_RNA$Genome

# Subfamily
subfamily_RNA = read.table("RNAseq/features/TEother_subfamily_merge_average.txt",sep='\t')
subfamily_RNA_chrY = read.table("RNAseq/features/TEother_subfamily_merge_average_chrY.txt",sep='\t')
subfamily_RNA = merge(subfamily_RNA,subfamily_RNA_chrY,by=c("V1","V2"),all.x=TRUE)
subfamily_RNA$Coverage = ifelse(is.na(subfamily_RNA$V3.y),subfamily_RNA$V3.x,subfamily_RNA$V3.x-subfamily_RNA$V3.y)
subfamily_RNA$Cov_length = ifelse(is.na(subfamily_RNA$V4.y),subfamily_RNA$V4.x,subfamily_RNA$V4.x-subfamily_RNA$V4.y)
subfamily_RNA = subfamily_RNA[,c(1:2,7:8)]
colnames(subfamily_RNA) = c("Sample","subfamily","Coverage","Cov_length")
subfamily_RNA = merge(subfamily_RNA,expand.grid(subfamily=rmsk_TE_subfamily$subfamily,Sample=metadata[which(!is.na(metadata$RNA)),]$Sample),
                      by=c("subfamily","Sample"),all=TRUE)
subfamily_RNA[is.na(subfamily_RNA)] = 0
subfamily_RNA$Length = ifelse(metadata[match(subfamily_RNA$Sample,metadata$Sample),]$chrY == "Yes",
                          rmsk_TE_subfamily[match(subfamily_RNA$subfamily,rmsk_TE_subfamily$subfamily),]$Total_length,
                          rmsk_TE_subfamily[match(subfamily_RNA$subfamily,rmsk_TE_subfamily$subfamily),]$Total_length_noY)
subfamily_RNA$Average = subfamily_RNA$Coverage/subfamily_RNA$Length

rm(list=c("class_RNA_chrY","subfamily_RNA_chrY"))