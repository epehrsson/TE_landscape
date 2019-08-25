# Creates a dataframe with summary statistics for each TE class ("rmsk_TE_class")

# From the dataframe of all individual TEs, the number of TEs, families, and subfamilies,
# Median/sd of length (bp), median mappability, and median age (Jukes-Cantor evolutionary distance)
rmsk_TE_class = ddply(rmsk_TE,~class_update,summarize,
                      Count = length(Length),
                      Families = length(unique(family)), 
                      Subfamilies = length(unique(subfamily)), 
                      Median_length = median(Length), 
                      SD_length = sd(Length), 
                      Mappability = median(mappability), 
                      Age = median(JC_distance))

# Number of TEs, excluding chrY
rmsk_TE_class = merge(rmsk_TE_class,ddply(rmsk_TE[which(rmsk_TE$chromosome != "chrY"),],~class_update,summarize,Count_noY = length(Length)),by="class_update")

# Proportion of class members overlapping each genic feature
rmsk_TE_class = merge(rmsk_TE_class,ddply(rmsk_TE,~class_update,function(y) apply(y[,15:35],2,function(x) length(na.omit(x))/length(x))),by="class_update")

# Total number of CpGs per class
TE_class_CpG_count = read.table("WGBS/class/TE_CpG_class.txt",sep='\t',row.names=1)
colnames(TE_class_CpG_count) = c("CpGs")
rownames(TE_class_CpG_count)[c(12,14)] = c("SVA","Other")
rmsk_TE_class$CpGs = TE_class_CpG_count[as.vector(rmsk_TE_class$class_update),]/2
rm(TE_class_CpG_count)

# Proportion of all CpGs and CpGs overlapping TEs in each class
rmsk_TE_class$Percent_TE_CpGs = rmsk_TE_class$CpGs/TE_CPGS
rmsk_TE_class$Percent_all_CpGs = rmsk_TE_class$CpGs/ALL_CPGS

# Number of TEs with CpGs
TE_CpG_count = read.table("WGBS/TE_CpG_count.txt",sep='\t')
colnames(TE_CpG_count) = c(TE_coordinates[c(1:4,6,5,7)],"CpGs")
test = table(TE_CpG_count$class)
names(test)[7] = "SVA"
test[13] = sum(test[c(2,4,6,8,10:12)])
names(test)[13] = "Other"
rmsk_TE_class$TEs_wCpG = test[as.vector(rmsk_TE_class$class_update)]
rm(list=c("TE_CpG_count","test"))

# Proportion of TEs with CpGs
rmsk_TE_class$TEs_wCpG_per = rmsk_TE_class$TEs_wCpG/rmsk_TE_class$Count

# Total length of class (unique bp), with and without chrY
class_lengths = cbind(read.table("features/TEs/class/class_lengths.txt",sep='\t'),read.table("features/TEs/class/class_lengths_noY.txt",sep='\t'))[c(1:3,5,7:8),c(1:2,4)]
colnames(class_lengths) = c("class","Total_length","Total_length_noY")
class_lengths$class = factor(class_lengths$class,levels=c(levels(class_lengths$class),"SVA"))
class_lengths[which(class_lengths$class == "Other"),]$class = "SVA"
class_lengths[which(class_lengths$class == "Unconfident_RC"),]$class = "Other"
rmsk_TE_class$Total_length = class_lengths[match(rmsk_TE_class$class_update,class_lengths$class),]$Total_length
rmsk_TE_class$Total_length_noY = class_lengths[match(rmsk_TE_class$class_update,class_lengths$class),]$Total_length_noY

# Total length of each class across all samples with each type of epigenetic data
rmsk_TE_class$chromHMM_total_width = as.numeric(sample_counts["chrY","chromHMM"])*rmsk_TE_class$Total_length + as.numeric(sample_counts["All","chromHMM"]-sample_counts["chrY","chromHMM"])*rmsk_TE_class$Total_length_noY
rmsk_TE_class$WGBS_total_width = as.numeric(sample_counts["chrY","WGBS"])*rmsk_TE_class$Total_length + as.numeric(sample_counts["All","WGBS"]-sample_counts["chrY","WGBS"])*rmsk_TE_class$Total_length_noY
rmsk_TE_class$DNase_total_width = as.numeric(sample_counts["chrY","DNase"])*rmsk_TE_class$Total_length + as.numeric(sample_counts["All","DNase"]-sample_counts["chrY","DNase"])*rmsk_TE_class$Total_length_noY
rmsk_TE_class$H3K27ac_total_width = as.numeric(sample_counts["chrY","H3K27ac"])*rmsk_TE_class$Total_length + as.numeric(sample_counts["All","H3K27ac"]-sample_counts["chrY","H3K27ac"])*rmsk_TE_class$Total_length_noY
rm(class_lengths)

# Rearrange class order
rmsk_TE_class = rmsk_TE_class[c(1:3,5,6,4),]