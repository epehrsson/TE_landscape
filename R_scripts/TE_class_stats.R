# TE statistics by class
# See 4/19/2016, 4/25/2016, 8/24/2016, 8/25/2016, 9/20/2016, 9/21/2016, 9/23/2016, 9/27/2016, 9/28/2016, 11/4/2016, 
# 11/5/2016, 11/7/2016, 11/18/2016, 12/16/2016, 1/31/2017, 2/1/2017, 2/3/2017, 2/6/2017, 2/9/2017, 2/10/2017, 2/25/2017, 
# 2/27/2017, 2/28/2017, 3/5/2017, 3/8/2017, 5/12/2017, 5/14/2017, 5/15/2017, 5/16/2017, 5/17/2017, 6/7/2017, 6/14/2017, 6/15/2017, 7/21/2017, 7/24/2017, 8/1/2017, 8/2/2017

#load("R_datasets/rmsk_TE.RData")

# Number of instances, median/sd of length, mean/sd mappability, mean age/sd for each class
rmsk_TE_class = ddply(rmsk_TE,~class_update,summarize,
                      Count = length(Length),
                      Families = length(unique(family)), 
                      Subfamilies = length(unique(subfamily)), 
                      Median_length = median(Length), 
                      SD_length = sd(Length), 
                      Mappability = mean(mappability), 
                      Mappability_SD = sd(mappability), 
                      Age = mean(JC_distance), 
                      Age_SD = sd(JC_distance))
rmsk_TE_class = merge(rmsk_TE_class,ddply(rmsk_TE[which(rmsk_TE$chromosome != "chrY"),],~class_update,summarize,Count_noY = length(Length)),by="class_update")

# Proportion of each class overlapping each feature
rmsk_TE_class = merge(rmsk_TE_class,ddply(rmsk_TE,~class_update,function(y) apply(y[,13:33],2,function(x) length(na.omit(x))/length(x))),by="class_update")

# Number of CpGs per class
TE_class_CpG_count = read.table("WGBS/class/TE_CpG_class.txt",sep='\t',row.names=1)
colnames(TE_class_CpG_count) = c("CpGs")
rownames(TE_class_CpG_count)[c(12,14)] = c("SVA","Other")
rmsk_TE_class$CpGs = TE_class_CpG_count[as.vector(rmsk_TE_class$class_update),]
rm(TE_class_CpG_count)

# Proportion of CpGs (all, overlapping TEs) in each class
rmsk_TE_class$Percent_TE_CpGs = rmsk_TE_class$CpGs/TE_CPGS
rmsk_TE_class$Percent_all_CpGs = rmsk_TE_class$CpGs/ALL_CPGS

# Number of CpGs per TE
TE_CpG_count = read.table("WGBS/TE_CpG_count.txt",sep='\t')
colnames(TE_CpG_count) = c("chromosome","start","stop","subfamily","class","family","strand","CpGs")
test = table(TE_CpG_count$class)
names(test)[7] = "SVA"
test[13] = sum(test[c(2,4,6,8,10:12)])
names(test)[13] = "Other"
rmsk_TE_class$TEs_wCpG = test[as.vector(rmsk_TE_class$class_update)]
rm(list=c("TE_CpG_count","test"))

# CpGs per TE, number/proportion of TEs with CpGs
rmsk_TE_class$TEs_wCpG_per = rmsk_TE_class$TEs_wCpG/rmsk_TE_class$Count
rmsk_TE_class$Mean_CpG = rmsk_TE_class$CpGs/rmsk_TE_class$Count
rmsk_TE_class$Mean_CpG_wCpG = rmsk_TE_class$CpGs/rmsk_TE_class$TEs_wCpG

# Total length per class
class_lengths = cbind(read.table("features/TEs/class/class_lengths.txt",sep='\t'),read.table("features/TEs/class/class_lengths_noY.txt",sep='\t'))[c(1:3,5,7:8),c(1:2,4)]
colnames(class_lengths) = c("class","Total_length","Total_length_noY")
class_lengths$class = factor(class_lengths$class,levels=c(levels(class_lengths$class),"SVA"))
class_lengths[which(class_lengths$class == "Other"),]$class = "SVA"
class_lengths[which(class_lengths$class == "Unconfident_RC"),]$class = "Other"
rmsk_TE_class$Total_length = class_lengths[match(rmsk_TE_class$class_update,class_lengths$class),]$Total_length
rmsk_TE_class$Total_length_noY = class_lengths[match(rmsk_TE_class$class_update,class_lengths$class),]$Total_length_noY
rmsk_TE_class$chromHMM_total_width = as.numeric(sample_counts["chrY","chromHMM"])*rmsk_TE_class$Total_length + as.numeric(sample_counts["All","chromHMM"]-sample_counts["chrY","chromHMM"])*rmsk_TE_class$Total_length_noY
rmsk_TE_class$WGBS_total_width = as.numeric(sample_counts["chrY","WGBS"])*rmsk_TE_class$Total_length + as.numeric(sample_counts["All","WGBS"]-sample_counts["chrY","WGBS"])*rmsk_TE_class$Total_length_noY
rmsk_TE_class$DNase_total_width = as.numeric(sample_counts["chrY","DNase"])*rmsk_TE_class$Total_length + as.numeric(sample_counts["All","DNase"]-sample_counts["chrY","DNase"])*rmsk_TE_class$Total_length_noY
rmsk_TE_class$H3K27ac_total_width = as.numeric(sample_counts["chrY","H3K27ac"])*rmsk_TE_class$Total_length + as.numeric(sample_counts["All","H3K27ac"]-sample_counts["chrY","H3K27ac"])*rmsk_TE_class$Total_length_noY
rm(class_lengths)

rmsk_TE_class = rmsk_TE_class[c(1:3,5,6,4),]