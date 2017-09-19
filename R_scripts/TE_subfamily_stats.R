# TE statistics by subfamily
# See 4/19/2016, 4/25/2016, 8/24/2016, 8/25/2016, 9/20/2016, 9/21/2016, 9/27/2016, 9/28/2016, 11/4/2016, 11/5/2016, 11/7/2016, 11/18/2016, 12/16/2016, 1/31/2017, 2/1/2017, 2/3/2017, 2/6/2017, 2/9/2017, 2/10/2017, 2/25/2017, 2/27/2017, 2/28/2017, 3/5/2017, 3/8/2017, 5/14/2017, 5/15/2017, 5/16/2017, 5/17/2017, 6/7/2017, 6/14/2017, 6/15/2017, 7/21/2017, 7/24/2017, 8/1/2017, 8/2/2017

load("R_datasets/rmsk_TE.RData")

# Number of instances, median/sd of length, mean/sd mappability and age by subfamily
rmsk_TE_subfamily = ddply(rmsk_TE,~subfamily+family+class_update,summarize,Count = length(Length),Median_length = median(Length), SD_length = sd(Length), Mappability = mean(mappability), Mappability_SD = sd(mappability), Age = mean(JC_distance), Age_SD = sd(JC_distance))
rmsk_TE_subfamily = merge(rmsk_TE_subfamily,ddply(rmsk_TE[which(rmsk_TE$chromosome != "chrY"),],~subfamily+family+class_update,summarize,Count_noY = length(Length)),by=c("subfamily","family","class_update"))

# Proportion of each subfamily in each feature
rmsk_TE_subfamily = merge(rmsk_TE_subfamily,ddply(rmsk_TE,~subfamily,function(y) apply(y[,13:32],2,function(x) length(na.omit(x))/length(x))),by="subfamily")

# Number of CpGs per subfamily
TE_subfamily_CpG_count = read.table("WGBS/subfamily/TE_CpG_subfamily.txt",sep='\t')
colnames(TE_subfamily_CpG_count) = c("subfamily","CpGs")
rmsk_TE_subfamily = merge(rmsk_TE_subfamily,TE_subfamily_CpG_count,by="subfamily",all.x=TRUE)
rmsk_TE_subfamily[which(is.na(rmsk_TE_subfamily$CpGs)),]$CpGs = 0
rm(TE_subfamily_CpG_count)

# Total length of subfamily
subfamily_length = cbind(read.table("features/TEs/subfamily/subfamily_lengths.txt",sep='\t'),read.table("features/TEs/subfamily/subfamily_lengths_noY.txt",sep='\t'))[,c(1:2,4)]
colnames(subfamily_length) = c("subfamily","Total_length","Total_length_noY")
rmsk_TE_subfamily = merge(rmsk_TE_subfamily,subfamily_length,by="subfamily")
rm(subfamily_length)

rmsk_TE_subfamily$CpGs_per_TE = rmsk_TE_subfamily$CpGs/rmsk_TE_subfamily$Count
rmsk_TE_subfamily$CpGs_per_kbp = rmsk_TE_subfamily$CpGs/(rmsk_TE_subfamily$Total_length/1000)