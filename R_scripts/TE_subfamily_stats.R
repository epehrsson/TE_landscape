# TE statistics by subfamily

# Number of instances, median/sd of length, mean/sd mappability and age by subfamily
rmsk_TE_subfamily = ddply(rmsk_TE,~subfamily+family+class_update,summarize,
                          Count = length(Length),
                          Length = median(Length), 
                          Mappability = median(mappability), 
                          Age = median(JC_distance))
rmsk_TE_subfamily = merge(rmsk_TE_subfamily,
                          ddply(rmsk_TE[which(rmsk_TE$chromosome != "chrY"),],~subfamily+family+class_update,
                                summarize,Count_noY = length(Length)),by=c("subfamily","family","class_update"))

# Total length of subfamily
subfamily_length = cbind(read.table("features/TEs/subfamily/subfamily_lengths.txt",sep='\t'),
                         read.table("features/TEs/subfamily/subfamily_lengths_noY.txt",sep='\t'))[,c(1:2,4)]
colnames(subfamily_length) = c("subfamily","Total_length","Total_length_noY")
rmsk_TE_subfamily = merge(rmsk_TE_subfamily,subfamily_length,by="subfamily")
rm(subfamily_length)

# Number of CpGs per subfamily
TE_subfamily_CpG_count = read.table("WGBS/subfamily/TE_CpG_subfamily.txt",sep='\t')
colnames(TE_subfamily_CpG_count) = c("subfamily","CpGs")
rmsk_TE_subfamily = merge(rmsk_TE_subfamily,TE_subfamily_CpG_count,by="subfamily",all.x=TRUE)
rmsk_TE_subfamily[which(is.na(rmsk_TE_subfamily$CpGs)),]$CpGs = 0
rmsk_TE_subfamily$CpGs = rmsk_TE_subfamily$CpGs/2
rm(TE_subfamily_CpG_count)

rmsk_TE_subfamily$CpGs_per_TE = rmsk_TE_subfamily$CpGs/rmsk_TE_subfamily$Count
rmsk_TE_subfamily$CpGs_per_kbp = rmsk_TE_subfamily$CpGs/(rmsk_TE_subfamily$Total_length/1000)

# Number of TEs with CpGs
TE_CpG_count = read.table("WGBS/TE_CpG_count.txt",sep='\t')
colnames(TE_CpG_count) = c(TE_coordinates[c(1:4,6,5,7)],"CpGs")
subfamily_wCpG_count = aggregate(data=TE_CpG_count,CpGs~subfamily,length)
colnames(subfamily_wCpG_count)[2] = "Count_CpGs"
rmsk_TE_subfamily = merge(rmsk_TE_subfamily,subfamily_wCpG_count,by="subfamily",all.x=TRUE)
rmsk_TE_subfamily[is.na(rmsk_TE_subfamily)] = 0
rm(list=c("TE_CpG_count","subfamily_wCpG_count"))

# Proportion of each subfamily in each feature
rmsk_TE_subfamily = merge(rmsk_TE_subfamily,ddply(rmsk_TE,~subfamily,function(y) 
  apply(y[,cohorts],2,function(x) length(na.omit(x))/length(x))),by="subfamily")

# Proportion of each subfamily on each chromosome
rmsk_TE_chromosome = ddply(rmsk_TE,.(subfamily,chromosome),summarise,Chr=length(chromosome))
rmsk_TE_chromosome = ddply(rmsk_TE_chromosome,.(subfamily),transform,Chr_percent=Chr/sum(Chr))
rmsk_TE_chromosome = dcast(rmsk_TE_chromosome,subfamily~chromosome,value.var="Chr_percent")
rmsk_TE_chromosome[is.na(rmsk_TE_chromosome)] = 0
rmsk_TE_subfamily = merge(rmsk_TE_subfamily,rmsk_TE_chromosome,by="subfamily")
rm(rmsk_TE_chromosome)
