# TE statistics
# See 4/19/2016, 4/25/2016, 8/24/2016, 9/21/2016, 9/27/2016, 2/3/2017, 2/9/2017, 5/16/2017, 7/21/2017

# RepeatMasker file restricted to four TE classes
rmsk_TE_strand = read.table(file="TE_landscape/rmsk_TE.txt",sep='\t')
colnames(rmsk_TE_strand) = c("Chromosome","Start","Stop","Subfamily","Class","Family","Strand")
rmsk_TE_merge = merge(rmsk_TE,rmsk_TE_strand,by=c("Chromosome","Start","Stop","Subfamily","Class","Family"))
rmsk_TE = rmsk_TE_merge

# Median TE lengths
median(rmsk_TE$Length)

# Number of instances, median length, sd of length for each class
merge(merge(aggregate(data=rmsk_TE,Length ~ Class,FUN=length),aggregate(data=rmsk_TE,Length ~ Class,FUN=median),by=c("Class")),aggregate(data=rmsk_TE,Length ~ Class,FUN=sd),by=c("Class"))
colnames(test) = c("Class","Count","Median","SD")

# Number of instances, median length, sd of length for each subfamily
merge(merge(aggregate(data=rmsk_TE,Length ~ Class+Family+Subfamily,FUN=length),aggregate(data=rmsk_TE,Length ~ Class+Family+Subfamily,FUN=median),by=c("Class","Family","Subfamily")),aggregate(data=rmsk_TE,Length ~ Class+Family+Subfamily,FUN=sd),by=c("Class","Family","Subfamily"))
colnames(test)[4:6] = c("Count","Median","SD")

# Number of instances for each subfamily without Y chromosome
rmsk_TE_stats_subfamily = merge(rmsk_TE_stats_subfamily,aggregate(data=rmsk_TE[which(rmsk_TE$Chromosome != "chrY"),],Length ~ Class+Family+Subfamily,FUN=length),by=c("Class","Family","Subfamily"))
colnames(rmsk_TE_stats_subfamily)[7] = "Count_noY"

# Total length of subfamily
rmsk_TE_stats_subfamily$Total_length = test[match(rmsk_TE_stats_subfamily$Subfamily,test$Subfamily),]$Length_ik

# RepeatMasker file restricted to other TE classes
rmsk_other = read.table(file='other/rmsk_other.txt',sep='\t')
colnames(rmsk_other) = colnames(rmsk_TE)[c(1:6,8)]
rmsk_other$Length = rmsk_other$Stop-rmsk_other$Start

# Median TE lengths
median(rbind(rmsk_TE,rmsk_other)$Length)

# Number of instances, median length, sd of length for each class
rmsk_other_stats_class = merge(merge(aggregate(data=rmsk_other,Length ~ Class,FUN=length),aggregate(data=rmsk_other,Length ~ Class,FUN=median),by=c("Class")),aggregate(data=rmsk_other,Length ~ Class,FUN=sd),by=c("Class"))
colnames(rmsk_other_stats_class) = c("Class","Count","Median","SD")

# Number of instances, median length, sd of length for each subfamily
rmsk_other_stats_subfamily = merge(merge(merge(aggregate(data=rmsk_other,Length ~ Class+Family+Subfamily,FUN=length),aggregate(data=rmsk_other,Length ~ Class+Family+Subfamily,FUN=median),by=c("Class","Family","Subfamily")),aggregate(data=rmsk_other,Length ~ Class+Family+Subfamily,FUN=sd),by=c("Class","Family","Subfamily")),aggregate(data=rmsk_other[which(rmsk_other$Chromosome != "chrY"),],Length ~ Class+Family+Subfamily,FUN=length),by=c("Class","Family","Subfamily"))
colnames(rmsk_other_stats_subfamily)[4:7] = c("Count","Median","SD","Count_noY")

# Total length of subfamily
rmsk_other_stats_subfamily$Total_length = test[match(rmsk_other_stats_subfamily$Subfamily,test$Subfamily),]$Length_ik

# RepeatMasker file restricted to all TE classes
rmsk_TEother = rbind(rmsk_TE,rmsk_other)
rmsk_TEother$Class_update = rmsk_TEother$Class
rmsk_TEother$Class_update = factor(rmsk_TEother$Class_update,levels=c("DNA","LINE","SINE","LTR","RC","Other","Unconfident"))
rmsk_TEother[which(rmsk_TEother$Class %in% c("DNA?","LINE?","SINE?","LTR?","Unknown","Unknown?")),]$Class_update = "Unconfident"

# Number of instances, median length, sd of length for each class, all TEs
rmsk_TEother_stats_class = rbind(rmsk_TE_stats_class,rmsk_other_stats_class)

# Number of instances, median length, sd of length for each subfamily, all TEs
rmsk_TEother_stats_subfamily = rbind(rmsk_TE_stats_subfamily,rmsk_other_stats_subfamily)
