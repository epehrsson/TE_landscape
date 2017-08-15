# Mouse TE statistics
# See 5/19/2016, 2/9/2017, 5/17/2017

# RepeatMasker file for mm9
mm9_rmsk_TE = read.table("mm9_rmsk_TE.txt",sep='\t')
colnames(mm9_rmsk_TE) = c("Chromosome","Start","Stop","Subfamily","Class","Family","Strand")
mm9_rmsk_TE$Length = mm9_rmsk_TE$Stop - mm9_rmsk_TE$Start

# RepeatMasker file for mm9, restricted to other TEs
mm9_rmsk_other = read.table("Mouse/mm9_rmsk_other.txt",sep='\t')
colnames(mm9_rmsk_other) = c("Chromosome","Start","Stop","Subfamily","Class","Family","Strand")
mm9_rmsk_other$Length = mm9_rmsk_other$Stop - mm9_rmsk_other$Start

# Number of instances, median, sd by class
mm9_rmsk_TE_stats = merge(merge(aggregate(data=mm9_rmsk_TE,Length ~ Class,FUN=length),aggregate(data=mm9_rmsk_TE,Length ~ Class,FUN=median),by=c("Class")),aggregate(data=mm9_rmsk_TE,Length ~ Class,FUN=sd),by=c("Class"))
colnames(mm9_rmsk_TE_stats) = c("Class","Count","Median","SD")

# Number of instances, median, sd by class
mm9_rmsk_other_stats = merge(merge(aggregate(data=mm9_rmsk_other,Length ~ Class,FUN=length),aggregate(data=mm9_rmsk_other,Length ~ Class,FUN=median),by=c("Class")),aggregate(data=mm9_rmsk_other,Length ~ Class,FUN=sd),by=c("Class"))
colnames(mm9_rmsk_other_stats) = c("Class","Count","Median","SD")

# RepeatMasker file for mm10
mm10_rmsk_TE = read.table("mm10_rmsk_TE.txt",sep='\t')
colnames(mm10_rmsk_TE) = c("chromosome","start","stop","subfamily","class","family","strand")
mm10_rmsk_TE$Length = mm10_rmsk_TE$stop - mm10_rmsk_TE$start
