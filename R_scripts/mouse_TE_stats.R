# Mouse TE statistics
# See 5/19/2016, 9/23/2016, 2/9/2017, 5/17/2017

# RepeatMasker file for mm10
mm10_rmsk_TE = read.table("mm10_rmsk_TE.txt",sep='\t')
colnames(mm10_rmsk_TE) = c("chromosome","start","stop","subfamily","class","family","strand")
mm10_rmsk_TE$Length = mm10_rmsk_TE$stop - mm10_rmsk_TE$start

# Number of TEs, subfamilies, and length by class, mm10
mm10_rmsk_TE_class_pi = merge(merge(aggregate(data=mm10_rmsk_TE,Length ~ class,FUN=length),aggregate(data=mm10_rmsk_TE,subfamily ~ class,function(x) length(unique(x))),by=c("class")),aggregate(data=mm10_rmsk_TE,Length ~ class,FUN=sum),by=c("class"))
colnames(mm10_rmsk_TE_class_pi) = c("Class","Elements","Subfamilies","Length")
test = c("Unconfident",colSums(mm10_rmsk_TE_class_pi[which(mm10_rmsk_TE_class_pi$Class %in% c("DNA?","LINE?","LTR?","SINE?","Unknown","Unknown?")),2:4]))
names(test)[1] = "Class"
mm10_rmsk_TE_class_pi$Class = factor(mm10_rmsk_TE_class_pi$Class,levels=c("DNA","LINE","LTR","SINE","RC","Unconfident"))
mm10_rmsk_TE_class_pi = rbind(mm10_rmsk_TE_class_pi[c(1,3,5,7,9),],test)
mm10_rmsk_TE_class_pi[,2:4] = apply(mm10_rmsk_TE_class_pi[,2:4],2,function(x) as.numeric(x))
mm10_rmsk_TE_class_pi = mm10_rmsk_TE_class_pi[c(1:3,5,4,6),]