# Mouse TE pi chart summaries 
# See 9/23/2016, 2/9/2017, 5/17/2017

# Number of TEs, subfamilies, and length by class
mm9_rmsk_TE_class_pi = merge(merge(aggregate(data=mm9_rmsk_TE,Length ~ Class,FUN=length),aggregate(data=mm9_rmsk_TE,Subfamily ~ Class,function(x) length(unique(x))),by=c("Class")),aggregate(data=mm9_rmsk_TE,Length ~ Class,FUN=sum),by=c("Class"))
colnames(mm9_rmsk_TE_class_pi) = c("Class","Elements","Subfamilies","Length")

# Number of other TEs, subfamilies, and length by class
mm9_rmsk_TEother_class_pi = merge(merge(aggregate(data=rbind(mm9_rmsk_TE,mm9_rmsk_other),Length ~ Class,FUN=length),aggregate(data=rbind(mm9_rmsk_TE,mm9_rmsk_other),Subfamily ~ Class,function(x) length(unique(x))),by=c("Class")),aggregate(data=rbind(mm9_rmsk_TE,mm9_rmsk_other),Length ~ Class,FUN=sum),by=c("Class"))
colnames(mm9_rmsk_TEother_class_pi) = c("Class","Elements","Subfamilies","Length")

# Number of TEs, subfamilies, and length by class, mm10
mm10_rmsk_TE_class_pi = merge(merge(aggregate(data=mm10_rmsk_TE,Length ~ class,FUN=length),aggregate(data=mm10_rmsk_TE,subfamily ~ class,function(x) length(unique(x))),by=c("class")),aggregate(data=mm10_rmsk_TE,Length ~ class,FUN=sum),by=c("class"))
colnames(mm10_rmsk_TE_class_pi) = c("Class","Elements","Subfamilies","Length")
test = c("Unconfident",colSums(mm10_rmsk_TE_class_pi[which(mm10_rmsk_TE_class_pi$Class %in% c("DNA?","LINE?","LTR?","SINE?","Unknown","Unknown?")),2:4]))
names(test)[1] = "Class"
mm10_rmsk_TE_class_pi$Class = factor(mm10_rmsk_TE_class_pi$Class,levels=c("DNA","LINE","LTR","SINE","RC","Unconfident"))
mm10_rmsk_TE_class_pi = rbind(mm10_rmsk_TE_class_pi[c(1,3,5,7,9),],test)
mm10_rmsk_TE_class_pi[,2:4] = apply(mm10_rmsk_TE_class_pi[,2:4],2,function(x) as.numeric(x))
mm10_rmsk_TE_class_pi = mm10_rmsk_TE_class_pi[c(1:3,5,4,6),]
