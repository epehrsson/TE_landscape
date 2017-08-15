# TE pi chart summaries
# See 9/23/2016, 9/28/2016, 2/9/2017, 3/8/2017, 5/12/2017, 5/17/2017

# Number of TEs, subfamilies, and length by class
rmsk_TE_class_pi = merge(merge(aggregate(data=rmsk_TE,Length ~ Class,FUN=length),aggregate(data=rmsk_TE,Subfamily ~ Class,function(x) length(unique(x))),by=c("Class")),aggregate(data=rmsk_TE,Length ~ Class,FUN=sum),by=c("Class"))
colnames(rmsk_TE_class_pi) = c("Class","Elements","Subfamilies","Length")
rmsk_TE_class_pi$Mouse_ortholog = table(human_mouse_orthologs$Class_human)

# Number of all TEs, subfamilies, and length by class
rmsk_TEother_class_pi = merge(merge(aggregate(data=rmsk_TEother,Length ~ Class,FUN=length),aggregate(data=rmsk_TEother,Subfamily ~ Class,function(x) length(unique(x))),by=c("Class")),aggregate(data=rmsk_TEother,Length ~ Class,FUN=sum),by=c("Class"))
colnames(rmsk_TEother_class_pi) = c("Class","Elements","Subfamilies","Length") 
rmsk_TEother_class_pi = merge(rmsk_TEother_class_pi,as.data.frame(table(human_mouse_orthologs_mm10$human_class)),by.x=c("Class"),by.y=c("Var1"),all.x=TRUE)

# Adding mouse orthologs
colnames(rmsk_TEother_class_pi)[5] = "Mouse_ortholog"
rmsk_TEother_class_pi[which(is.na(rmsk_TEother_class_pi$Mouse_ortholog)),]$Mouse_ortholog = 0
rmsk_TEother_class_pi = merge(rmsk_TEother_class_pi,aggregate(data=rmsk_TEother_stats_subfamily[which(rmsk_TEother_stats_subfamily$Subfamily %in% mm10_rmsk_TE$subfamily),],Subfamily~Class,length),by=c("Class"),all.x=TRUE)
colnames(rmsk_TEother_class_pi)[6] = "Mouse_ortholog_subfamily"
rmsk_TEother_class_pi[which(is.na(rmsk_TEother_class_pi$Mouse_ortholog_subfamily)),]$Mouse_ortholog_subfamily = 0

# Adding Unconfident class
test = c("Unconfident",apply(rmsk_TEother_class_pi[which(rmsk_TEother_class_pi$Class %in% c("DNA?","LINE?","SINE?","LTR?","Unknown","Unknown?")),2:6],2,sum))
names(test)[1] = "Class"
rmsk_TEother_class_pi$Class = factor(rmsk_TEother_class_pi$Class,levels=c(levels(rmsk_TEother_class_pi$Class),"Unconfident"))
rmsk_TEother_class_pi = rbind(rmsk_TEother_class_pi,test)
