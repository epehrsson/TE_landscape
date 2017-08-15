# Age analyses
# See 11/5/2016, 11/7/2016, 11/18/2016, 12/16/2016, 2/1/2017, 2/6/2017, 2/9/2017, 2/10/2017, 2/27/2017, 2/28/2017, 3/5/2017, 3/8/2017, 5/14/2017, 5/17/2017, 6/7/2017, 6/14/2017, 6/15/2017, 7/24/2017, 8/1/2017, 8/2/2017

# Individual
# JC evolutionary distances
rmsk_TE_age = read.table("rmsk_TE_JCage.txt",sep='\t')
colnames(rmsk_TE_age) = c("chromosome","start","stop","subfamily","family","class","strand","substitutions","JC_distance")

# Combined with number of samples in each state, hypomethylated
rmsk_TE_age = merge(rmsk_TE_age,potential_TE_state[,1:22],by=c("chromosome","start","stop","subfamily","family","class","strand"))
rmsk_TE_age = merge(rmsk_TE_age,TE_meth_wCpG_average[,c(1:7,45)],by=c("chromosome","start","stop","subfamily","family","class","strand"),all.x=TRUE)
rmsk_TE_age$Hypo <- NULL
rmsk_TE_age = merge(rmsk_TE_age,rmsk_TE_feature[,c(1:7,9:15,17)],by=c("chromosome","start","stop","subfamily","family","class","strand"))

# Correlation between individual TE age, number of samples hypomethylated
unlist(cor.test(rmsk_TE_age$JC_distance,rmsk_TE_age$Hypomethylated))[c("p.value","estimate.cor")]
by(rmsk_TE_age,rmsk_TE_age$class,function(x) unlist(cor.test(x$JC_distance,x$Hypomethylated))[c("p.value","estimate.cor")])

# Subfamilies enriched in 3'UTR and hypomethylated
sort(table(droplevels(rmsk_TE_age[which(!(is.na(rmsk_TE_age$'3UTR')) & rmsk_TE_age$Hypomethylated > 10),]$subfamily)))

# JC evolutionary distances, other TEs
rmsk_other_age = read.table("other/rmsk_other_JCage.txt",sep='\t')
colnames(rmsk_other_age) = c("chromosome","start","stop","subfamily","family","class","strand","substitutions","JC_distance")
rmsk_other_age = merge(rmsk_other_age,potential_TEother_state[,1:22],by=c("chromosome","start","stop","subfamily","family","class","strand"))
rmsk_other_age = merge(rmsk_other_age,other_meth_wCpG_average[,c(1:7,45)],by=c("chromosome","start","stop","subfamily","family","class","strand"),all.x=TRUE)
rmsk_other_age = merge(rmsk_other_age,rmsk_other_feature[,c(1:7,9:15,17)],by=c("chromosome","start","stop","subfamily","family","class","strand"))

# JC evolutionary distances, all TEs
rmsk_TEother_age = rbind(rmsk_TE_age,rmsk_other_age)

# Proportion of TEs ever in state by genic feature
apply(rmsk_TEother_age[,25:32],2,function(x) apply(rmsk_TEother_age[which(!(is.na(x))),10:24],2,mean))/127

apply(rmsk_TEother_age[,c(10:24,33)],2,function(x) unlist(cor.test(rmsk_TEother_age$JC_distance,x))[c("p.value","estimate.cor")])
by(rmsk_TEother_age,rmsk_TEother_age$class,function(y) apply(y[,c(10:24,33)],2,function(x) unlist(cor.test(y$JC_distance,x))[c("p.value","estimate.cor")]))
apply(rmsk_TEother_age[,25:32],2,function(x) unlist(wilcox.test(rmsk_TEother_age[which(is.na(x)),]$JC_distance,rmsk_TEother_age[which(!(is.na(x))),]$JC_distance))["p.value"])
by(rmsk_TEother_age,rmsk_TEother_age$class,function(y) apply(y[,25:32],2,function(x) unlist(wilcox.test(y[which(is.na(x)),]$JC_distance,y[which(!(is.na(x))),]$JC_distance))["p.value"]))
apply(rmsk_TEother_age[,25:32],2,function(x) mean(rmsk_TEother_age[which(is.na(x)),]$JC_distance) - mean(rmsk_TEother_age[which(!(is.na(x))),]$JC_distance))

# Adding Unconfident class
rmsk_TEother_age$class_update = rmsk_TEother_age$class
rmsk_TEother_age$class_update = factor(rmsk_TEother_age$class_update,levels=c(levels(rmsk_TEother_age$class)[c(1:4,8:9)],"Unconfident"))
rmsk_TEother_age[which(rmsk_TEother_age$class %in% c("DNA?","LINE?","LTR?","SINE?","Unknown","Unknown?")),]$class_update = "Unconfident"

# Adding number of samples overlapping DNase peak
rmsk_TEother_age = merge(rmsk_TEother_age,TE_DNase_peaks[,c(1:7,61)],by=c("chromosome","start","stop","subfamily","family","class","strand"),all.x=TRUE)
colnames(rmsk_TEother_age)[35] = "DNase"
rmsk_TEother_age[which(is.na(rmsk_TEother_age$DNase)),]$DNase = 0

# Correlation between Dnase, age
cor.test(rmsk_TEother_age$JC_distance,rmsk_TEother_age$DNase)
by(rmsk_TEother_age,rmsk_TEother_age$class_update,function(y) cor.test(y$JC_distance,y$DNase))

# Updating features
rmsk_TEother_age = merge(rmsk_TEother_age[,c(1:24,33:35)],rmsk_TEother_feature,by=c("chromosome","start","stop","subfamily","family","class","strand","class_update"))

# Correlation of age and feature overlap
apply(rmsk_TEother_age[,28:35],2,function(x) unlist(wilcox.test(rmsk_TEother_age[which(is.na(x)),]$JC_distance,rmsk_TEother_age[which(!(is.na(x))),]$JC_distance))["p.value"])
by(rmsk_TEother_age,rmsk_TEother_age$class_update,function(y) apply(y[,28:35],2,function(x) unlist(wilcox.test(y[which(is.na(x)),]$JC_distance,y[which(!(is.na(x))),]$JC_distance))["p.value"]))

# Average TE, Alu, SVA age
mean(rmsk_TEother_age$JC_distance)
mean(rmsk_TEother_age[which(rmsk_TEother_age$family == "Alu"),]$JC_distance)
mean(rmsk_TEother_age[which(rmsk_TEother_age$family == "Other"),]$JC_distance)

# Average age of TEs with mouse orthologs
mean(merge(rmsk_TEother_age,unique(human_mouse_orthologs_mm10[,8:14]),by.x=c("chromosome","start","stop","subfamily","family","class","strand"),by.y=c("human_chr_hg19","human_start_hg19","human_stop_hg19","human_subfamily","human_family","human_class","human_strand_hg19"))$JC_distance)

# Adding H3K27ac
rmsk_TEother_age = merge(rmsk_TEother_age,TE_H3K27ac_peaks[,c(1:7,106)],by=c("chromosome","start","stop","subfamily","family","class","strand"),all.x=TRUE)
colnames(rmsk_TEother_age)[36] = "H3K27ac"
rmsk_TEother_age[which(is.na(rmsk_TEother_age$H3K27ac)),]$H3K27ac = 0

# Updating hypomethylation
rmsk_TEother_age = merge(rmsk_TEother_age,TE_meth_average[,c(1:7,46)],by=c("chromosome","start","stop","subfamily","family","class","strand"),all.x=TRUE)
rmsk_TEother_age$Hypomethylated.x <- NULL
colnames(rmsk_TEother_age)[36] = "Hypomethylated"

# Correlation between number of samples hypomethylated or overlapping H3K27ac, age
apply(rmsk_TEother_age[,35:36],2,function(x) unlist(cor.test(rmsk_TEother_age$JC_distance,x))[c("p.value","estimate.cor")])
by(rmsk_TEother_age,rmsk_TEother_age$class_update,function(y) apply(y[,35:36],2,function(x) unlist(cor.test(y$JC_distance,x))[c("p.value","estimate.cor")]))

# Correlation between number of subfamilies hypomethylated (besides IMR90) and age
rmsk_TEother_age = merge(rmsk_TEother_age,TE_meth_average[,c(1:7,50)],by=c("chromosome","start","stop","subfamily","family","class","strand"),all.x=TRUE)
cor.test(rmsk_TEother_age$JC_distance,rmsk_TEother_age$Hypomethylated_noIMR90)
by(rmsk_TEother_age,rmsk_TEother_age$class_update,function(y) unlist(cor.test(y$JC_distance,y$Hypomethylated_noIMR90))[c("p.value","estimate.cor")])

# Subfamily
# Mean JC evolutionary distance per subfamily
rmsk_TE_age_subfamily = merge(aggregate(data=rmsk_TE_age,JC_distance~subfamily+family+class,mean),aggregate(data=rmsk_TE_age,JC_distance~subfamily+family+class,sd),by=c("subfamily","family","class"))
colnames(rmsk_TE_age_subfamily)[4:5] = c("JC_distance_mean","JC_distance_sd")

# Combined with mean methylation, CpGs, CpGs per TE, CpGs per length
rmsk_TE_age_subfamily = merge(rmsk_TE_age_subfamily,TE_meth_subfamily[,c(1:3,41:44)],by=c("subfamily","family","class"),all.x=TRUE)
rmsk_TE_age_subfamily$Meth_mean_noIMR90 <- NULL
rmsk_TE_age_subfamily = merge(rmsk_TE_age_subfamily,TE_meth_subfamily_CpG_count,by=c("subfamily","family","class"))
colnames(rmsk_TE_age_subfamily)[6:7] = c("Meth_mean_noIMR90","Total_CpGs")
rmsk_TE_age_subfamily$CpGs_per_TE = rmsk_TE_age_subfamily$Total_CpGs/rmsk_TE_age_subfamily$TEs
test = subfamily_sample[which(subfamily_sample$Sample == "E001"),c(1:3,5)]
colnames(test) = c("class","family","subfamily","Length_ik")
rmsk_TE_age_subfamily = merge(rmsk_TE_age_subfamily,test,by=c("subfamily","family","class"))
rmsk_TE_age_subfamily$CpGs_per_kbp = rmsk_TE_age_subfamily$Total_CpGs/(rmsk_TE_age_subfamily$Length_ik/1000)

# Correlation between age, subfamily average methylation (all samples)
unlist(cor.test(rmsk_TE_age_subfamily$JC_distance_mean,rmsk_TE_age_subfamily$Mean))[c("p.value","estimate.cor")]
by(rmsk_TE_age_subfamily,rmsk_TE_age_subfamily$class,function(x) unlist(cor.test(x$JC_distance_mean,x$Mean))[c("p.value","estimate.cor")])

# Mean methylation of LTR subfamilies
mean(rmsk_TE_age_subfamily[which(rmsk_TE_age_subfamily$class == "LTR"),]$Mean)
sd(rmsk_TE_age_subfamily[which(rmsk_TE_age_subfamily$class == "LTR"),]$Mean)

# Correlation between age, subfamily average methylation for Alu subfamilies
cor.test(rmsk_TE_age_subfamily[which(rmsk_TE_age_subfamily$family == "Alu"),]$JC_distance_mean,rmsk_TE_age_subfamily[which(rmsk_TE_age_subfamily$family == "Alu"),]$Mean)

# Correlation between age, number of samples enriched in state
cor.test(rmsk_TE_age_subfamily$JC_distance_mean,subfamily_state_sample_counts[match(rmsk_TE_age_subfamily$subfamily,subfamily_state_sample_counts[which(subfamily_state_sample_counts$State== "7_Enh"),]$Subfamily),]$V1)

# Correlation between age, likelihood enriched in state
wilcox.test(rmsk_TE_age_subfamily[which(subfamily_state_sample_counts[match(rmsk_TE_age_subfamily$subfamily,subfamily_state_sample_counts[which(subfamily_state_sample_counts$State== "6_EnhG"),]$Subfamily),]$V1 == 0),]$JC_distance_mean,rmsk_TE_age_subfamily[which(subfamily_state_sample_counts[match(rmsk_TE_age_subfamily$subfamily,subfamily_state_sample_counts[which(subfamily_state_sample_counts$State== "6_EnhG"),]$Subfamily),]$V1 >0),]$JC_distance_mean)

# Correlation between age, number of members
cor.test(rmsk_TE_age_subfamily$JC_distance_mean,rmsk_TE_age_subfamily$TEs)
by(data=rmsk_TE_age_subfamily,rmsk_TE_age_subfamily$class,function(x) cor.test(x$JC_distance_mean,x$TEs))

# Average JC evolutionary distance by subfamily, other TEs
rmsk_other_age_subfamily = merge(aggregate(data=rmsk_other_age,JC_distance~subfamily+family+class,mean),aggregate(data=rmsk_other_age,JC_distance~subfamily+family+class,sd),by=c("subfamily","family","class"))
colnames(rmsk_other_age_subfamily)[4:5] = c("JC_distance_mean","JC_distance_sd")
rmsk_other_age_subfamily = merge(rmsk_other_age_subfamily,other_meth_subfamily_CpG_count[,1:5],by=c("subfamily","family","class"))
colnames(rmsk_other_age_subfamily)[6] = c("Total_CpGs")
rmsk_other_age_subfamily$CpGs_per_TE = rmsk_other_age_subfamily$Total_CpGs/rmsk_other_age_subfamily$TEs
rmsk_other_age_subfamily = merge(rmsk_other_age_subfamily,other_meth_subfamily[,c(1:3,41:44)],by=c("subfamily","family","class"),all.x=TRUE)
test = subfamily_sample[which(subfamily_sample$Sample == "E001"),c(1:3,5)]
colnames(test) = c("class","family","subfamily","Length_ik")
rmsk_other_age_subfamily = merge(rmsk_other_age_subfamily,test,by=c("subfamily","family","class"))
rmsk_other_age_subfamily$CpGs_per_kbp = rmsk_other_age_subfamily$Total_CpGs/(rmsk_other_age_subfamily$Length_ik/1000)

# Repeated matrix creation without error
rmsk_TEother_age_subfamily = rbind(rmsk_TE_age_subfamily[,c(1:7,9,11:13)],rmsk_other_age_subfamily[,c(1:10,14)])

# Adding length of subfamily
rmsk_TEother_age_subfamily = merge(rmsk_TEother_age_subfamily,rmsk_TEother_stats_subfamily[,c(3,8)],by.x=c("subfamily"),by.y=c("Subfamily"))
 
# Adding proportion of subfamily in genic feature
rmsk_TEother_age_subfamily = merge(rmsk_TEother_age_subfamily,aggregate(data=rmsk_TEother_feature,genic~subfamily,function(x) sum(x > 0)/length(x)),by=c("subfamily"))

# Adding proportion of subfamily in promoter
rmsk_TEother_age_subfamily = merge(rmsk_TEother_age_subfamily,aggregate(data=rmsk_TEother_feature,Promoter~subfamily,function(x) length(na.omit(x))/length(x),na.action=na.pass),by=c("subfamily"))

# Adding number of samples enriched in chromHMM states
rmsk_TEother_age_subfamily = merge(rmsk_TEother_age_subfamily,dcast(subfamily_state_sample_counts,Class+Family+Subfamily~State,value.var=c("V1"))[,3:18],by.x="subfamily",by.y="Subfamily",all.x=TRUE)

# Adding mean hypomethylation proportion
rmsk_TEother_age_subfamily = merge(rmsk_TEother_age_subfamily,TEother_meth_wCpG_subfamily_hypo[,c(1:3,44)],by=c("subfamily","family","class"),all.x=TRUE)
colnames(rmsk_TEother_age_subfamily)[30] = "Mean_hypomethylated"
colnames(rmsk_TEother_age_subfamily)[10] = "Mean"
rmsk_TEother_age_subfamily[is.na(rmsk_TEother_age_subfamily)] = 0

# Adding Unconfident class
rmsk_TEother_age_subfamily$class_update = rmsk_TEother_age_subfamily$class
rmsk_TEother_age_subfamily$class_update = factor(rmsk_TEother_age_subfamily$class_update,levels=c(levels(rmsk_TEother_age_subfamily$class)[c(1:4,8:9)],"Unconfident"))
rmsk_TEother_age_subfamily[which(rmsk_TEother_age_subfamily$class %in% c("DNA?","LINE?","LTR?","SINE?","Unknown","Unknown?")),]$class_update = "Unconfident"
rmsk_TEother_age_subfamily[is.na(rmsk_TEother_age_subfamily)] = 0

# Adding Dnase enrichment
rmsk_TEother_age_subfamily = merge(rmsk_TEother_age_subfamily,subfamily_DNase_sample_counts[,3:4],by.x="subfamily",by.y="Subfamily",all.x=TRUE)
colnames(rmsk_TEother_age_subfamily)[32] = "DNase"

# Correlation between age, number of Dnase enrichments
cor.test(rmsk_TEother_age_subfamily$JC_distance_mean,rmsk_TEother_age_subfamily$DNase) 
by(rmsk_TEother_age_subfamily,rmsk_TEother_age_subfamily$class_update,function(x) unlist(cor.test(x$JC_distance_mean,x$DNase))[c("p.value","estimate.cor")])
wilcox.test(rmsk_TEother_age_subfamily[which(rmsk_TEother_age_subfamily$DNase == 0),]$JC_distance_mean,rmsk_TEother_age_subfamily[which(rmsk_TEother_age_subfamily$DNase > 0),]$JC_distance_mean)
rmsk_TEother_age_subfamily = rmsk_TEother_age_subfamily[,c(1:12,15:32)]

# Average age of subfamilies also in mouse 
mean(rmsk_TEother_age_subfamily[which(rmsk_TEother_age_subfamily$subfamily %in% mm10_rmsk_TE$subfamily),]$JC_distance_mean)

# Adding max hypomethylated proportion
rmsk_TEother_age_subfamily = merge(rmsk_TEother_age_subfamily,TEother_meth_wCpG_subfamily_hypo[,c(1:3,45)],by=c("subfamily","family","class"),all.x=TRUE)
colnames(rmsk_TEother_age_subfamily)[31] = "Max_hypo"

# Correlation between age, max hypomethylated proportion
cor.test(rmsk_TEother_age_subfamily$JC_distance_mean,rmsk_TEother_age_subfamily$Max_hypo)
by(rmsk_TEother_age_subfamily,rmsk_TEother_age_subfamily$class_update,function(x) cor.test(x$JC_distance_mean,x$Max_hypo))

# Adding max hypomethylated proportion, no IMR90
rmsk_TEother_age_subfamily = merge(rmsk_TEother_age_subfamily,TEother_meth_wCpG_subfamily_hypo[,c(1:3,46)],by=c("subfamily","family","class"),all.x=TRUE)
colnames(rmsk_TEother_age_subfamily)[31:32] = c("Max_hypo","Max_noIMR90_hypo")

# Correlation between age, max hypomethylated proportion
cor.test(rmsk_TEother_age_subfamily$JC_distance_mean,rmsk_TEother_age_subfamily$Max_noIMR90_hypo)
by(rmsk_TEother_age_subfamily,rmsk_TEother_age_subfamily$class_update,function(x) cor.test(x$JC_distance_mean,x$Max_noIMR90_hypo))

# Updating hypomethylated proportion, CpGs per subfamily
rmsk_TEother_age_subfamily$Mean <- NULL
rmsk_TEother_age_subfamily$Mean_noIMR90 <- NULL
rmsk_TEother_age_subfamily = rmsk_TEother_age_subfamily[,c(1:5,11:25,27:28)]
rmsk_TEother_age_subfamily = merge(rmsk_TEother_age_subfamily,TE_subfamily_CpG_count,by.x=c("subfamily","family","class"),by.y=c("Subfamily","Family","Class"))
rmsk_TEother_age_subfamily = merge(rmsk_TEother_age_subfamily,TE_meth_subfamily_hypo[,c(1:3,41:42,45:46)],by=c("subfamily","family","class"),all.x=TRUE)
colnames(rmsk_TEother_age_subfamily)[29:30] = c("Mean_hypo","Mean_hypo_noIMR90")
colnames(rmsk_TEother_age_subfamily)[31:32] = c("Max_hypo","Max_hypo_noIMR90")

# Correlation between max/mean hypo, age
cor.test(rmsk_TEother_age_subfamily$JC_distance_mean,rmsk_TEother_age_subfamily$Max_hypo)
by(rmsk_TEother_age_subfamily,rmsk_TEother_age_subfamily$class_update,function(x) cor.test(x$JC_distance_mean,x$Max_hypo))
cor.test(rmsk_TEother_age_subfamily$JC_distance_mean,rmsk_TEother_age_subfamily$Max_hypo_noIMR90)
by(rmsk_TEother_age_subfamily,rmsk_TEother_age_subfamily$class_update,function(x) cor.test(x$JC_distance_mean,x$Max_hypo_noIMR90))
cor.test(rmsk_TEother_age_subfamily$JC_distance_mean,rmsk_TEother_age_subfamily$Mean_hypo)
by(rmsk_TEother_age_subfamily,rmsk_TEother_age_subfamily$class_update,function(x) cor.test(x$JC_distance_mean,x$Mean_hypo))
cor.test(rmsk_TEother_age_subfamily$JC_distance_mean,rmsk_TEother_age_subfamily$Mean_hypo_noIMR90)
by(rmsk_TEother_age_subfamily,rmsk_TEother_age_subfamily$class_update,function(x) cor.test(x$JC_distance_mean,x$Mean_hypo_noIMR90))

# Adding H3K27ac
rmsk_TEother_age_subfamily = merge(rmsk_TEother_age_subfamily,subfamily_H3K27ac_sample_counts[,3:4],by.x="subfamily",by.y="Subfamily",all.x=TRUE)
colnames(rmsk_TEother_age_subfamily)[28] = "H3K27ac"

# Correlation between H3K27ac enrichment, age
cor.test(rmsk_TEother_age_subfamily$JC_distance_mean,rmsk_TEother_age_subfamily$H3K27ac) 
by(rmsk_TEother_age_subfamily,rmsk_TEother_age_subfamily$class_update,function(x) unlist(cor.test(x$JC_distance_mean,x$H3K27ac))[c("p.value","estimate.cor")])
wilcox.test(rmsk_TEother_age_subfamily[which(rmsk_TEother_age_subfamily$H3K27ac == 0),]$JC_distance_mean,rmsk_TEother_age_subfamily[which(rmsk_TEother_age_subfamily$H3K27ac > 0),]$JC_distance_mean)

# Correlation between hypomethylated CpG enrichment, age
rmsk_TEother_age_subfamily = merge(rmsk_TEother_age_subfamily,subfamily_hypo_sample_counts[,3:4],by=c("subfamily"))
colnames(rmsk_TEother_age_subfamily)[33] = "CpG_Hypo"
rmsk_TEother_age_subfamily[which(rmsk_TEother_age_subfamily$JC_distance_mean < mean(rmsk_TEother_age_subfamily$JC_distance_mean)-2*sd(rmsk_TEother_age_subfamily$JC_distance_mean)),]