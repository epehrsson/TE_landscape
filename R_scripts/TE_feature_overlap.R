# TE overlap with features
# See 8/24/2016, 8/25/2016, 9/20/2016, 9/28/2016, 1/31/2017, 2/6/2017, 2/9/2017, 2/10/2017, 2/25/2017, 2/27/2017, 3/8/2017, 5/15/2017, 6/14/2017, 8/1/2017

# Overlap files
feature_files = c("rmsk_TEother_refseq_promoters_merge.txt","rmsk_TEother_refseq_5UTR_merge.txt","rmsk_TEother_refseq_coding_exon_merge.txt","rmsk_TEother_refseq_3UTR_merge.txt","rmsk_TEother_refseq_exons_merge.txt","rmsk_TEother_refseq_introns_merge.txt","rmsk_TEother_refseq_intergenic.txt")

# Overlap between TEs and Refseq features
features = lapply(feature_files,function(x) read.table(x,sep='\t'))
features = lapply(features, setNames, nm =c("Chromosome","Start","Stop","Subfamily","Class","Family","Strand"))
colnames(features[[1]])[8] = "Promoter"
colnames(features[[2]])[8] = "5UTR"
colnames(features[[3]])[8] = "CDS"
colnames(features[[4]])[8] = "3UTR"
colnames(features[[5]])[8] = "Exon"
colnames(features[[6]])[8] = "Intron"
colnames(features[[7]])[8] = "Intergenic"
features[[8]] = read.table("rmsk_TEother_cpgIsland.txt",sep='\t')
colnames(features[[8]]) = c("Chromosome","Start","Stop","Subfamily","Class","Family","Strand","CGI")

# Overlap between all TEs and Refseq features
test_feature = merge(rmsk_TEother[,c(1:6,8:9)],features[[1]],by=c("Chromosome","Start","Stop","Subfamily","Family","Class","Strand"),all.x=TRUE)
test_feature = merge(test_feature,features[[2]],by=c("Chromosome","Start","Stop","Subfamily","Family","Class","Strand"),all.x=TRUE)
test_feature = merge(test_feature,features[[3]],by=c("Chromosome","Start","Stop","Subfamily","Family","Class","Strand"),all.x=TRUE)
test_feature = merge(test_feature,features[[4]],by=c("Chromosome","Start","Stop","Subfamily","Family","Class","Strand"),all.x=TRUE)
test_feature = merge(test_feature,features[[5]],by=c("Chromosome","Start","Stop","Subfamily","Family","Class","Strand"),all.x=TRUE)
test_feature = merge(test_feature,features[[6]],by=c("Chromosome","Start","Stop","Subfamily","Family","Class","Strand"),all.x=TRUE)
test_feature = merge(test_feature,features[[7]],by=c("Chromosome","Start","Stop","Subfamily","Family","Class","Strand"),all.x=TRUE)
rmsk_TEother_feature = test_feature
rmsk_TEother_feature = merge(rmsk_TEother_feature,features[[8]],by=c("Chromosome","Start","Stop","Subfamily","Family","Class","Strand"),all.x=TRUE)
colnames(rmsk_TEother_feature)[1:8] = colnames(rmsk_TEother_age)[c(1:7,34)]
rmsk_TEother_feature$Feature_count = apply(rmsk_TEother_feature[,9:15],1,function(x) 7-sum(is.na(x)))

# Proportion of class overlapping each feature
rmsk_TEother_feature_class = ddply(rmsk_TEother_feature,"Class_update",function(y) apply(y[,9:16],2,function(x) length(na.omit(x))/length(x)))
rmsk_TEother_feature_class$Class_update = factor(rmsk_TEother_feature_class$Class_update,levels=c("All",levels(rmsk_TEother_feature_class$Class_update)))
rmsk_TEother_feature_class = rbind(rmsk_TEother_feature_class,c("All",apply(rmsk_TEother_feature[,9:16],2,function(x) length(na.omit(x))/4430788)))
rownames(rmsk_TEother_feature_class) = rmsk_TEother_feature_class[,1]
rmsk_TEother_feature_class = rmsk_TEother_feature_class[,2:9]

# Enrichment of TE classes overlapping features (LOR)
t(apply(rmsk_TEother_feature_class,1,function(x) log2(as.numeric(x)/as.numeric(rmsk_TEother_feature_class[8,]))))

# Proportion of subfamily in each genic feature
rmsk_TEother_feature_subfamily = ddply(rmsk_TEother_feature,"subfamily",function(y) apply(y[,9:15],2,function(x) dim(y[which(!(is.na(x))),])[1]/dim(y)[1]))
rmsk_TEother_feature_subfamily = merge(rmsk_TEother_feature_subfamily,rmsk_TEother_stats_subfamily[,c(1:4,8)],by.x="subfamily",by.y="Subfamily")
rmsk_TEother_feature_subfamily =rmsk_TEother_feature_subfamily[,c(1,9:12,2:8)]

# Subfamilies with much higher proportion in feature
apply(rmsk_TEother_feature_subfamily[,6:12],2,function(x) rmsk_TEother_feature_subfamily[which(x > mean(x)+2*sd(x) & rmsk_TEother_feature_subfamily$Count*x > 2),])

# TEs overlapping the hg19 blacklist
TE_blacklist = read.table("TE_blacklist.bed",sep='\t')
colnames(TE_blacklist) = c("Chromosome","Start","Stop","Subfamily","Class","Family","Strand","Chromosome_blacklist","Start_blacklist","Stop_blacklist","Overlap")

