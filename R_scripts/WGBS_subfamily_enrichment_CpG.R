# Enrichment of hypomethylated CpGs in a subfamily
# See 2/6/2017, 7/21/2017, 7/23/2017, 7/26/2017, 8/2/2017

# Proportion of subfamily CpGs in methylation state by sample
subfamily_CpG_meth = read.table("WGBS/subfamily_CpG_Meth_states.txt",sep='\t')
subfamily_CpG_meth$V1 = mapvalues(subfamily_CpG_meth$V1,seq(4,40,1),as.vector(EID_metadata_meth$Sample))
subfamily_CpG_meth = subfamily_CpG_meth[,c(1,6,2:5)]
colnames(subfamily_CpG_meth) = c("Sample","subfamily",meth_states)
subfamily_CpG_meth[is.na(subfamily_CpG_meth)] = 0

subfamily_CpG_meth = merge(subfamily_CpG_meth,TE_subfamily_CpG_count[,1:4],by.x=c("subfamily"),by.y=c("Subfamily"),all.x=TRUE)

# Proportion of all CpGs hypomethylated
subfamily_CpG_meth$Hypomethylated_all = apply(subfamily_CpG_meth,1,function(x) all_CpG_meth[x[2],]$Hypomethylated)

# Enrichment of hypomethylated CpGs in sample x subfamily
subfamily_CpG_meth$Enrichment = log2((subfamily_CpG_meth$Hypomethylated/subfamily_CpG_meth$CpGs)/(subfamily_CpG_meth$Hypomethylated_all/56434896))

# With thresholds 
subfamily_CpG_meth[which(subfamily_CpG_meth$Enrichment > 1.5 & subfamily_CpG_meth$CpGs >= 25 & subfamily_CpG_meth$Hypomethylated >= 6),]

# Adding Unconfident class
subfamily_CpG_meth$class_update = subfamily_CpG_meth$Class
subfamily_CpG_meth$class_update = factor(subfamily_CpG_meth$class_update,levels=c("DNA","LINE","LTR","SINE","RC","Other","Unconfident"))
subfamily_CpG_meth[which(subfamily_CpG_meth$Class %in% c("DNA?","LINE?","LTR?","SINE?","Unknown","Unknown?")),]$class_update = "Unconfident"

# Proprotion of all hypomethylated CpGs in subfamily
subfamily_CpG_meth$Hypo_percent = subfamily_CpG_meth$Hypomethylated/subfamily_CpG_meth$Hypomethylated_all

# Number of hypomethylation enrichments per subfamily
subfamily_hypo_sample_counts = ddply(subfamily_CpG_meth,.(Class,Family,subfamily),function(x) sum(x$Enrichment > 1.5 & x$CpGs >= 25 & x$Hypomethylated >= 6))

# Adding metadata
subfamily_CpG_meth = merge(subfamily_CpG_meth,EID_metadata[,c(1,5,7:11)],by=c("Sample"),all.x=TRUE)

# Finding proportion of sample groupings enriched
enrichment_proportion_age = aggregate(data=subfamily_CpG_meth[which(subfamily_CpG_meth$CpGs >= 25 & subfamily_CpG_meth$Hypomethylated >= 6),],Enrichment~subfamily+Age,function(x) sum(x > 1.5))
enrichment_proportion_anatomy = aggregate(data=subfamily_CpG_meth[which(subfamily_CpG_meth$CpGs >= 25 & subfamily_CpG_meth$Hypomethylated >= 6),],Enrichment~subfamily+Anatomy,function(x) sum(x > 1.5))
enrichment_proportion_germline = aggregate(data=subfamily_CpG_meth[which(subfamily_CpG_meth$CpGs >= 25 & subfamily_CpG_meth$Hypomethylated >= 6),],Enrichment~subfamily+Germline,function(x) sum(x > 1.5))
enrichment_proportion_group = aggregate(data=subfamily_CpG_meth[which(subfamily_CpG_meth$CpGs >= 25 & subfamily_CpG_meth$Hypomethylated >= 6),],Enrichment~subfamily+Group,function(x) sum(x > 1.5))
enrichment_proportion_type = aggregate(data=subfamily_CpG_meth[which(subfamily_CpG_meth$CpGs >= 25 & subfamily_CpG_meth$Hypomethylated >= 6),],Enrichment~subfamily+Type,function(x) sum(x > 1.5))
enrichment_proportion_age$Metadata = rep("Age",dim(enrichment_proportion_age)[1])
enrichment_proportion_anatomy$Metadata = rep("Anatomy",dim(enrichment_proportion_anatomy)[1])
enrichment_proportion_germline$Metadata = rep("Germline",dim(enrichment_proportion_germline)[1])
enrichment_proportion_group$Metadata = rep("Group",dim(enrichment_proportion_group)[1])
enrichment_proportion_type$Metadata = rep("Type",dim(enrichment_proportion_type)[1])
colnames(enrichment_proportion_age)[2] = "Category"
colnames(enrichment_proportion_anatomy)[2] = "Category"
colnames(enrichment_proportion_germline)[2] = "Category"
colnames(enrichment_proportion_group)[2] = "Category"
colnames(enrichment_proportion_type)[2] = "Category"
enrichment_proportion_test = rbind(enrichment_proportion_age,enrichment_proportion_anatomy,enrichment_proportion_germline,enrichment_proportion_group,enrichment_proportion_type)
enrichment_proportion_test$Proportion = apply(enrichment_proportion_test,1,function(x) as.numeric(x[3])/length(EID_metadata_meth[which(EID_metadata_meth[,x[4]] == x[2]),x[4]]))
enrichment_proportion_test[which(enrichment_proportion_test$Proportion > 0.5 & enrichment_proportion_test$Enrichment > 1),] #95 entries