# Enrichment by subfamily x state x sample with threshold
subfamily_state_sample_filter = subfamily_state_sample[which(subfamily_state_sample$Length_ijk >= 600 & subfamily_state_sample$Length_ik > 5000),]

# Proportion of sample grouping where subfamily is enriched
enrichment_proportion_type = aggregate(data=subfamily_state_sample_filter,Enrichment~Subfamily+State+Type,function(x) sum(x > 1.5))
enrichment_proportion_type$Metadata = rep("Type",dim(enrichment_proportion_type)[1])
enrichment_proportion_test = rbind(enrichment_proportion_age,enrichment_proportion_anatomy,enrichment_proportion_cancer,enrichment_proportion_germline,enrichment_proportion_group,enrichment_proportion_type)
enrichment_proportion_test$Proportion = apply(enrichment_proportion_test,1,function(x) as.numeric(x[4])/length(EID_metadata[which(EID_metadata[,x[5]] == x[3]),x[5]]))

# Proportion of sample grouping where subfamily is >1% of state
enrichment_proportion_cancer = aggregate(data=subfamily_state_sample[which(subfamily_state_sample$Length_ijk >= 600),],Length_percent_jk~Subfamily+State+Cancer,function(x) sum(x > 0.01))
enrichment_proportion_percent = rbind(enrichment_proportion_age,enrichment_proportion_anatomy,enrichment_proportion_cancer,enrichment_proportion_germline,enrichment_proportion_group,enrichment_proportion_type)
enrichment_proportion_percent$Proportion = apply(enrichment_proportion_percent,1,function(x) as.numeric(x[4])/length(EID_metadata[which(EID_metadata[,x[5]] == x[3]),x[5]]))

# With thresholds 
subfamily_CpG_meth[which(subfamily_CpG_meth$Enrichment > 1.5 & subfamily_CpG_meth$CpGs >= 25 & subfamily_CpG_meth$Hypomethylated >= 6),]

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

# Filtered with thresholds
subfamily_DNase_sample_filter = subfamily_DNase_sample[which(subfamily_DNase_sample$Length_ijk >= 600 & subfamily_DNase_sample$Length_ik > 5000),]

# Enrichment in Roadmap sample groups
enrichment_proportion_test$Proportion = apply(enrichment_proportion_test,1,function(x) as.numeric(x[3])/length(EID_metadata[which(EID_metadata[,x[4]] == x[2]),x[4]]))

# Filtered with thresholds
subfamily_H3K27ac_sample_filter = subfamily_H3K27ac_sample[which(subfamily_H3K27ac_sample$Length_ijk >= 600 & subfamily_H3K27ac_sample$Length_ik > 5000),]