# Subfamily chromHMM enrichment
# See 4/27/2016, 4/28/2016, 5/4/2016, 5/24/2016, 6/29/2016, 7/6/2016, 7/7/2016, 7/8/2016, 7/10/2016, 7/11/2016, 8/28/2016, 8/29/2016, 9/8/2016, 9/9/2016, 9/17/2016, 9/18/2016, 9/27/2016, 9/29/2016, 11/9/2016, 11/10/2016, 11/15/2016, 11/16/2016, 11/22/2016, 11/27/2016, 
# 1/13/2017, 1/19/2017, 1/20/2017, 1/26/2017, 2/6/2017, 2/9/2017, 2/10/2017, 2/13/2017, 2/14/2017, 2/16/2017, 2/21/2017, 3/16/2017, 5/15/2017, 5/16/2017, 5/23/2017, 5/24/2017, 5/29/2017, 5/30/2017, 6/15/2017

# Number of bp in each subfamily x sample x state
subfamily_state_sample = read.table("TE_subfamilies/subfamily_state_sample.txt",sep='\t')
colnames(subfamily_state_sample) = c("Subfamily","State","Sample","Length_ijk")

# Adding subfamily x state x sample to matrix
subfamily_state_sample_expand = expand.grid(Subfamily = levels(subfamily_state_sample$Subfamily),Sample = levels(subfamily_state_sample$Sample),State = levels(subfamily_state_sample$State))
subfamily_state_sample_expand = join(subfamily_state_sample_expand,rbind(rmsk_TE_stats_subfamily,rmsk_other_stats_subfamily)[,1:3],by=c("Subfamily"))
subfamily_state_sample = join(subfamily_state_sample_expand,subfamily_state_sample,by=c("Subfamily","Sample","State"),type="left")[,c(1,5,4,2,3,6)]
subfamily_state_sample[is.na(subfamily_state_sample$Length_ijk),]$Length_ijk = 0

# Number of bp in each subfamily by sample
subfamily_sample = aggregate(data=subfamily_state_sample,Length_ijk ~ Class+Family+Subfamily+Sample,FUN=sum)
colnames(subfamily_sample)[5] = "Length_ik"
subfamily_state_sample = join(subfamily_state_sample,subfamily_sample,by=c("Subfamily","Family","Class","Sample"),type="left")

# Number of bp in each state by sample
subfamily_state_sample = join(subfamily_state_sample,mnemonics_states_long,by=c("State","Sample"),type="left")
colnames(subfamily_state_sample)[8] = "Length_jk"

# Number of bp in each sample
subfamily_state_sample = join(subfamily_state_sample,mnemonics_states_long[which(mnemonics_states_long$State == "Total"),],by=c("Sample"),type="left")
colnames(subfamily_state_sample)[10] = "Length_k"
subfamily_state_sample = subfamily_state_sample[,c(1:8,10)]

# Enrichment for subfamily x state x sample
subfamily_state_sample$Enrichment = log2((subfamily_state_sample$Length_ijk/subfamily_state_sample$Length_ik)/(subfamily_state_sample$Length_jk/subfamily_state_sample$Length_k))

# Adding proportion of subfamily members in state
subfamily_state_sample_members = read.table("TE_landscape/subfamily_state_sample.txt",sep='\t')
colnames(subfamily_state_sample_members) = c("Subfamily","Class","Family","State","Sample","Members")

# Expanded to all combinations
subfamily_state_sample_members = join(subfamily_state_sample_expand,subfamily_state_sample_members,by=c("Subfamily","Family","Class","Sample","State"),type="left")[,c(1,5,4,3,2,6)]
subfamily_state_sample_members[is.na(subfamily_state_sample_members$Members),]$Members = 0

# Addition proportion of members
subfamily_state_sample_members$Percent = subfamily_state_sample_members$Members/rmsk_TE_stats_subfamily[match(subfamily_state_sample_members$Subfamily,rmsk_TE_stats_subfamily$Subfamily),]$Count

# Addition of proportion of members, Y chromosome
subfamily_state_sample_members[which(subfamily_state_sample_members$Sample %in% c("E116","E117","E123","E124","E126","E127")),]$Percent = subfamily_state_sample_members[which(subfamily_state_sample_members$Sample %in% c("E116","E117","E123","E124","E126","E127")),]$Members/rmsk_TE_stats_subfamily[match(subfamily_state_sample_members[which(subfamily_state_sample_members$Sample %in% c("E116","E117","E123","E124","E126","E127")),]$Subfamily,rmsk_TE_stats_subfamily$Subfamily),]$Count_noY

subfamily_state_sample_members_other = read.table("other/other_subfamily_state_sample.txt",sep='\t')
colnames(subfamily_state_sample_members_other) = c("Subfamily","Class","Family","State","Sample","Members")
subfamily_state_sample_members_other = join(join(expand.grid(Subfamily = levels(subfamily_state_sample_members_other$Subfamily),Sample = levels(subfamily_state_sample_members_other$Sample),State = levels(subfamily_state_sample_members_other$State)),rmsk_other_stats_subfamily[,1:3],by=c("Subfamily")),subfamily_state_sample_members_other,by=c("Subfamily","Family","Class","Sample","State"),type="left")
subfamily_state_sample_members_other[is.na(subfamily_state_sample_members_other$Members),]$Members = 0
subfamily_state_sample_members_other$Percent = subfamily_state_sample_members_other$Members/rmsk_other_stats_subfamily[match(subfamily_state_sample_members_other$Subfamily,rmsk_other_stats_subfamily$Subfamily),]$Count
subfamily_state_sample_members_other[which(subfamily_state_sample_members_other$Sample %in% c("E116","E117","E123","E124","E126","E127")),]$Percent = subfamily_state_sample_members_other[which(subfamily_state_sample_members_other$Sample %in% c("E116","E117","E123","E124","E126","E127")),]$Members/rmsk_other_stats_subfamily[match(subfamily_state_sample_members_other[which(subfamily_state_sample_members_other$Sample %in% c("E116","E117","E123","E124","E126","E127")),]$Subfamily,rmsk_other_stats_subfamily$Subfamily),]$Count_noY

subfamily_state_sample = merge(subfamily_state_sample,rbind(subfamily_state_sample_members,subfamily_state_sample_members_other),by=c("Subfamily","Family","Class","State","Sample"))

# Add metadata for samples
subfamily_state_sample = join(subfamily_state_sample,EID_metadata[,c(1,3,5,7:11)])
subfamily_state_sample = subfamily_state_sample[,c(1:5,13:19,6:12)]

# Proportion of state in subfamily in sample
subfamily_state_sample$Length_percent_jk = subfamily_state_sample$Length_ijk/subfamily_state_sample$Length_jk

# Subfamilies with >1% bp of state, not enriched
table(droplevels(subfamily_state_sample[which(subfamily_state_sample$Length_percent_jk > 0.01 & subfamily_state_sample$Enrichment <= 1.5 & subfamily_state_sample$State %in% chromHMM_states[1:7]),c(1,4)]))

# Enrichment by subfamily x state x sample with threshold
subfamily_state_sample_filter = subfamily_state_sample[which(subfamily_state_sample$Length_ijk >= 600 & subfamily_state_sample$Length_ik > 5000),]
subfamily_state_sample_filter$Class_update = subfamily_state_sample_filter$Class
subfamily_state_sample_filter$Class_update = factor(subfamily_state_sample_filter$Class_update,levels=c("DNA","LINE","LTR","SINE","Other","RC","Unconfident"))
subfamily_state_sample_filter[which(subfamily_state_sample_filter$Class %in% c("DNA?","LINE?","SINE?","LTR?","Unknown","Unknown?")),]$Class_update = "Unconfident"

# Average proportion of members in a state for subfamilies enriched in a state
by(subfamily_state_sample_filter[which(subfamily_state_sample_filter$Enrichment > 1.5),],subfamily_state_sample_filter[which(subfamily_state_sample_filter$Enrichment > 1.5),]$State,function(x) mean(x$Percent))

# Number of enrichments per subfamily x state
subfamily_state_sample_counts = ddply(subfamily_state_sample,.(Class,Family,Subfamily,State),function(x) sum(x$Enrichment > 1.5 & x$Length_ijk >= 600 & x$Length_ik > 5000))
subfamily_state_sample_counts$Class_update = subfamily_state_sample_counts$Class
subfamily_state_sample_counts$Class_update = factor(subfamily_state_sample_counts$Class_update,levels=c("DNA","LINE","LTR","SINE","RC","Other","Unconfident"))
subfamily_state_sample_counts[which(subfamily_state_sample_counts$Class %in% c("DNA?","LINE?","LTR?","SINE?","Unknown","Unknown?")),]$Class_update = "Unconfident"

# Proportion of sample grouping where subfamily is enriched
enrichment_proportion_type = aggregate(data=subfamily_state_sample_filter,Enrichment~Subfamily+State+Type,function(x) sum(x > 1.5))
enrichment_proportion_type$Metadata = rep("Type",dim(enrichment_proportion_type)[1])
enrichment_proportion_test = rbind(enrichment_proportion_age,enrichment_proportion_anatomy,enrichment_proportion_cancer,enrichment_proportion_germline,enrichment_proportion_group,enrichment_proportion_type)
enrichment_proportion_test$Proportion = apply(enrichment_proportion_test,1,function(x) as.numeric(x[4])/length(EID_metadata[which(EID_metadata[,x[5]] == x[3]),x[5]]))

# Proportion of sample grouping where subfamily is >1% of state
enrichment_proportion_cancer = aggregate(data=subfamily_state_sample[which(subfamily_state_sample$Length_ijk >= 600),],Length_percent_jk~Subfamily+State+Cancer,function(x) sum(x > 0.01))
enrichment_proportion_percent = rbind(enrichment_proportion_age,enrichment_proportion_anatomy,enrichment_proportion_cancer,enrichment_proportion_germline,enrichment_proportion_group,enrichment_proportion_type)
enrichment_proportion_percent$Proportion = apply(enrichment_proportion_percent,1,function(x) as.numeric(x[4])/length(EID_metadata[which(EID_metadata[,x[5]] == x[3]),x[5]]))

# Plotting instances where subfamily is >1% of state
percent_all = subfamily_state_sample[which(subfamily_state_sample$Length_ijk >= 600 & subfamily_state_sample$Length_percent_jk > 0.01 & subfamily_state_sample$State %in% chromHMM_states[1:7] & subfamily_state_sample$Subfamily %in% c("L2a","L2b","L2c","MIR","MIRb")),]
percent_all$Subfamily_State = paste(percent_all$Subfamily,percent_all$State,sep=":")
percent_all = dcast(percent_all[,c(5,20,21)],Subfamily_State~Sample,value.var = c("Length_percent_jk"))
rownames(percent_all) = percent_all[,1]
percent_all = percent_all[,2:128]
percent_all[is.na(percent_all)] = 0
percent_all[percent_all > 0.01] = 1

# Statistics and plots
# Subfamilies with at least one enrichment by state
subfamily_state_sample_enriched = sapply(chromHMM_states,function(x) as.vector(unique(subfamily_state_sample_filter[which(subfamily_state_sample_filter$State == x & subfamily_state_sample_filter$Enrichment > 1.5),]$Subfamily)))
names(subfamily_state_sample_enriched) = chromHMM_states

# Total number of enriched subfamilies
unique(unlist(subfamily_state_sample_enriched))

# Number of enriched subfamilies by class and state
by(subfamily_state_sample_filter[which(subfamily_state_sample_filter$Enrichment > 1.5),],subfamily_state_sample_filter[which(subfamily_state_sample_filter$Enrichment > 1.5),]$State,function(x) table(unique(x[1:3])$Class))
sapply(enriched_subfamily,function(x) table(rmsk_TE_stats_subfamily[which(rmsk_TE_stats_subfamily$Subfamily %in% x),]$Class))

# Matrix of subfamily x sample for each state
subfamily_sample_byState = by(subfamily_state_sample,subfamily_state_sample$State,function(x) dcast(x,Subfamily ~ Sample,value.var="Enrichment")) 

# Filling out enrichment table
table(subfamily_state_sample_filter[which(subfamily_state_sample_filter$Enrichment > 1.5),]$State) #Enrichments by state
lapply(subfamily_state_sample_enriched,length) #Enriched subfamilies by state
lapply(subfamily_state_sample_enriched,function(x) table(rmsk_TEother_stats_subfamily[match(x,rmsk_TEother_stats_subfamily$Subfamily),]$Class)) #Enriched subfamilies by state and class

# Updating enrichment spreadsheet
TEother_meth_wCpG_subfamily_hypo[which(!(TEother_meth_wCpG_subfamily_hypo$subfamily %in% TEother_meth_exclude) & TEother_meth_wCpG_subfamily_hypo$Samples_hypo_10 > 0),c(1:3,43)]
c(as.vector(TE_meth_wCpG_subfamily_hypo_var_20$subfamily),as.vector(other_meth_wCpG_subfamily_hypo_var_20$subfamily))
table(droplevels(subfamily_state_sample[which(subfamily_state_sample$Length_percent_jk > 0.01  & subfamily_state_sample$State %in% chromHMM_states[1:7] & subfamily_state_sample$Length_ijk >= 600),c(1,4)]))
sort(unlist(lapply(TEother_meth_hypo_prop,function(x) sum(x > 0.01))))
       