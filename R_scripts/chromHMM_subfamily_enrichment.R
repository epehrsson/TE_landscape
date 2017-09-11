# Subfamily chromHMM enrichment
# See 4/27/2016, 4/28/2016, 5/4/2016, 5/24/2016, 6/29/2016, 7/6/2016, 7/7/2016, 7/8/2016, 7/10/2016, 7/11/2016, 8/28/2016, 8/29/2016, 9/8/2016, 9/9/2016, 9/17/2016, 9/18/2016, 9/27/2016, 9/29/2016, 11/9/2016, 11/10/2016, 11/15/2016, 11/16/2016, 11/22/2016, 11/27/2016, 
# 1/13/2017, 1/19/2017, 1/20/2017, 1/26/2017, 2/6/2017, 2/9/2017, 2/10/2017, 2/13/2017, 2/14/2017, 2/16/2017, 2/21/2017, 3/16/2017, 5/15/2017, 5/16/2017, 5/23/2017, 5/24/2017, 5/29/2017, 5/30/2017, 6/15/2017

library(plyr)
source("R_scripts/TE_subfamily_stats.R")

# Number of bp in each subfamily x sample x state
subfamily_state_sample = read.table("chromHMM/subfamily/subfamily_state_sample.txt",sep='\t')
colnames(subfamily_state_sample) = c("subfamily","State","Sample","Length_ijk")

# Adding subfamily x state x sample to matrix
subfamily_state_sample_expand = expand.grid(subfamily = levels(subfamily_state_sample$subfamily),Sample = levels(subfamily_state_sample$Sample),State = levels(subfamily_state_sample$State))
subfamily_state_sample_expand = join(subfamily_state_sample_expand,rmsk_TE_subfamily[,1:3],by=c("subfamily"))
subfamily_state_sample = join(subfamily_state_sample_expand,subfamily_state_sample,by=c("subfamily","Sample","State"),type="left")[,c(1,5,4,2,3,6)]
subfamily_state_sample[is.na(subfamily_state_sample$Length_ijk),]$Length_ijk = 0
rm(subfamily_state_sample_expand)

# Number of bp in each subfamily by sample
subfamily_sample = aggregate(data=subfamily_state_sample,Length_ijk ~ class_update+family+subfamily+Sample,FUN=sum)
colnames(subfamily_sample)[5] = "Length_ik"
subfamily_state_sample = join(subfamily_state_sample,subfamily_sample,by=c("subfamily","family","class_update","Sample"),type="left")
rm(subfamily_sample)

# Number of bp in each state by sample
mnemonics_states_genome = melt(as.matrix(t(read.table("chromHMM/genome/mnemonics_state.txt",sep='\t',header=TRUE,row.names=1))))
colnames(mnemonics_states_genome) = c("Sample","State","Length_jk")

subfamily_state_sample = join(subfamily_state_sample,mnemonics_states_genome,by=c("State","Sample"),type="left")
colnames(subfamily_state_sample)[8] = "Length_jk"

# Number of bp in each sample
subfamily_state_sample = join(subfamily_state_sample,mnemonics_states_genome[which(mnemonics_states_genome$State == "Total"),],by=c("Sample"),type="left")
colnames(subfamily_state_sample)[10] = "Length_k"
subfamily_state_sample = subfamily_state_sample[,c(1:8,10)]
rm(mnemonics_states_genome)

# Enrichment for subfamily x state x sample
subfamily_state_sample$Enrichment = log2((subfamily_state_sample$Length_ijk/subfamily_state_sample$Length_ik)/(subfamily_state_sample$Length_jk/subfamily_state_sample$Length_k))

# Proportion of state in subfamily in sample
subfamily_state_sample$Length_percent_jk = subfamily_state_sample$Length_ijk/subfamily_state_sample$Length_jk

# Number of subfamily members in state
subfamily_state_sample_members = rbind(read.table("chromHMM/subfamily/subfamily_state_sample_old.txt",sep='\t'),read.table("chromHMM/subfamily/other_subfamily_state_sample.txt",sep='\t'))
colnames(subfamily_state_sample_members) = c("subfamily","class","family","State","Sample","Members")

# Proportion of members in state
subfamily_state_sample_members$Percent = subfamily_state_sample_members$Members/rmsk_TE_subfamily[match(subfamily_state_sample_members$subfamily,rmsk_TE_subfamily$subfamily),]$Count
subfamily_state_sample_members[which(subfamily_state_sample_members$Sample %in% c("E116","E117","E123","E124","E126","E127")),]$Percent = subfamily_state_sample_members[which(subfamily_state_sample_members$Sample %in% c("E116","E117","E123","E124","E126","E127")),]$Members/rmsk_TE_subfamily[match(subfamily_state_sample_members[which(subfamily_state_sample_members$Sample %in% c("E116","E117","E123","E124","E126","E127")),]$subfamily,rmsk_TE_subfamily$subfamily),]$Count_noY

subfamily_state_sample = merge(subfamily_state_sample,subfamily_state_sample_members,by=c("subfamily","family","State","Sample"),all.x=TRUE)
subfamily_state_sample[which(is.na(subfamily_state_sample$Members)),]$Members = 0
subfamily_state_sample[which(is.na(subfamily_state_sample$Percent)),]$Percent = 0
rm(subfamily_state_sample_members)

# Add metadata for samples
subfamily_state_sample = join(subfamily_state_sample,metadata[,c(1,4:9)])
subfamily_state_sample = subfamily_state_sample[,c(1:2,5,3,4,15:20,6:10,11,13:14)]

# Number of enrichments per subfamily x state
subfamily_state_sample_counts = ddply(subfamily_state_sample,.(class_update,family,subfamily,State),function(x) sum(x$Enrichment > 1.5 & x$Length_ijk >= 600 & x$Length_ik > 5000))

# Number of >1% per subfamily x state
subfamily_state_sample_counts_pc = ddply(subfamily_state_sample,.(class_update,family,subfamily,State),function(x) sum(x$Length_percent_jk > 0.01 & x$Length_ijk >= 600))

# Matrix of subfamily x sample for each state
subfamily_sample_byState = by(subfamily_state_sample,subfamily_state_sample$State,function(x) dcast(x,subfamily ~ Sample,value.var="Enrichment")) 