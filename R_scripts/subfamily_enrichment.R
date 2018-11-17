# Subfamily enrichment
# See 4/27/2016, 4/28/2016, 5/4/2016, 5/24/2016, 6/29/2016, 7/6/2016, 7/7/2016, 7/8/2016, 7/10/2016, 7/11/2016, 8/28/2016, 8/29/2016, 9/8/2016, 9/9/2016, 9/17/2016, 9/18/2016, 9/27/2016, 9/29/2016, 11/9/2016, 11/10/2016, 11/15/2016, 11/16/2016, 11/22/2016, 11/27/2016, 
# 1/13/2017, 1/19/2017, 1/20/2017, 1/26/2017, 2/6/2017, 2/9/2017, 2/10/2017, 2/13/2017, 2/14/2017, 2/16/2017, 2/21/2017, 3/16/2017, 5/15/2017, 5/16/2017, 5/23/2017, 5/24/2017, 5/29/2017, 5/30/2017, 6/5/17, 6/6/17, 6/7/17, 6/9/17, 6/15/2017, 7/4/2017, 7/21/2017, 7/23/2017, 7/24/2017, 7/26/2017, 8/2/2017

# Updated with peak summit overlap 5/2/18

# chromHMM/DNase/H3K27ac

## Number of unique chromHMM blocks/peaks per sample per subfamily

# chromHMM
subfamily_chromHMM_sample = read.table("chromHMM/subfamily/subfamily_state_sample.txt",sep='\t')
colnames(subfamily_chromHMM_sample) = c("subfamily","State","Sample","Length_ijk")

subfamily_chromHMM_sample_expand = expand.grid(subfamily = rmsk_TE_subfamily$subfamily,Sample = metadata$Sample,State = chromHMM_states)
subfamily_chromHMM_sample = merge(subfamily_chromHMM_sample,subfamily_chromHMM_sample_expand,by=c("subfamily","Sample","State"),all.y=TRUE)
subfamily_chromHMM_sample[is.na(subfamily_chromHMM_sample)] = 0
rm(subfamily_chromHMM_sample_expand)

# DNase
subfamily_DNase_sample = read.table("DNase/subfamily_DNase_sample_summit.txt",sep='\t')
colnames(subfamily_DNase_sample) = c("Sample","subfamily","Length_ijk")
subfamily_DNase_sample$State = rep("DNase",dim(subfamily_DNase_sample)[1])

# H3K27ac
subfamily_H3K27ac_sample = read.table("H3K27ac/subfamily_H3K27ac_sample_summit.txt",sep='\t')
colnames(subfamily_H3K27ac_sample) = c("Sample","subfamily","Length_ijk")
subfamily_H3K27ac_sample$State = rep("H3K27ac",dim(subfamily_H3K27ac_sample)[1])

# Combine
subfamily_state_sample = rbind(subfamily_chromHMM_sample,subfamily_DNase_sample,subfamily_H3K27ac_sample)
rm(list=c("subfamily_chromHMM_sample","subfamily_DNase_sample","subfamily_H3K27ac_sample"))

## Number of unique chromHMM blocks/peaks per sample

# chromHMM (total length)
chromHMM_state_sample = mnemonics_states_genome[,1:3]
colnames(chromHMM_state_sample)[3] = "Length_jk"

# DNase
DNase_peaks_sample = DNase_stats[,c("Sample","Peaks")]
DNase_peaks_sample$State = rep("DNase",dim(DNase_peaks_sample)[1])
colnames(DNase_peaks_sample)[2] = "Length_jk"

# H3K27ac
H3K27ac_peaks_sample = H3K27ac_stats[,c("Sample","Peaks")]
H3K27ac_peaks_sample$State = rep("H3K27ac",dim(H3K27ac_peaks_sample)[1])
colnames(H3K27ac_peaks_sample)[2] = "Length_jk"

# Combine 
peaks_sample = rbind(chromHMM_state_sample,DNase_peaks_sample,H3K27ac_peaks_sample)
rm(list=c("chromHMM_state_sample","DNase_peaks_sample","H3K27ac_peaks_sample"))

## Join
subfamily_state_sample = merge(subfamily_state_sample,peaks_sample,by=c("Sample","State"),all.x=TRUE)
rm(peaks_sample)

## Length of subfamily
subfamily_state_sample$Length_ik = ifelse(metadata[match(subfamily_state_sample$Sample,metadata$Sample),]$chrY == "Yes",
                                          rmsk_TE_subfamily[match(subfamily_state_sample$subfamily,rmsk_TE_subfamily$subfamily),]$Total_length,
                                          rmsk_TE_subfamily[match(subfamily_state_sample$subfamily,rmsk_TE_subfamily$subfamily),]$Total_length_noY)

## Length of genome
subfamily_state_sample$Length_k = ifelse(metadata[match(subfamily_state_sample$Sample,metadata$Sample),]$chrY == "Yes",GENOME_WIDTH,GENOME_WIDTH_noY)

# LOR Enrichment for subfamily x state x sample
subfamily_state_sample$Enrichment = log2(1e-20+((subfamily_state_sample$Length_ijk/subfamily_state_sample$Length_jk)/(subfamily_state_sample$Length_ik/subfamily_state_sample$Length_k)))

# Add class and family
subfamily_state_sample = merge(subfamily_state_sample,rmsk_TE_subfamily[,c("subfamily","family","class_update")],by="subfamily",all.x=TRUE)

# Add metadata for samples
subfamily_state_sample = merge(subfamily_state_sample,metadata[,c(1,4:9)],by=c("Sample"),all.x=TRUE)

# Proportion of peaks in subfamily in sample
subfamily_state_sample$Length_percent_jk = subfamily_state_sample$Length_ijk/subfamily_state_sample$Length_jk

## Number of subfamily members in state

# chromHMM
subfamily_chromHMM_sample_members = read.table("chromHMM/subfamily/subfamily_state_sample_summit.txt",sep='\t')
colnames(subfamily_chromHMM_sample_members) = c("subfamily","Sample","State","Members")

# DNase
TE_DNase_peaks_members = melt(aggregate(data=TE_DNase_peaks[,c(4,8:60)],.~subfamily,function(x) sum(as.numeric(x) > 0)),id.vars=c("subfamily"))
colnames(TE_DNase_peaks_members)[2:3] = c("Sample","Members")
TE_DNase_peaks_members$State = rep("DNase",dim(TE_DNase_peaks_members)[1])

# H3K27ac
TE_H3K27ac_peaks_members = melt(aggregate(data=TE_H3K27ac_peaks[,c(4,8:105)],.~subfamily,function(x) sum(as.numeric(x) > 0)),id.vars=c("subfamily"))
colnames(TE_H3K27ac_peaks_members)[2:3] = c("Sample","Members")
TE_H3K27ac_peaks_members$State = rep("H3K27ac",dim(TE_H3K27ac_peaks_members)[1])

# Combine
subfamily_state_sample_members = rbind(subfamily_chromHMM_sample_members,TE_DNase_peaks_members,TE_H3K27ac_peaks_members)
rm(list=c("subfamily_chromHMM_sample_members","TE_DNase_peaks_members","TE_H3K27ac_peaks_members"))

# Join
subfamily_state_sample = merge(subfamily_state_sample,subfamily_state_sample_members,by=c("subfamily","State","Sample"),all.x=TRUE)
subfamily_state_sample[which(is.na(subfamily_state_sample$Members)),]$Members = 0
rm(subfamily_state_sample_members)

## Proportion of members in state
subfamily_state_sample$Count = ifelse(metadata[match(subfamily_state_sample$Sample,metadata$Sample),]$chrY == "Yes",
                                                                                       rmsk_TE_subfamily[match(subfamily_state_sample$subfamily,rmsk_TE_subfamily$subfamily),]$Count,
                                                                                       rmsk_TE_subfamily[match(subfamily_state_sample$subfamily,rmsk_TE_subfamily$subfamily),]$Count_noY)
subfamily_state_sample$Percent = subfamily_state_sample$Members/subfamily_state_sample$Count

# WGBS (CpGs)
# Proportion of subfamily CpGs in methylation state by sample
subfamily_CpG_meth = read.table("WGBS/subfamily_CpG_Meth_states.txt",sep='\t')
subfamily_CpG_meth$V1 = mapvalues(subfamily_CpG_meth$V1,seq(4,40,1),as.vector(metadata[which(!is.na(metadata$WGBS)),]$Sample))
colnames(subfamily_CpG_meth) = c("Sample",meth_states,"subfamily")
subfamily_CpG_meth[is.na(subfamily_CpG_meth)] = 0
subfamily_CpG_meth[,meth_states] = subfamily_CpG_meth[,meth_states]/2
subfamily_CpG_meth = melt(subfamily_CpG_meth,id.vars=c("Sample","subfamily"))
colnames(subfamily_CpG_meth)[3:4] = c("State","CpG_ijk")

# CpGs per subfamily
subfamily_CpG_meth = merge(subfamily_CpG_meth,rmsk_TE_subfamily[,c("subfamily","family","class_update","CpGs")],by=c("subfamily"),all.x=TRUE)
colnames(subfamily_CpG_meth)[7] = "CpG_ik"

# Proportion of all CpGs in methylation state by sample
subfamily_CpG_meth$CpG_jk = apply(subfamily_CpG_meth,1,function(x) all_CpG_meth[x[2],x[3]])

# Enrichment of state CpGs in sample x subfamily
subfamily_CpG_meth$Enrichment = log2(1e-20+((subfamily_CpG_meth$CpG_ijk/subfamily_CpG_meth$CpG_ik)/(subfamily_CpG_meth$CpG_jk/ALL_CPGS)))

# Proportion of all state CpGs in subfamily
subfamily_CpG_meth$CpG_ijk_jk = subfamily_CpG_meth$CpG_ijk/subfamily_CpG_meth$CpG_jk

# Adding metadata
subfamily_CpG_meth = merge(subfamily_CpG_meth,metadata[,c(1,4:9)],by=c("Sample"),all.x=TRUE)

# Members in state (average methylation)
subfamily_CpG_meth = merge(subfamily_CpG_meth,TE_meth_subfamily[,c("subfamily","State","Sample","Members")],by=c("subfamily","Sample","State"))

# Divide by number of members with CpGs in subfamily
subfamily_CpG_meth$Count = rmsk_TE_subfamily[match(subfamily_CpG_meth$subfamily,rmsk_TE_subfamily$subfamily),]$Count_CpGs
subfamily_CpG_meth$Percent = subfamily_CpG_meth$Members/subfamily_CpG_meth$Count

# Combine matrices (no filtering)
columns = c("subfamily","family","class_update","State","Sample",sample_categories,enrichment_names[c(1:2,5:9)])
subfamily_state_sample_combined = rbind(subfamily_state_sample[,columns],
                                        rename(subfamily_CpG_meth,c("CpG_ijk"="Length_ijk","CpG_ik"="Length_ik","CpG_ijk_jk"="Length_percent_jk"))[,columns])

# Combine filtered matrices
subfamily_state_sample_filter = subfamily_state_sample_combined[which(subfamily_state_sample_combined$Members > THRESHOLD_IJK_MEMBER & subfamily_state_sample_combined$Count > THRESHOLD_IK_MEMBER),]
subfamily_state_sample_filter$State = factor(subfamily_state_sample_filter$State,levels=states[1:21])

# Number of enrichments per subfamily x state
subfamily_state_sample_counts = ddply(subfamily_state_sample_filter,.(class_update,family,subfamily,State),function(x) sum(x$Enrichment > THRESHOLD_LOR))
subfamily_state_expand = expand.grid(subfamily = levels(subfamily_state_sample$subfamily),State = levels(subfamily_state_sample_combined$State))
subfamily_state_expand = join(subfamily_state_expand,rmsk_TE_subfamily[,c("subfamily","family","class_update")],by=c("subfamily"))
subfamily_state_sample_counts = join(subfamily_state_expand,subfamily_state_sample_counts,by=c("subfamily","family","class_update","State"),type="left")[,c(1,3,4,2,5)]
rm(subfamily_state_expand)
subfamily_state_sample_counts[is.na(subfamily_state_sample_counts)] = 0
subfamily_state_sample_counts$State = factor(subfamily_state_sample_counts$State,levels=states[1:21])

subfamily_state_sample_counts$Sample.Proportion = ifelse(subfamily_state_sample_counts$State %in% chromHMM_states,subfamily_state_sample_counts$V1/sample_counts["All","chromHMM"],
                                                         ifelse(subfamily_state_sample_counts$State %in% meth_states,subfamily_state_sample_counts$V1/sample_counts["All","WGBS"],
                                                                ifelse(subfamily_state_sample_counts$State == "DNase",subfamily_state_sample_counts$V1/sample_counts["All","DNase"],
                                                                       ifelse(subfamily_state_sample_counts$State == "H3K27ac",subfamily_state_sample_counts$V1/sample_counts["All","H3K27ac"],"NA"))))
subfamily_state_sample_counts$Sample.Proportion = as.numeric(subfamily_state_sample_counts$Sample.Proportion)

subfamily_state_sample_counts_combine = merge(aggregate(data=subfamily_state_sample_counts,V1~State,function(x) sum(x > 0)),
                                              dcast(aggregate(data=subfamily_state_sample_counts,V1~State+class_update,function(x) sum(x > 0)),State~class_update,value.var = "V1"),
                                              by=c("State"))
subfamily_state_sample_counts_combine = subfamily_state_sample_counts_combine[match(levels(subfamily_state_sample_counts$State),subfamily_state_sample_counts_combine$State),]

# Number of >1% per subfamily x state
subfamily_state_sample_counts_pc = ddply(subfamily_state_sample_combined,.(class_update,family,subfamily,State),function(x) sum(x$Length_percent_jk > THRESHOLD_PC))
subfamily_state_sample_counts_pc_combine = merge(aggregate(data=subfamily_state_sample_counts_pc,V1~State,function(x) sum(x > 0)),
                                              dcast(aggregate(data=subfamily_state_sample_counts_pc,V1~State+class_update,function(x) sum(x > 0)),State~class_update,value.var = "V1"),
                                              by=c("State"))

subfamily_state_sample_counts_pc$Sample.Proportion = ifelse(subfamily_state_sample_counts_pc$State %in% chromHMM_states,subfamily_state_sample_counts_pc$V1/sample_counts["All","chromHMM"],
                                                         ifelse(subfamily_state_sample_counts_pc$State %in% meth_states,subfamily_state_sample_counts_pc$V1/sample_counts["All","WGBS"],
                                                                ifelse(subfamily_state_sample_counts_pc$State == "DNase",subfamily_state_sample_counts_pc$V1/sample_counts["All","DNase"],
                                                                       ifelse(subfamily_state_sample_counts_pc$State == "H3K27ac",subfamily_state_sample_counts_pc$V1/sample_counts["All","H3K27ac"],"NA"))))
subfamily_state_sample_counts_pc$Sample.Proportion = as.numeric(subfamily_state_sample_counts_pc$Sample.Proportion)

rm(subfamily_state_sample)