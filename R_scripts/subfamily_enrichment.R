# Subfamily enrichment
# See 4/27/2016, 4/28/2016, 5/4/2016, 5/24/2016, 6/29/2016, 7/6/2016, 7/7/2016, 7/8/2016, 7/10/2016, 7/11/2016, 8/28/2016, 8/29/2016, 9/8/2016, 9/9/2016, 9/17/2016, 9/18/2016, 9/27/2016, 9/29/2016, 11/9/2016, 11/10/2016, 11/15/2016, 11/16/2016, 11/22/2016, 11/27/2016, 
# 1/13/2017, 1/19/2017, 1/20/2017, 1/26/2017, 2/6/2017, 2/9/2017, 2/10/2017, 2/13/2017, 2/14/2017, 2/16/2017, 2/21/2017, 3/16/2017, 5/15/2017, 5/16/2017, 5/23/2017, 5/24/2017, 5/29/2017, 5/30/2017, 6/5/17, 6/6/17, 6/7/17, 6/9/17, 6/15/2017, 7/4/2017, 7/21/2017, 7/23/2017, 7/24/2017, 7/26/2017, 8/2/2017

#source("R_scripts/TE_subfamily_stats.R")
#source("R_scripts/WGBS_sample_CpG_state.R")
#source("R_scripts/DNase_overlap.R")
#source("R_scripts/H3K27ac_overlap.R")
#load("R_datasets/TE_DNase_peaks.RData")
#load("R_datasets/TE_H3K27ac_peaks.RData")

# chromHMM
# Number of bp in each subfamily x sample x state
subfamily_state_sample = read.table("chromHMM/subfamily/subfamily_state_sample.txt",sep='\t')
colnames(subfamily_state_sample) = c("subfamily","State","Sample","Length_ijk")

# Adding subfamily x state x sample to matrix
subfamily_state_sample_expand = expand.grid(subfamily = levels(subfamily_state_sample$subfamily),Sample = levels(subfamily_state_sample$Sample),State = levels(subfamily_state_sample$State))
subfamily_state_sample_expand = join(subfamily_state_sample_expand,rmsk_TE_subfamily[,c("subfamily","family","class_update")],by=c("subfamily"))
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

# Add metadata for samples
subfamily_state_sample = join(subfamily_state_sample,metadata[,c(1,4:9)])

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

subfamily_state_sample = subfamily_state_sample[,c(1:2,5,3,4,12:17,6:11,19:20)]


# WGBS (CpGs)
# Proportion of subfamily CpGs in methylation state by sample
subfamily_CpG_meth = read.table("WGBS/subfamily_CpG_Meth_states.txt",sep='\t')
subfamily_CpG_meth$V1 = mapvalues(subfamily_CpG_meth$V1,seq(4,40,1),as.vector(metadata[which(!is.na(metadata$WGBS)),]$Sample))
subfamily_CpG_meth = subfamily_CpG_meth[,c(1,6,2:5)]
colnames(subfamily_CpG_meth) = c("Sample","subfamily",meth_states)
subfamily_CpG_meth[is.na(subfamily_CpG_meth)] = 0
subfamily_CpG_meth = melt(subfamily_CpG_meth,id.vars=c("Sample","subfamily"))
colnames(subfamily_CpG_meth)[3:4] = c("State","CpG_ijk")

# CpGs per subfamily
subfamily_CpG_meth = merge(subfamily_CpG_meth,rmsk_TE_subfamily[,c(1:3,33)],by=c("subfamily"),all.x=TRUE)
colnames(subfamily_CpG_meth)[7] = "CpG_ik"

# Proportion of all CpGs in methylation state by sample
subfamily_CpG_meth$CpG_jk = apply(subfamily_CpG_meth,1,function(x) all_CpG_meth[x[2],match(x[3],colnames(all_CpG_meth))])

# Enrichment of state CpGs in sample x subfamily
subfamily_CpG_meth$Enrichment = log2((subfamily_CpG_meth$CpG_ijk/subfamily_CpG_meth$CpG_ik)/(subfamily_CpG_meth$CpG_jk/ALL_CPGS))

# Proportion of all state CpGs in subfamily
subfamily_CpG_meth$CpG_ijk_jk = subfamily_CpG_meth$CpG_ijk/subfamily_CpG_meth$CpG_jk

# Adding metadata
subfamily_CpG_meth = merge(subfamily_CpG_meth,metadata[,c(1,4:9)],by=c("Sample"),all.x=TRUE)

# Adding members in state
subfamily_CpG_members = read.table("WGBS/subfamily_CpG_state_members.txt",sep='\t')
colnames(subfamily_CpG_members) = c("subfamily","Sample",meth_states)
subfamily_CpG_members = melt(subfamily_CpG_members,id.vars=c("subfamily","Sample"))
colnames(subfamily_CpG_members)[3:4] = c("State","Members")
subfamily_CpG_members[is.na(subfamily_CpG_members)] = 0

subfamily_CpG_meth = merge(subfamily_CpG_meth,subfamily_CpG_members,by=c("subfamily","Sample","State"),all.x=TRUE)
rm(subfamily_CpG_members)
subfamily_CpG_meth$Percent = apply(subfamily_CpG_meth,1,function(x) as.numeric(x[17])/rmsk_TE_subfamily[match(x[1],rmsk_TE_subfamily$subfamily),]$Count_CpGs)

subfamily_CpG_meth = subfamily_CpG_meth[,c(1,5:6,3,2,11:16,4,7:10,17:18)]


# DNase 
# Length of DNase peak overlapping subfamily per sample
subfamily_DNase_sample = read.table("DNase/subfamily_DNase_sample.txt",sep='\t')
colnames(subfamily_DNase_sample) = c("subfamily","Sample","Length_ijk")

# Length of subfamily
subfamily_DNase_sample_expand = expand.grid(subfamily = levels(rmsk_TE_subfamily$subfamily),Sample = levels(subfamily_DNase_sample$Sample))
subfamily_DNase_sample_expand = join(subfamily_DNase_sample_expand,rmsk_TE_subfamily[,c("subfamily","family","class_update","Total_length")],by=c("subfamily"))
subfamily_DNase_sample = join(subfamily_DNase_sample_expand,subfamily_DNase_sample,by=c("subfamily","Sample"),type="left")
colnames(subfamily_DNase_sample)[5] = "Length_ik"
subfamily_DNase_sample[which(metadata[match(subfamily_DNase_sample$Sample,metadata$Sample),]$chrY == "No"),]$Length_ik = apply(subfamily_DNase_sample[which(metadata[match(subfamily_DNase_sample$Sample,metadata$Sample),]$chrY == "No"),],1,function(x) rmsk_TE_subfamily[match(x[1],rmsk_TE_subfamily$subfamily),]$Total_length_noY)
subfamily_DNase_sample[is.na(subfamily_DNase_sample$Length_ijk),]$Length_ijk = 0
rm(subfamily_DNase_sample_expand)

# Total length of DNase peaks per sample
subfamily_DNase_sample = merge(subfamily_DNase_sample,DNase_stats[,c(1,4)],by=c("Sample"),all.x=TRUE)
colnames(subfamily_DNase_sample)[7] = "Length_jk"

# Total length of sample
mnemonics_states_genome = melt(as.matrix(t(read.table("chromHMM/genome/mnemonics_state.txt",sep='\t',header=TRUE,row.names=1)[1,])))[,c(1,3)]
colnames(mnemonics_states_genome) = c("Sample","Length_k")
subfamily_DNase_sample = join(subfamily_DNase_sample,mnemonics_states_genome,by=c("Sample"),type="left")
rm(mnemonics_states_genome)

# LOR enrichment
subfamily_DNase_sample$Enrichment = log2((subfamily_DNase_sample$Length_ijk/subfamily_DNase_sample$Length_ik)/(subfamily_DNase_sample$Length_jk/subfamily_DNase_sample$Length_k))

# Total proportion of DNase peaks in subfamily
subfamily_DNase_sample$Length_percent_jk = subfamily_DNase_sample$Length_ijk/subfamily_DNase_sample$Length_jk

# Add metadata
subfamily_DNase_sample = merge(subfamily_DNase_sample,metadata[,c(1,4:9)],by=c("Sample"),all.x=TRUE)

# Number of members overlapping DNase peak
TE_DNase_peaks_members = merge(melt(aggregate(data=TE_DNase_peaks[,c(4,8:60)],.~subfamily,function(x) sum(as.numeric(x))),id.vars=c("subfamily")),melt(aggregate(data=TE_DNase_peaks[,c(4,8:60)],.~subfamily,function(x) sum(as.numeric(x) > 0)),id.vars=c("subfamily")),by=c("subfamily","variable"))
colnames(TE_DNase_peaks_members)[2:4] = c("Sample","DNase_peaks","Members")
TE_DNase_peaks_members = merge(TE_DNase_peaks_members,rmsk_TE_subfamily[,c("subfamily","Count")],by="subfamily")
TE_DNase_peaks_members[which(TE_DNase_peaks_members$Sample %in% c("E116","E117","E123","E124","E126","E127")),]$Count = rmsk_TE_subfamily[match(TE_DNase_peaks_members[which(TE_DNase_peaks_members$Sample %in% c("E116","E117","E123","E124","E126","E127")),]$subfamily,rmsk_TE_subfamily$subfamily),]$Count_noY

# Proportion of members overlapping DNase peak
TE_DNase_peaks_members$Percent = TE_DNase_peaks_members$Members/TE_DNase_peaks_members$Count

# Added to enrichment matrix
subfamily_DNase_sample = merge(subfamily_DNase_sample,TE_DNase_peaks_members,by=c("subfamily","Sample"),all.x=TRUE)
subfamily_DNase_sample[which(is.na(subfamily_DNase_sample$Percent)),17:20] = 0
subfamily_DNase_sample[which(subfamily_DNase_sample$subfamily == "Tigger2a_Car"),]$Count = 2
subfamily_DNase_sample$State = rep("DNase",dim(subfamily_DNase_sample)[1])

subfamily_DNase_sample = subfamily_DNase_sample[,c(1,3:4,21,2,11:16,5:10,18,20,17,19)]


# H3K27ac 
# Length of H3K27ac peak overlapping subfamily per sample
subfamily_H3K27ac_sample = read.table("H3K27ac/subfamily_H3K27ac_sample.txt",sep='\t')
colnames(subfamily_H3K27ac_sample) = c("subfamily","Sample","Length_ijk")

# Length of subfamily
subfamily_H3K27ac_sample_expand = expand.grid(subfamily = levels(rmsk_TE_subfamily$subfamily),Sample = levels(subfamily_H3K27ac_sample$Sample))
subfamily_H3K27ac_sample_expand = join(subfamily_H3K27ac_sample_expand,rmsk_TE_subfamily[,c("subfamily","family","class_update","Total_length")],by=c("subfamily"))
subfamily_H3K27ac_sample = join(subfamily_H3K27ac_sample_expand,subfamily_H3K27ac_sample,by=c("subfamily","Sample"),type="left")
colnames(subfamily_H3K27ac_sample)[5] = "Length_ik"
subfamily_H3K27ac_sample[which(metadata[match(subfamily_H3K27ac_sample$Sample,metadata$Sample),]$chrY == "No"),]$Length_ik = apply(subfamily_H3K27ac_sample[which(metadata[match(subfamily_H3K27ac_sample$Sample,metadata$Sample),]$chrY == "No"),],1,function(x) rmsk_TE_subfamily[match(x[1],rmsk_TE_subfamily$subfamily),]$Total_length_noY)
subfamily_H3K27ac_sample[is.na(subfamily_H3K27ac_sample$Length_ijk),]$Length_ijk = 0
rm(subfamily_H3K27ac_sample_expand)

# Total length of H3K27ac peaks per sample
subfamily_H3K27ac_sample = merge(subfamily_H3K27ac_sample,H3K27ac_stats[,c(1,4)],by=c("Sample"),all.x=TRUE)
colnames(subfamily_H3K27ac_sample)[7] = "Length_jk"

# Total length of sample
mnemonics_states_genome = melt(as.matrix(t(read.table("chromHMM/genome/mnemonics_state.txt",sep='\t',header=TRUE,row.names=1)[1,])))[,c(1,3)]
colnames(mnemonics_states_genome) = c("Sample","Length_k")
subfamily_H3K27ac_sample = join(subfamily_H3K27ac_sample,mnemonics_states_genome,by=c("Sample"),type="left")
rm(mnemonics_states_genome)

# LOR enrichment
subfamily_H3K27ac_sample$Enrichment = log2((subfamily_H3K27ac_sample$Length_ijk/subfamily_H3K27ac_sample$Length_ik)/(subfamily_H3K27ac_sample$Length_jk/subfamily_H3K27ac_sample$Length_k))

# Total proportion of H3K27ac peaks in subfamily
subfamily_H3K27ac_sample$Length_percent_jk = subfamily_H3K27ac_sample$Length_ijk/subfamily_H3K27ac_sample$Length_jk

# Add metadata
subfamily_H3K27ac_sample = merge(subfamily_H3K27ac_sample,metadata[,c(1,4:9)],by=c("Sample"),all.x=TRUE)

# Number of members overlapping H3K27ac peak
TE_H3K27ac_peaks_members = merge(melt(aggregate(data=TE_H3K27ac_peaks[,c(4,8:105)],.~subfamily,function(x) sum(as.numeric(x))),id.vars=c("subfamily")),melt(aggregate(data=TE_H3K27ac_peaks[,c(4,8:105)],.~subfamily,function(x) sum(as.numeric(x) > 0)),id.vars=c("subfamily")),by=c("subfamily","variable"))
colnames(TE_H3K27ac_peaks_members)[2:4] = c("Sample","H3K27ac_peaks","Members")
TE_H3K27ac_peaks_members = merge(TE_H3K27ac_peaks_members,rmsk_TE_subfamily[,c("subfamily","Count")],by="subfamily")
TE_H3K27ac_peaks_members[which(TE_H3K27ac_peaks_members$Sample %in% c("E116","E117","E123","E124","E126","E127")),]$Count = rmsk_TE_subfamily[match(TE_H3K27ac_peaks_members[which(TE_H3K27ac_peaks_members$Sample %in% c("E116","E117","E123","E124","E126","E127")),]$subfamily,rmsk_TE_subfamily$subfamily),]$Count_noY

# Proportion of members overlapping H3K27ac peak
TE_H3K27ac_peaks_members$Percent = TE_H3K27ac_peaks_members$Members/TE_H3K27ac_peaks_members$Count

# Added to enrichment matrix
subfamily_H3K27ac_sample = merge(subfamily_H3K27ac_sample,TE_H3K27ac_peaks_members,by=c("subfamily","Sample"),all.x=TRUE)
subfamily_H3K27ac_sample[which(is.na(subfamily_H3K27ac_sample$Percent)),17:20] = 0
subfamily_H3K27ac_sample[which(subfamily_H3K27ac_sample$subfamily == "Tigger2a_Car"),]$Count = 2
subfamily_H3K27ac_sample$State = rep("H3K27ac",dim(subfamily_H3K27ac_sample)[1])

subfamily_H3K27ac_sample = subfamily_H3K27ac_sample[,c(1,3:4,21,2,11:16,5:10,18,20,17,19)]

# Combine matrices (no filtering)
test = subfamily_CpG_meth[,c(1:13,15:18)]
colnames(test)[c(12:13,15)] = c("Length_ijk","Length_ik","Length_percent_jk")
subfamily_state_sample_combined = rbind(subfamily_state_sample[,c(1:13,16:19)],test,subfamily_DNase_sample[,c(1:13,16:19)],subfamily_H3K27ac_sample[,c(1:13,16:19)])

# Combine filtered matrices
test2 = subfamily_CpG_meth[which(subfamily_CpG_meth$CpG_ijk >= THRESHOLD_IJK_CPG & subfamily_CpG_meth$CpG_ik > THRESHOLD_IK_CPG & subfamily_CpG_meth$Members >= THRESHOLD_IJK_MEMBERS),]
colnames(test2)[c(12:13,16)] = c("Length_ijk","Length_ik","Length_percent_jk")
subfamily_state_sample_filter = rbind(subfamily_state_sample[which(subfamily_state_sample$Length_ijk >= THRESHOLD_IJK_BASE & subfamily_state_sample$Length_ik > THRESHOLD_IK_BASE & subfamily_state_sample$Members >= THRESHOLD_IJK_MEMBERS),c(1:13,16:19)],
                                      test2[,c(1:13,15:18)],
                                      subfamily_DNase_sample[which(subfamily_DNase_sample$Length_ijk >= THRESHOLD_IJK_BASE & subfamily_DNase_sample$Length_ik > THRESHOLD_IK_BASE & subfamily_DNase_sample$Members >= THRESHOLD_IJK_MEMBERS),c(1:13,16:19)],
                                      subfamily_H3K27ac_sample[which(subfamily_H3K27ac_sample$Length_ijk >= THRESHOLD_IJK_BASE & subfamily_H3K27ac_sample$Length_ik > THRESHOLD_IK_BASE & subfamily_H3K27ac_sample$Members >= THRESHOLD_IJK_MEMBERS),c(1:13,16:19)])
subfamily_state_sample_filter$State = factor(subfamily_state_sample_filter$State,levels=c(chromHMM_states,meth_states,"DNase","H3K27ac"))

# Number of enrichments per subfamily x state
subfamily_state_sample_counts = ddply(subfamily_state_sample_filter,.(class_update,family,subfamily,State),function(x) sum(x$Enrichment > THRESHOLD_LOR))
subfamily_state_expand = expand.grid(subfamily = levels(subfamily_state_sample$subfamily),State = levels(subfamily_state_sample_combined$State))
subfamily_state_expand = join(subfamily_state_expand,rmsk_TE_subfamily[,c("subfamily","family","class_update")],by=c("subfamily"))
subfamily_state_sample_counts = join(subfamily_state_expand,subfamily_state_sample_counts,by=c("subfamily","family","class_update","State"),type="left")[,c(1,3,4,2,5)]
rm(subfamily_state_expand)
subfamily_state_sample_counts[is.na(subfamily_state_sample_counts)] = 0
subfamily_state_sample_counts$State = factor(subfamily_state_sample_counts$State,levels=c(chromHMM_states,meth_states,"DNase","H3K27ac"))

subfamily_state_sample_counts_combine = merge(aggregate(data=subfamily_state_sample_counts,V1~State,function(x) sum(x > 0)),
                                              dcast(aggregate(data=subfamily_state_sample_counts,V1~State+class_update,function(x) sum(x > 0)),State~class_update,value.var = "V1"),
                                              by=c("State"))
subfamily_state_sample_counts_combine = subfamily_state_sample_counts_combine[match(levels(subfamily_state_sample_counts$State),subfamily_state_sample_counts_combine$State),]

# Number of >1% per subfamily x state
subfamily_state_sample_counts_pc = ddply(subfamily_state_sample_combined,.(class_update,family,subfamily,State),function(x) sum(x$Length_percent_jk > THRESHOLD_PC))
subfamily_state_sample_counts_pc_combine = merge(aggregate(data=subfamily_state_sample_counts_pc,V1~State,function(x) sum(x > 0)),
                                              dcast(aggregate(data=subfamily_state_sample_counts_pc,V1~State+class_update,function(x) sum(x > 0)),State~class_update,value.var = "V1"),
                                              by=c("State"))
