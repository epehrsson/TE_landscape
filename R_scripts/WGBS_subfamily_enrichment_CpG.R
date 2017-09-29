# Enrichment of hypomethylated CpGs in a subfamily
# See 2/6/2017, 7/21/2017, 7/23/2017, 7/26/2017, 8/2/2017

source("R_scripts/TE_subfamily_stats.R")
source("R_scripts/WGBS_sample_CpG_state.R")

# Proportion of subfamily CpGs in methylation state by sample
subfamily_CpG_meth = read.table("WGBS/subfamily_CpG_Meth_states.txt",sep='\t')
subfamily_CpG_meth$V1 = mapvalues(subfamily_CpG_meth$V1,seq(4,40,1),as.vector(metadata[which(!is.na(metadata$WGBS)),]$Sample))
subfamily_CpG_meth = subfamily_CpG_meth[,c(1,6,2:5)]
colnames(subfamily_CpG_meth) = c("Sample","subfamily",meth_states)
subfamily_CpG_meth[is.na(subfamily_CpG_meth)] = 0
subfamily_CpG_meth = melt(subfamily_CpG_meth,id.vars=c("Sample","subfamily"))
colnames(subfamily_CpG_meth)[3:4] = c("State","CpG_ijk")

# CpGs per subfamily
subfamily_CpG_meth = merge(subfamily_CpG_meth,rmsk_TE_subfamily[,c(1:3,32)],by=c("subfamily"),all.x=TRUE)
colnames(subfamily_CpG_meth)[7] = "CpG_ik"

# Proportion of all CpGs in methylation state by sample
subfamily_CpG_meth$CpG_jk = apply(subfamily_CpG_meth,1,function(x) all_CpG_meth[x[2],match(x[3],colnames(all_CpG_meth))])

# Enrichment of state CpGs in sample x subfamily
subfamily_CpG_meth$Enrichment = log2((subfamily_CpG_meth$CpG_ijk/subfamily_CpG_meth$CpG_ik)/(subfamily_CpG_meth$CpG_jk/56434896))

# Proportion of all state CpGs in subfamily
subfamily_CpG_meth$CpG_ijk_jk = subfamily_CpG_meth$CpG_ijk/subfamily_CpG_meth$CpG_jk

# Adding metadata
subfamily_CpG_meth = merge(subfamily_CpG_meth,metadata[,c(1,4:9)],by=c("Sample"),all.x=TRUE)
subfamily_CpG_meth = subfamily_CpG_meth[,c(2,5:6,1,3,4,7:16)]

# Adding members in state
subfamily_CpG_members = read.table("WGBS/subfamily_CpG_state_members.txt",sep='\t')
colnames(subfamily_CpG_members) = c("subfamily","Sample",meth_states)
subfamily_CpG_members = melt(subfamily_CpG_members,id.vars=c("subfamily","Sample"))
colnames(subfamily_CpG_members)[3:4] = c("State","Members")
subfamily_CpG_members[is.na(subfamily_CpG_members)] = 0

subfamily_CpG_meth = merge(subfamily_CpG_meth,subfamily_CpG_members,by=c("subfamily","Sample","State"),all.x=TRUE)
subfamily_CpG_meth$Percent = apply(subfamily_CpG_meth,1,function(x) as.numeric(x[17])/rmsk_TE_subfamily[match(x[1],rmsk_TE_subfamily$subfamily),]$Count_CpGs)

# Number of hypomethylation enrichments per subfamily
subfamily_hypo_sample_counts = ddply(subfamily_CpG_meth,.(class_update,family,subfamily,State),function(x) sum(x$Enrichment > 1.5 & x$CpG_ik >= 25 & x$CpG_ijk >= 6))
subfamily_hypo_sample_counts_noIMR90 = ddply(subfamily_CpG_meth,.(class_update,family,subfamily,State),function(x) sum(x$Enrichment > 1.5 & x$CpG_ik >= 25 & x$CpG_ijk >= 6 & x$Sample != "E017"))    

# Number of >1% per subfamily x state
subfamily_hypo_sample_counts_pc = ddply(subfamily_CpG_meth,.(class_update,family,subfamily,State),function(x) sum(x$CpG_ijk_jk > 0.01 & x$CpG_ijk >= 6))
subfamily_hypo_sample_counts_pc_IMR90 = ddply(subfamily_CpG_meth,.(class_update,family,subfamily,State),function(x) sum(x$CpG_ijk_jk > 0.01 & x$CpG_ijk >= 6))