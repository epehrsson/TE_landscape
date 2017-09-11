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

# CpGs per subfamily
subfamily_CpG_meth = merge(subfamily_CpG_meth,rmsk_TE_subfamily[,c(1:3,32)],by=c("subfamily"),all.x=TRUE)

# Proportion of all CpGs in methylation state by sample
subfamily_CpG_meth$Hypomethylated_all = apply(subfamily_CpG_meth,1,function(x) all_CpG_meth[x[2],]$Hypomethylated)
subfamily_CpG_meth$Intermediate_all = apply(subfamily_CpG_meth,1,function(x) all_CpG_meth[x[2],]$Intermediate)
subfamily_CpG_meth$Hypermethylated_all = apply(subfamily_CpG_meth,1,function(x) all_CpG_meth[x[2],]$Hypermethylated)
subfamily_CpG_meth$Missing_all = apply(subfamily_CpG_meth,1,function(x) all_CpG_meth[x[2],]$Missing)

# Enrichment of hypomethylated CpGs in sample x subfamily
subfamily_CpG_meth$Hypo_enrichment = log2((subfamily_CpG_meth$Hypomethylated/subfamily_CpG_meth$CpGs)/(subfamily_CpG_meth$Hypomethylated_all/56434896))
subfamily_CpG_meth$Inter_enrichment = log2((subfamily_CpG_meth$Intermediate/subfamily_CpG_meth$CpGs)/(subfamily_CpG_meth$Intermediate_all/56434896))
subfamily_CpG_meth$Hyper_enrichment = log2((subfamily_CpG_meth$Hypermethylated/subfamily_CpG_meth$CpGs)/(subfamily_CpG_meth$Hypermethylated_all/56434896))
subfamily_CpG_meth$Missing_enrichment = log2((subfamily_CpG_meth$Missing/subfamily_CpG_meth$CpGs)/(subfamily_CpG_meth$Missing_all/56434896))

# Proportion of all hypomethylated CpGs in subfamily
subfamily_CpG_meth$Hypo_percent = subfamily_CpG_meth$Hypomethylated/subfamily_CpG_meth$Hypomethylated_all
subfamily_CpG_meth$Inter_percent = subfamily_CpG_meth$Intermediate/subfamily_CpG_meth$Intermediate_all
subfamily_CpG_meth$Hyper_percent = subfamily_CpG_meth$Hypermethylated/subfamily_CpG_meth$Hypermethylated_all
subfamily_CpG_meth$Missing_percent = subfamily_CpG_meth$Missing/subfamily_CpG_meth$Missing_all

# Adding metadata
subfamily_CpG_meth = merge(subfamily_CpG_meth,metadata[,c(1,4:9)],by=c("Sample"),all.x=TRUE)
subfamily_CpG_meth = subfamily_CpG_meth[,c(2,7:8,1,3:6,9:27)]

# Number of hypomethylation enrichments per subfamily
subfamily_hypo_sample_counts = rbind(ddply(subfamily_CpG_meth,.(class_update,family,subfamily),function(x) sum(x$Hypo_enrichment > 1.5 & x$CpGs >= 25 & x$Hypomethylated >= 6)),
                                     ddply(subfamily_CpG_meth,.(class_update,family,subfamily),function(x) sum(x$Inter_enrichment > 1.5 & x$CpGs >= 25 & x$Intermediate >= 6)),
                                     ddply(subfamily_CpG_meth,.(class_update,family,subfamily),function(x) sum(x$Hyper_enrichment > 1.5 & x$CpGs >= 25 & x$Hypermethylated >= 6)),
                                     ddply(subfamily_CpG_meth,.(class_update,family,subfamily),function(x) sum(x$Missing_enrichment > 1.5 & x$CpGs >= 25 & x$Missing >= 6)),
                                     ddply(subfamily_CpG_meth,.(class_update,family,subfamily),function(x) sum(x$Hypo_enrichment > 1.5 & x$CpGs >= 25 & x$Hypomethylated >= 6 & x$Sample != "E017")),
                                     ddply(subfamily_CpG_meth,.(class_update,family,subfamily),function(x) sum(x$Inter_enrichment > 1.5 & x$CpGs >= 25 & x$Intermediate >= 6 & x$Sample != "E017")),
                                     ddply(subfamily_CpG_meth,.(class_update,family,subfamily),function(x) sum(x$Hyper_enrichment > 1.5 & x$CpGs >= 25 & x$Hypermethylated >= 6 & x$Sample != "E017")),
                                     ddply(subfamily_CpG_meth,.(class_update,family,subfamily),function(x) sum(x$Missing_enrichment > 1.5 & x$CpGs >= 25 & x$Missing >= 6 & x$Sample != "E017")))
subfamily_hypo_sample_counts$State = c(rep("Hypomethylated",968),rep("Intermediate",968),rep("Hypermethylated",968),rep("Missing",968),rep("Hypomethylated_noIMR90",968),rep("Intermediate_noIMR90",968),rep("Hypermethylated_noIMR90",968),rep("Missing_noIMR90",968))

# Number of >1% per subfamily x state
subfamily_hypo_sample_counts_pc = rbind(ddply(subfamily_CpG_meth,.(class_update,family,subfamily),function(x) sum(x$Hypo_percent > 0.01 & x$Hypomethylated >= 6)),
                                     ddply(subfamily_CpG_meth,.(class_update,family,subfamily),function(x) sum(x$Inter_percent > 0.01 & x$Intermediate >= 6)),
                                     ddply(subfamily_CpG_meth,.(class_update,family,subfamily),function(x) sum(x$Hyper_percent > 0.01 & x$Hypermethylated >= 6)),
                                     ddply(subfamily_CpG_meth,.(class_update,family,subfamily),function(x) sum(x$Missing_percent > 0.01 & x$Missing >= 6)),
                                     ddply(subfamily_CpG_meth,.(class_update,family,subfamily),function(x) sum(x$Hypo_percent > 0.01 & x$Hypomethylated >= 6 & x$Sample != "E017")),
                                     ddply(subfamily_CpG_meth,.(class_update,family,subfamily),function(x) sum(x$Inter_percent > 0.01 & x$Intermediate >= 6 & x$Sample != "E017")),
                                     ddply(subfamily_CpG_meth,.(class_update,family,subfamily),function(x) sum(x$Hyper_percent > 0.01 & x$Hypermethylated >= 6 & x$Sample != "E017")),
                                     ddply(subfamily_CpG_meth,.(class_update,family,subfamily),function(x) sum(x$Missing_percent > 0.01 & x$Missing >= 6 & x$Sample != "E017")))
subfamily_hypo_sample_counts_pc$State = c(rep("Hypomethylated",968),rep("Intermediate",968),rep("Hypermethylated",968),rep("Missing",968),rep("Hypomethylated_noIMR90",968),rep("Intermediate_noIMR90",968),rep("Hypermethylated_noIMR90",968),rep("Missing_noIMR90",968))
