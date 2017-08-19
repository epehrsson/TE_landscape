# Correlation of TE subfamily features with number of samples enriched in state

source("R_scripts/TE_subfamily_stats.R")

# Hypomethylated proportion mean, max, range per subfamily
source('~/TE_landscape/R_scripts/WGBS_subfamily_enrichment_TE.R')

rmsk_TE_subfamily = merge(rmsk_TE_subfamily,TE_meth_subfamily_hypo[,c(1:3,42:47)],by=c("subfamily","family","class_update"),all.x=TRUE)

# Number of samples enriched by subfamily
# Adding number of samples enriched in chromHMM states
rmsk_TE_subfamily = merge(rmsk_TE_subfamily,dcast(subfamily_state_sample_counts,Class+Family+Subfamily~State,value.var=c("V1"))[,3:18],by.x="subfamily",by.y="Subfamily",all.x=TRUE)

# Correlation between hypomethylated CpG enrichment, age
rmsk_TE_subfamily = merge(rmsk_TE_subfamily,subfamily_hypo_sample_counts[,3:4],by=c("subfamily"))
colnames(rmsk_TE_subfamily)[33] = "CpG_Hypo"

# Adding Dnase enrichment
rmsk_TE_subfamily = merge(rmsk_TE_subfamily,subfamily_DNase_sample_counts[,3:4],by.x="subfamily",by.y="Subfamily",all.x=TRUE)
colnames(rmsk_TE_subfamily)[32] = "DNase"

# Adding H3K27ac
rmsk_TE_subfamily = merge(rmsk_TE_subfamily,subfamily_H3K27ac_sample_counts[,3:4],by.x="subfamily",by.y="Subfamily",all.x=TRUE)
colnames(rmsk_TE_subfamily)[28] = "H3K27ac"