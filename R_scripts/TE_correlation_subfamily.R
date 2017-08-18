# Hypomethylated proportion mean, max per subfamily
rmsk_TEother_age_subfamily = merge(rmsk_TEother_age_subfamily,TE_meth_subfamily_hypo[,c(1:3,41:42,45:46)],by=c("subfamily","family","class"),all.x=TRUE)
colnames(rmsk_TEother_age_subfamily)[29:30] = c("Mean_hypo","Mean_hypo_noIMR90")
colnames(rmsk_TEother_age_subfamily)[31:32] = c("Max_hypo","Max_hypo_noIMR90")

# Number of samples enriched by subfamily
# Adding number of samples enriched in chromHMM states
rmsk_TEother_age_subfamily = merge(rmsk_TEother_age_subfamily,dcast(subfamily_state_sample_counts,Class+Family+Subfamily~State,value.var=c("V1"))[,3:18],by.x="subfamily",by.y="Subfamily",all.x=TRUE)

# Correlation between hypomethylated CpG enrichment, age
rmsk_TEother_age_subfamily = merge(rmsk_TEother_age_subfamily,subfamily_hypo_sample_counts[,3:4],by=c("subfamily"))
colnames(rmsk_TEother_age_subfamily)[33] = "CpG_Hypo"

# Adding Dnase enrichment
rmsk_TEother_age_subfamily = merge(rmsk_TEother_age_subfamily,subfamily_DNase_sample_counts[,3:4],by.x="subfamily",by.y="Subfamily",all.x=TRUE)
colnames(rmsk_TEother_age_subfamily)[32] = "DNase"

# Adding H3K27ac
rmsk_TEother_age_subfamily = merge(rmsk_TEother_age_subfamily,subfamily_H3K27ac_sample_counts[,3:4],by.x="subfamily",by.y="Subfamily",all.x=TRUE)
colnames(rmsk_TEother_age_subfamily)[28] = "H3K27ac"