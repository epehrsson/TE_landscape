# Correlation of TE subfamily features with number of samples enriched in state

library(reshape2)

source("R_scripts/TE_subfamily_stats.R")

# Number of samples enriched by subfamily
source("R_scripts/chromHMM_subfamily_enrichment.R")
source("R_scripts/WGBS_subfamily_enrichment_CpG.R")
source("R_scripts/DNase_subfamily_enrichment.R")
source("R_scripts/H3K27ac_subfamily_enrichment.R")

# Adding number of samples enriched in chromHMM states
rmsk_TE_subfamily_measure = merge(rmsk_TE_subfamily,dcast(subfamily_state_sample_counts,class_update+family+subfamily~State,value.var=c("V1"))[,c(3,10:17,4:9)],by="subfamily",all.x=TRUE)

# How to add chromHMM stats?
ddply(subfamily_state_sample,~subfamily+State,summarise,State_Max = max(Percent),State_Mean = mean(Percent),State_Range = (max(Percent)-min(Percent)))

# Adding hypomethylated CpG enrichment
rmsk_TE_subfamily_measure = merge(rmsk_TE_subfamily_measure,dcast(subfamily_hypo_sample_counts,class_update+family+subfamily~State,value.var=c("V1"))[,3:4],by="subfamily",all.x=TRUE)

# Hypomethylated proportion mean, max, range per subfamily
source('~/TE_landscape/R_scripts/WGBS_subfamily_enrichment_TE.R')

rmsk_TE_subfamily_measure = merge(rmsk_TE_subfamily_measure,TE_meth_subfamily_hypo[,c(1:3,42:47)],by=c("subfamily","family","class_update"),all.x=TRUE)

# Adding Dnase enrichment
rmsk_TE_subfamily_measure = merge(rmsk_TE_subfamily_measure,dcast(subfamily_DNase_sample_counts,class_update+family+subfamily~State,value.var=c("V1"))[,3:4],by="subfamily",all.x=TRUE)

# Add DNase proportion mean, max, range
rmsk_TE_subfamily_measure = merge(rmsk_TE_subfamily_measure,ddply(subfamily_DNase_sample,~subfamily,summarise,DNase_Max = max(Percent),DNase_Mean = mean(Percent),DNase_Range = (max(Percent)-min(Percent))),by="subfamily")

# Adding H3K27ac enrichment
rmsk_TE_subfamily_measure = merge(rmsk_TE_subfamily_measure,dcast(subfamily_H3K27ac_sample_counts,class_update+family+subfamily~State,value.var=c("V1"))[,3:4],by="subfamily",all.x=TRUE)

# Add H3K27ac proportion mean, max, range
rmsk_TE_subfamily_measure = merge(rmsk_TE_subfamily_measure,ddply(subfamily_H3K27ac_sample,~subfamily,summarise,H3K27ac_Max = max(Percent),H3K27ac_Mean = mean(Percent),H3K27ac_Range = (max(Percent)-min(Percent))),by="subfamily")