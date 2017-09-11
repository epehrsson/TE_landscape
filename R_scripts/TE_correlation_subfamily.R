# Correlation of TE subfamily features with number of samples enriched in state

library(reshape2)

source("R_scripts/TE_subfamily_stats.R")
load("R_datasets/rna.RData")

# Number of samples enriched by subfamily
source("R_scripts/chromHMM_subfamily_enrichment.R")
source("R_scripts/WGBS_subfamily_enrichment_CpG.R")
source('~/TE_landscape/R_scripts/WGBS_subfamily_enrichment_TE.R')
source("R_scripts/DNase_subfamily_enrichment.R")
source("R_scripts/H3K27ac_subfamily_enrichment.R")

# Adding number of samples enriched in chromHMM states
rmsk_TE_subfamily_measure = merge(rmsk_TE_subfamily,dcast(subfamily_state_sample_counts,class_update+family+subfamily~State,value.var=c("V1"))[,c(3:4,11:18,5:10)],by="subfamily",all.x=TRUE)

# Adding number of samples >1% in chromHMM states
test = dcast(subfamily_state_sample_counts_pc,class_update+family+subfamily~State,value.var=c("V1"))[,c(3:4,11:18,5:10)]
colnames(test)[2:16] = lapply(colnames(test)[2:16],function(x) paste(x,"PC",sep="_"))
rmsk_TE_subfamily_measure = merge(rmsk_TE_subfamily_measure,test,by="subfamily",all.x=TRUE)
rm(test)

# chromHMM proportion mean, max, range per subfamily
test = melt(ddply(subfamily_state_sample,~subfamily+State,summarise,Max = max(Percent),Mean = mean(Percent),Range = (max(Percent)-min(Percent))),id.vars=c("subfamily","State"))
test$Header = apply(test,1,function(x) paste(x[2],x[3],sep="_"))
test = dcast(test[,c(1,4:5)],subfamily~Header)
rmsk_TE_subfamily_measure = merge(rmsk_TE_subfamily_measure,test,by="subfamily")
rm(test)

# Adding methylation state CpG enrichment
rmsk_TE_subfamily_measure = merge(rmsk_TE_subfamily_measure,dcast(subfamily_hypo_sample_counts,class_update+family+subfamily~State,value.var=c("V1"))[,3:11],by="subfamily",all.x=TRUE)

# Adding number of samples >1% in methylation state
test = dcast(subfamily_hypo_sample_counts_pc,class_update+family+subfamily~State,value.var=c("V1"))[,3:11]
colnames(test)[2:9] = lapply(colnames(test)[2:9],function(x) paste(x,"PC",sep="_"))
rmsk_TE_subfamily_measure = merge(rmsk_TE_subfamily_measure,test,by="subfamily",all.x=TRUE)
rm(test)

# Methylation state proportion mean, max, range per subfamily (should be CpGs!)
for (i in 1:4){
  rmsk_TE_subfamily_measure = merge(rmsk_TE_subfamily_measure,TE_meth_subfamily_hypo[[i]][,c(1:3,41:46)],by=c("subfamily","family","class_update"),all.x=TRUE)
  z = length(colnames(rmsk_TE_subfamily_measure))
  colnames(rmsk_TE_subfamily_measure)[(z-5):z] = lapply(colnames(rmsk_TE_subfamily_measure)[(z-5):z],function(x) paste(names(TE_meth_subfamily_hypo)[i],x,sep="_"))
}

# Adding Dnase enrichment
rmsk_TE_subfamily_measure = merge(rmsk_TE_subfamily_measure,dcast(subfamily_DNase_sample_counts,class_update+family+subfamily~State,value.var=c("V1"))[,3:4],by="subfamily",all.x=TRUE)

# Adding DNase >1%
rmsk_TE_subfamily_measure = merge(rmsk_TE_subfamily_measure,dcast(subfamily_DNase_sample_counts_pc,class_update+family+subfamily~State,value.var=c("V1"))[,3:4],by="subfamily",all.x=TRUE)
z = length(colnames(rmsk_TE_subfamily_measure))
colnames(rmsk_TE_subfamily_measure)[(z-1):z] = c("DNase","DNase_PC")

# Add DNase proportion mean, max, range
rmsk_TE_subfamily_measure = merge(rmsk_TE_subfamily_measure,ddply(subfamily_DNase_sample,~subfamily,summarise,DNase_Max = max(Percent),DNase_Mean = mean(Percent),DNase_Range = (max(Percent)-min(Percent))),by="subfamily")

# Adding H3K27ac enrichment
rmsk_TE_subfamily_measure = merge(rmsk_TE_subfamily_measure,dcast(subfamily_H3K27ac_sample_counts,class_update+family+subfamily~State,value.var=c("V1"))[,3:4],by="subfamily",all.x=TRUE)

# Adding H3K27ac >1%
rmsk_TE_subfamily_measure = merge(rmsk_TE_subfamily_measure,dcast(subfamily_H3K27ac_sample_counts_pc,class_update+family+subfamily~State,value.var=c("V1"))[,3:4],by="subfamily",all.x=TRUE)
z = length(colnames(rmsk_TE_subfamily_measure))
colnames(rmsk_TE_subfamily_measure)[(z-1):z] = c("H3K27ac","H3K27ac_PC")

# Add H3K27ac proportion mean, max, range
rmsk_TE_subfamily_measure = merge(rmsk_TE_subfamily_measure,ddply(subfamily_H3K27ac_sample,~subfamily,summarise,H3K27ac_Max = max(Percent),H3K27ac_Mean = mean(Percent),H3K27ac_Range = (max(Percent)-min(Percent))),by="subfamily")

# Adding RNA expression 
rmsk_TE_subfamily_measure = merge(rmsk_TE_subfamily_measure,RNA_TE_agnostic_subfamily[,c(3,57:58)],by="subfamily")
