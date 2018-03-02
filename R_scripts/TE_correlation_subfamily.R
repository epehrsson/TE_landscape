# Correlation of TE subfamily features with number of samples enriched in state

#source("R_scripts/TE_subfamily_stats.R")
#source("R_scripts/subfamily_enrichment.R")
#load("R_datasets/rna.RData")

# Adding number of samples enriched in each state
rmsk_TE_subfamily_measure = merge(rmsk_TE_subfamily,
                                  dcast(subfamily_state_sample_counts,class_update+family+subfamily~State,value.var=c("V1"))[,c("subfamily",chromHMM_states,meth_states,"DNase","H3K27ac")],
                                  by="subfamily",all.x=TRUE)

# Adding number of samples >1% in each state
test = dcast(subfamily_state_sample_counts_pc,class_update+family+subfamily~State,value.var=c("V1"))[,c("subfamily",chromHMM_states,meth_states,"DNase","H3K27ac")]
colnames(test)[2:22] = lapply(colnames(test)[2:22],function(x) paste(x,"PC",sep="_"))
rmsk_TE_subfamily_measure = merge(rmsk_TE_subfamily_measure,test,by="subfamily",all.x=TRUE)
rm(test)

# Adding C-GATE
rmsk_TE_subfamily_measure$CGate = rep("no",968)
rmsk_TE_subfamily_measure[which(rmsk_TE_subfamily_measure$subfamily %in% cgate_subfams),]$CGate = "yes"

contrasts(rmsk_TE_subfamily_measure$class_update) <- contr.sum
