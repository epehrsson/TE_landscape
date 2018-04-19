# Identify sample group-specific TEs within each subfamily.
# For every category grouping, counts the number of subfamily TEs that are in the state in at least xx% of that grouping, excluding those that are not very specific to the grouping. 
# Allows thresholding of the state coverage of individual TEs and the number of samples required to consider a sample grouping (>1).

subfamilies_1TssA = as.vector(unique(subfamily_state_sample[which(subfamily_state_sample$State == "1_TssA" 
                                                                  & subfamily_state_sample$Members > 0),]$subfamily))
clusters_1_TssA = lapply(subfamilies_1TssA,function(x) enrichment_clusters(x,'1_TssA','chromHMM',
                                                                           coverage_threshold=0,score_threshold=3,category_threshold=0.33,sample_threshold=1,TE_threshold=3))
names(clusters_1_TssA) = subfamilies_1TssA
clusters_1_TssA = ldply(clusters_1_TssA)
colnames(clusters_1_TssA)[1] = "subfamily"
clusters_1_TssA = merge(clusters_1_TssA,rmsk_TE_subfamily[,c(1,4)],by="subfamily",all.x=TRUE)
clusters_1_TssA$TE_percent = clusters_1_TssA$TEs/clusters_1_TssA$Count
write.table(clusters_1_TssA,file="enrichment/clusters_1_TssA.txt",sep='\t',row.names = FALSE,quote=FALSE)

# Printed out clusters for 1_TssA (all) and 4_Tx, 5_TxWk, 6_Enh, H3K27ac, 7_Enh (subfamilies with <5000 members)

# Min, max, and mean proportion of all TEs in the state for each sample grouping for all states of interest
test = merge(combine_boxplot,metadata[,c(1,4:9)],by="Sample",all.x=TRUE)
test = melt(test,id.vars=c("Sample","State","Proportion"))
colnames(test)[4:5] = c("Category","Grouping")
ddply(test[which(test$State %in% c(chromHMM_states[1:7],meth_states[1:2],"DNase","H3K27ac")),],.(State,Category,Grouping),summarise,Min=min(Proportion),Max=max(Proportion),Mean=mean(Proportion))
