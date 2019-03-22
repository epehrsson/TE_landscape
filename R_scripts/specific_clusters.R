# Identify sample group-specific TEs within each subfamily.
# For every category grouping, counts the number of subfamily TEs that are in the state in at least xx% of that grouping, excluding those that are not very specific to the grouping. 
# Allows thresholding of the state coverage of individual TEs and the number of samples required to consider a sample grouping (>1).

enrichment_clusters = function(subfamily,state,metric,coverage_threshold=0,score_threshold=1,category_threshold=0.5,sample_threshold=1,TE_threshold=3){
  print(subfamily)
  
  # Filter metadata
  metadata_matrix = filter_metadata(metadata,metric)
  metadata_table = ddply(melt(metadata_matrix[,c(1,4:9)],id.var="Sample"),.(variable,value),summarise,Num=length(Sample))
  colnames(metadata_table) = c("Category","Grouping","Category_samples")
  print("Filtered metadata")
  
  # Read in subfamily members ever in state
  subfamily_in_state = get_subfamily_in_state(subfamily,state,metric)
  print("Got individual TEs")
  
  # Add missing samples 
  subfamily_matrix = dcast(subfamily_in_state,chromosome+start+stop+subfamily+family+class+strand~Sample,value.var="Coverage")
  subfamily_matrix[,setdiff(metadata_matrix$Sample,colnames(subfamily_matrix))] = rep(NA,dim(subfamily_matrix)[1])
  subfamily_matrix[is.na(subfamily_matrix)] = 0
  subfamily_matrix = melt(subfamily_matrix,id.vars=c(TE_coordinates[c(1:4,6,5,7)]))
  colnames(subfamily_matrix)[8:9] = c("Sample","Coverage")
  print("Added missing samples")
  
  # Add metadata
  subfamily_matrix = merge(subfamily_matrix,metadata_matrix[,c("Sample","Age","Anatomy","Cancer","Germline","Group","Type")])
  subfamily_matrix = melt(subfamily_matrix,id.vars=c(TE_coordinates[c(1:4,6,5,7)],"Sample","Coverage"))
  colnames(subfamily_matrix)[10:11] = c("Category","Grouping")
  print("Added metadata")
  
  # Fraction of sample category each TE is in state
  TE_category_matrix = ddply(subfamily_matrix,.(chromosome,start,stop,subfamily,family,class,strand,Category,Grouping),here(summarise),Samples=sum(Coverage > coverage_threshold))
  print("Fraction of sample category each TE is in state")
  
  # Fraction of all samples each TE is in state
  TE_matrix = ddply(unique(subfamily_matrix[,1:9]),.(chromosome,start,stop,subfamily,family,class,strand),here(summarise),All=sum(Coverage > coverage_threshold))
  TE_category_matrix = merge(TE_category_matrix,TE_matrix,by=TE_coordinates[c(1:4,6,5,7)],all.x=TRUE)
  print("Fraction of all samples each TE is in state")
  
  # Score specificity of TE to that grouping - percent of group samples TE is in state vs all other samples
  TE_category_matrix = merge(TE_category_matrix,metadata_table,by=c("Category","Grouping"),all.x=TRUE)
  TE_category_matrix$Proportion = TE_category_matrix$Samples/TE_category_matrix$Category_samples
  TE_category_matrix$Score = TE_category_matrix$Proportion/((TE_category_matrix$All-TE_category_matrix$Samples)/(dim(metadata_matrix)[1]-TE_category_matrix$Category_samples))
  print("Calculating TE score")
  
  # Find categories with at least xx TEs with score xx
  category_matrix = ddply(TE_category_matrix,.(Category,Grouping,Category_samples),here(summarise),TEs=length(Category[which(Score > score_threshold & Proportion > category_threshold & Samples > sample_threshold)]))
  print("Finding categories over threshold")
  
  return(category_matrix[which(category_matrix$TEs >= TE_threshold),])
}

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
