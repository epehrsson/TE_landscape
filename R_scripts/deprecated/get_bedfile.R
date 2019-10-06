subfamily = "MER121"
state = "7_Enh"
metric = "chromHMM"
enrichment = THRESHOLD_LOR
category = ""
grouping = ""
  
investigate_candidate_indv = function(subfamily,state,metric,enrichment=THRESHOLD_LOR,print_fig=FALSE){
  # Get members ever in state
  subfamily_in_state = get_subfamily_in_state(subfamily,state,metric)
  
  # Samples where subfamily is enriched in state
  samples = get_enriched_samples(subfamily,state,enrichment)
  
  # Metadata for samples where subfamily is enriched
  print(metadata[which(metadata$Sample %in% samples),c(1:2,4:9)])
  
  # Plot subfamily
  if (length(samples) == 0) {
    samples = NULL
  }
  if (print_fig == "TRUE"){
    png(paste("enrichment/heatmaps/Heatmap_",subfamily,"_",state,"_",enrichment,".png",sep=""))
  }
  
  plot_binary_heatmap_indv(subfamily_in_state,metric=metric,highlight_samples=samples)
  
  if (print_fig == "TRUE"){
    dev.off()
  }
}

plot_binary_heatmap_indv = function(individual,metric="chromHMM",highlight_samples=NULL)
{
  #Input is individual TEs from subfamily ever in state, all samples
  candidate_indv = individual
  print("Got individual TEs") 
  
  #Filter metadata based on metric
  metadata_matrix = filter_metadata(metadata,metric)
  print("Filtered metadata")
  
  #Column metadata
  column_metadata = data.frame(Group=metadata_matrix$Group,Anatomy=metadata_matrix$Anatomy,Age=metadata_matrix$Age,Cancer=metadata_matrix$Cancer,Germline=metadata_matrix$Germline,Type=metadata_matrix$Type)
  column_colors = list(Age=age_colors,Cancer=cancer_colors,Germline=germline_colors,Type=type_colors,Group=group_colors,Anatomy=anatomy_colors,Class=class_colors[c(1:4,6:7)])
  if (!is.null(highlight_samples)){
    column_metadata$Enriched = factor(rep("No",dim(metadata_matrix)[1]),levels=c("No","Yes"))
    column_metadata[which(metadata_matrix$Sample %in% highlight_samples),]$Enriched = "Yes"
    column_colors$Enriched = c("white","black")
  }
  print("Assigned metadata colors")
  
  candidate_indv = dcast(candidate_indv,chromosome+start+stop+subfamily+family+class+strand~Sample,value.var="Coverage")
  rownames(candidate_indv) = apply(candidate_indv,1,function(x) paste(x[1],x[2],x[3],x[4],x[5],x[6],x[7],sep="_"))
  candidate_indv[,setdiff(metadata_matrix$Sample,colnames(candidate_indv))] = rep(NA,dim(candidate_indv)[1])
  candidate_indv = candidate_indv[,8:dim(candidate_indv)[2]]
  candidate_indv = candidate_indv[,as.vector(metadata_matrix$Sample)]
  print("Formatted matrix")
  
  #Convert NA to 0
  candidate_indv[is.na(candidate_indv)] = 0
  print("Removed NAs")
  
  # Plot
  aheatmap(candidate_indv,Rowv=FALSE,Colv=FALSE,distfun="binary",color=c("white","red"),annCol=column_metadata,annColors=column_colors,border_color="NA",annLegend=FALSE)
}

# Plot figure
investigate_candidate_indv(subfamily,state,metric,enrichment,print_fig=FALSE)

# Get members ever in state
subfamily_in_state = get_subfamily_in_state(subfamily,state,metric)
  
# Samples where subfamily is enriched in state
samples = get_enriched_samples(subfamily,state,enrichment)

# Subfamilies members in state when enriched
subfamily_enriched = get_subfamily_enriched(subfamily,state)
  
# Samples of interest
samples_filter = intersect(as.vector(metadata[which(metadata[[category]] == grouping),]$Sample),samples)
length(samples_filter)
  
# Subfamily members in the state in samples of interest, bed format with number of samples
subfamily_ubiq = aggregate(data=subfamily_in_state,Sample~chromosome+start+stop+subfamily+strand,length)
colnames(subfamily_ubiq)[6] = "Total_samples"
subfamily_bedfile = aggregate(data=subfamily_in_state[which(subfamily_in_state$Sample %in% samples_filter),],
                              Sample~chromosome+start+stop+subfamily+strand,length)
subfamily_bedfile = merge(subfamily_bedfile,subfamily_ubiq,by=c("chromosome","start","stop","subfamily","strand"))
subfamily_bedfile$Score = (subfamily_bedfile$Sample/length(samples_filter))-
  ((subfamily_bedfile$Total_samples-subfamily_bedfile$Sample)/(sample_counts["All",metric]-length(samples_filter)))
subfamily_bedfile = subfamily_bedfile[order(subfamily_bedfile$Score,subfamily_bedfile$chromosome,subfamily_bedfile$start),]

write.table(subfamily_bedfile[which(subfamily_bedfile$Sample > 1),c(1:4,6,5,7:8)],
            file=paste("enrichment/",subfamily,"_",state,"_enrich",enrichment,"_",category,"_",grouping,".bed",sep=""),
            quote=FALSE,row.names = FALSE,col.names=FALSE,sep='\t')
