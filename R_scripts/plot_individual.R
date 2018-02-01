get_subfamily_in_state = function(subfamily,state,metric){
  # Get subfamily members ever in state
  if (metric=="chromHMM"){
    subfamily_in_state = read.table(paste("chromHMM/subfamily/by_state/",subfamily,"_",state,".txt",sep=""),sep='\t')
    colnames(subfamily_in_state) = c("chromosome","start","stop","subfamily","class","family","strand","State","Overlap","Sample")
  } else if (metric=="WGBS") {
    subfamily_in_state = read.table(paste("WGBS/subfamily/by_state/",subfamily,"_",state,".txt",sep=""),sep='\t')
    colnames(subfamily_in_state) = c("chromosome","start","stop","subfamily","class","family","strand","Sample","Overlap","State")
  } else if (metric=="DNase" | metric=="H3K27ac"){
    matrix = get(paste("TE_",metric,"_overlap",sep=""))
    subfamily_in_state = melt(matrix[which(matrix$subfamily == subfamily),c(1:5,7:ncol(matrix))],
                              id.vars=c("chromosome","start","stop","subfamily","class_update","family","strand"))
    colnames(subfamily_in_state) = c("chromosome","start","stop","subfamily","class","family","strand","Sample","Overlap")
  }
  return(subfamily_in_state)
}

get_enriched_samples = function(subfamily,state,enrichment=THRESHOLD_LOR){
  samples = as.vector(unique(subfamily_state_sample_filter[which(subfamily_state_sample_filter$subfamily == subfamily 
                                                                 & subfamily_state_sample_filter$State == state
                                                                 & subfamily_state_sample_filter$Enrichment > enrichment),]$Sample))
  return(samples)
}

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
    png(paste("enrichment/heatmaps/Heatmap_",subfamily,"_",state,"_enriched",enrichment,".png",sep=""))
  }
  
  plot_binary_heatmap_indv(subfamily_in_state,metric=metric,highlight_samples=samples)
  
  if (print_fig == "TRUE"){
    dev.off()
  }
}