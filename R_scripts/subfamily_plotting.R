run_tsne = function(matrix, perplex=30, the = 0){
  tsne = Rtsne(matrix,perplexity = perplex,theta= the,max_iter = 5000,pca_scale=TRUE)
  tsne_plot = as.data.frame(tsne$Y)
  tsne_plot$object = rownames(matrix)
  return(tsne_plot)
}

plot_biplot = function(pca,axes=c(1,2),metric,category,colors,variables=NULL,guide=TRUE){
  # Filter metadata
  metadata_matrix = filter_metadata(metadata,metric)

  # Get loadings 
  loadings = as.data.frame(pca$rotation %*% diag(pca$sdev))
  loadings$variable = rownames(loadings)
  
  # Get scores scaled to unit variance
  scores = as.data.frame(pca$x %*% diag(1/pca$sdev))
  scores$observations = rownames(scores)
  
  # Get variance
  variance = 100*pca$sdev^2/sum(pca$sdev^2)
  
  # Calculate scaling factor (from biplot)
  unsigned.range = function(x) c(-abs(min(x)),abs(max(x)))
  scale = 1/(max(unsigned.range(loadings[,axes[1]])/unsigned.range(scores[,axes[1]]),
              unsigned.range(loadings[,axes[2]])/unsigned.range(scores[,axes[2]])))
  print(scale)
  
  # Filter loadings to variables of interest
  loadings_filter = loadings[which(loadings$variable %in% variables),]
  
  # Plot biplot 
  ggplot(scores,aes(x=scores[,axes[1]],y=scores[,axes[2]],color = metadata_matrix[,category]),environment = environment()) +
    geom_point() + theme(text=element_text(face="bold"),aspect.ratio = 1) +
    labs(x=paste("PC",axes[1]," (",round(variance[axes[1]],1),"%)",sep=""),
         y=paste("PC",axes[2]," (",round(variance[axes[2]],1),"%)",sep="")) +
    scale_colour_manual(name=category,values=colors) + 
    geom_segment(data=loadings_filter,aes(x=0,y=0,xend=loadings_filter[,axes[1]]*scale*0.95,yend=loadings_filter[,axes[2]]*scale*0.95),color="red",arrow = arrow(length = unit(0.01, "npc"))) + 
    geom_text(data=loadings_filter,aes(x=loadings_filter[,axes[1]]*scale,y=loadings_filter[,axes[2]]*scale,label=variable),color="red",size=3)
}

write_subfamily_candidates = function(candidate_list,state,print_coords=TRUE){
  # See 7/11/2016, 8/29/2016, 9/29/2016, 11/22/2016, 11/27/2016, 5/23/2017, 5/24/2017, 5/29/2017, 5/30/2017, 6/15/2017
  
  stateX = paste("X",state,sep="")
  if (state == "8_ZNF/Rpts"){
    print_state = "8_ZNF.Rpts"
  } else {
    print_state = state
  }
  
  # Statistics for candidate subfamilies
  test = ddply(subfamily_state_sample_filter[which(subfamily_state_sample_filter$subfamily %in% candidate_list & subfamily_state_sample_filter$Enrichment > THRESHOLD_LOR & subfamily_state_sample_filter$State == state),],.(subfamily),summarise,Min = min(Members),Max = max(Members),Median = median(Members))
  test = merge(test,rmsk_TE_subfamily_ever[,c("subfamily",state)],by="subfamily")
  colnames(test)[5] = "Members_ever"
  test = merge(test,subfamily_state_sample_counts[which(subfamily_state_sample_counts$State == state),c(1:3,5)],by="subfamily")
  colnames(test)[8] = "Samples_enriched"
  test = merge(test,rmsk_TE_subfamily[,c(1,4)],by="subfamily")
  colnames(test)[9] = "Members"
  write.table(test[,c(1,6:7,8,9,5,2:4)],row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t',file=paste("enrichment/candidate_",print_state,"_stat.txt",sep=""))
  
  # Write enriched subfamily coordinates
  if (print_coords == TRUE){
    if (state %in% chromHMM_states){
      # TEs ever in state, by subfamily	 
      lapply(candidate_list,function(x) write.table(rmsk_TE_measure[which(rmsk_TE_measure$subfamily == x & rmsk_TE_measure[[stateX]] > 0),c(colnames(rmsk_TE_measure)[1:4],stateX,"strand")],sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE,file=paste("enrichment/",x,"_",print_state,".bed",sep="")))
      # TEs never in state, by subfamily	 
      lapply(candidate_list,function(x) write.table(rmsk_TE_measure[which(rmsk_TE_measure$subfamily == x & rmsk_TE_measure[[stateX]] == 0),c(colnames(rmsk_TE_measure)[1:4],stateX,"strand")],sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE,file=paste("enrichment/",x,"_no",print_state,".bed",sep="")))
    } else {
      # TEs ever in state, by subfamily	 
      lapply(candidate_list,function(x) write.table(rmsk_TE_measure[which(rmsk_TE_measure$subfamily == x & rmsk_TE_measure[[state]] > 0),c(colnames(rmsk_TE_measure)[1:4],state,"strand")],sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE,file=paste("enrichment/",x,"_",print_state,".bed",sep="")))
      # TEs never in state, by subfamily	 
      lapply(candidate_list,function(x) write.table(rmsk_TE_measure[which(rmsk_TE_measure$subfamily == x & rmsk_TE_measure[[state]] == 0),c(colnames(rmsk_TE_measure)[1:4],state,"strand")],sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE,file=paste("enrichment/",x,"_no",print_state,".bed",sep="")))
    }
  }
  
  # Write samples where candidate subfamilies are enriched in state	 
  write.table(subfamily_state_sample_filter[which(subfamily_state_sample_filter$Enrichment > THRESHOLD_LOR & subfamily_state_sample_filter$State == state & subfamily_state_sample_filter$subfamily %in% candidate_list),c("subfamily","Sample","State")],row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t',file=paste("enrichment/candidate_",print_state,"_enriched.txt",sep=""))
}

enrichment_proportion = function(matrix,enrichment,threshold,metric,members_threshold=0){
  #Cannot process more than one state at a time
  metadata_matrix = filter_metadata(metadata,metric)
  
  proportions = list()
  
  for (i in 1:6) {
    filtered_matrix = matrix[,c("subfamily",enrichment,sample_categories[i],"Members")]
    colnames(filtered_matrix) = c("subfamily","Enrichment","Category","Members")
    aggregate_matrix = ddply(filtered_matrix,.(subfamily,Category),function(x) dim(x[which(x$Enrichment > threshold & x$Members >= members_threshold),])[1])
    colnames(aggregate_matrix)[3] = c("Enriched")
    aggregate_matrix$Metadata = rep(sample_categories[i],dim(aggregate_matrix)[1])
    proportions[[i]] = aggregate_matrix
  }
  proportions = ldply(proportions)
  proportions$Proportion = apply(proportions,1,function(x) as.numeric(x[3])/length(metadata_matrix[which(metadata_matrix[,x[4]] == x[2]),x[4]]))
  
  return(proportions)
}


