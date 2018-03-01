kruskal_if = function(x,variableA,variableB){ #Where x is data frame subset, variableA are values to be split, variableB is column to be tested
  if(length(unique(x[[variableB]])) > 1){
    p_value = unlist(kruskal.test(x[[variableA]] ~ x[[variableB]]))["p.value"]
  }
  else{
    p_value = NA
  }
  return(p_value)
}

wilcox_to_all = function(all,metadata){ #all and metadata are vectors, e.g., columns of data frames
  groups = levels(metadata)
  tests = p.adjust(as.numeric(unlist(lapply(groups,function(x) unlist(wilcox.test(all[which(metadata == x)],all[which(metadata != x)]))["p.value"]))),method="bonferroni")
  names(tests) = groups
  return(tests)
}

convert_class = function(class_vector){
  class_vector = factor(class_vector,levels=unique(c(levels(class_vector),"SVA")))
  class_vector[which(class_vector == "Other")] = "SVA"
  class_vector[which(class_vector %in% c("DNA?","LINE?","LTR?","SINE?","Unknown","Unknown?","RC","Unconfident","Unconfident_RC"))] = "Other"
  class_vector = factor(class_vector,levels=c("DNA","LINE","LTR","SINE","SVA","Other"))
  return(class_vector)
}

split_coding = function(feature_matrix,new_features=c()){
  to_split = gsub("coding_exon","CDS",feature_matrix$Cohort)
  feature_matrix$Feature = factor(unlist(lapply(to_split,function(x) unlist(strsplit(as.character(x),"_"))[1])),levels=c(new_features,features))
  feature_matrix$Coding = unlist(lapply(to_split,function(x) unlist(strsplit(as.character(x),"_"))[2]))
  feature_matrix[which(is.na(feature_matrix$Coding)),]$Coding = "All"
  feature_matrix$Coding = factor(feature_matrix$Coding,levels=c("All","pc","nc"))
  feature_matrix$Cohort = factor(feature_matrix$Cohort,levels=c(new_features,cohorts))
  
  return(feature_matrix)
}

# Potential
sample_distribution = function(TE_matrix,columns,samples){
  dist_TE = data.frame(Samples = seq(0,samples))
  for (i in columns){
    dist_TE = merge(dist_TE,as.data.frame(table(as.integer(TE_matrix[,i]))),by.x="Samples",by.y="Var1",all=TRUE)
    colnames(dist_TE)[length(colnames(dist_TE))] = colnames(TE_matrix)[i]
  }
  dist_TE[is.na(dist_TE)] = 0
  
  return(dist_TE)
}

cumulative_distribution = function(TE_matrix,columns,samples){ #TE_matrix is matrix of TEs x state listing number of samples, columns is state columns, samples is number of samples
  #Generating distribution matrix
  dist_TE = sample_distribution(TE_matrix,columns,samples)
  
  #Generating cumulative distribution matrix
  if (length(columns) > 1){
    dist_TE = dist_TE[which(dist_TE$Samples %in% seq(1,samples)),2:(length(columns)+1)]
  }
  else{
    dist_TE = as.data.frame(dist_TE[which(dist_TE$Samples %in% seq(1,samples)),2])
    colnames(dist_TE) = "State"
  }
  cum_dist_TE = as.data.frame(apply(dist_TE,2,function(x) cumsum(x)/sum(x)))
  rownames(cum_dist_TE) = seq(1:samples)
  cum_dist_TE = melt(as.matrix(cum_dist_TE))
  colnames(cum_dist_TE) = c("Samples","State","Proportion")
  cum_dist_TE$Sample_proportion = as.numeric(cum_dist_TE$Samples)/samples
  
  return(cum_dist_TE)
}

potential_stats = function(distribution,states,samples){
  if (states > 1){
    distribution = distribution[,2:(states+1)]
  }
  else{
    distribution = as.data.frame(distribution[,2])
  }
  stats = as.data.frame(t(rbind(apply(distribution,2,function(x) sum(x[2:(samples+1)])/(sum(x)/100)),
                                apply(distribution,2,function(x) sum(x*seq(0,samples))/sum(x))/(samples/100),
                                apply(as.data.frame(distribution[2:(samples+1),]),2,function(x) sum(x*seq(1,samples))/sum(x))/(samples/100))))
  colnames(stats) = c("Proportion_ever","Samples_avg_all","Samples_avg_ever")
  
  return(stats)
}

# TE features
correlate_spearman = function(matrix, indpt_var, response_vars){ #For all TEs, correlate TE features with number of samples in state
  # All TEs
  indv = as.data.frame(t(apply(matrix[,response_vars],2,function(x) unlist(cor.test(matrix[,indpt_var],x,method="spearman"))[c("p.value","estimate.rho")])))
  indv$State = rownames(indv)
  indv$class_update = rep("All",dim(indv)[1])
  
  # By class
  class = merge(melt(ddply(matrix,~class_update,function(y) apply(y[,response_vars],2,function(x) unlist(cor.test(y[,indpt_var],x,method="spearman"))["p.value"])),id.vars="class_update"),
                melt(ddply(matrix,~class_update,function(y) apply(y[,response_vars],2,function(x) unlist(cor.test(y[,indpt_var],x,method="spearman"))["estimate.rho"])),id.vars="class_update"),
                by=c("class_update","variable"))
  colnames(class) = c("class_update","State","p.value","estimate.rho")
  
  class = rbind(indv,class)
  return(class)
}

# Enrichment
enrichment_proportion = function(matrix,enrichment,threshold,metric,members_threshold=0){
  #Cannot process more than one state at a time
  metadata_matrix = filter_metadata(metadata,metric)

  categories = c("Age","Anatomy","Cancer","Germline","Group","Type")
  proportions = list()
  
  for (i in 1:6) {
    filtered_matrix = matrix[,c("subfamily",enrichment,categories[i],"Members")]
    colnames(filtered_matrix) = c("subfamily","Enrichment","Category","Members")
    aggregate_matrix = ddply(filtered_matrix,.(subfamily,Category),function(x) dim(x[which(x$Enrichment > threshold & x$Members >= members_threshold),])[1])
    colnames(aggregate_matrix)[3] = c("Enriched")
    aggregate_matrix$Metadata = rep(categories[i],dim(aggregate_matrix)[1])
    proportions[[i]] = aggregate_matrix
  }
  proportions = ldply(proportions)
  proportions$Proportion = apply(proportions,1,function(x) as.numeric(x[3])/length(metadata_matrix[which(metadata_matrix[,x[4]] == x[2]),x[4]]))
  
  return(proportions)
}

subfamily_member_to_matrix = function(in_list,state){
  new_matrix = dcast(in_list[which(in_list$State == state),],Class + Family + Subfamily ~ Sample)
  new_matrix[is.na(new_matrix)] = 0
  return(new_matrix)
}

# Plotting
get_legend = function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

gg_color_hue = function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

heatmap_color = function (n) {
  x <- ramp(seq.int(0, 1, length.out = n))
  if (ncol(x) == 4L) 
    rgb(x[, 1L], x[, 2L], x[, 3L], x[, 4L], maxColorValue = 255)
  else rgb(x[, 1L], x[, 2L], x[, 3L], maxColorValue = 255)
}

plot_pca = function(pca,axes,map,colorby,legend_title,level_colors,guide=TRUE){
  eigenvectors = as.data.frame(pca$x)
  eigenvalues = 100*pca$sdev^2/sum(pca$sdev^2)
  point_color = map[,colorby]
  #print(levels(point_color))
  colors = level_colors[levels(point_color)]
  #print(colors)
  
  if(guide == TRUE){
    ggplot(eigenvectors,aes(x=eigenvectors[,axes[1]],y=eigenvectors[,axes[2]]),environment = environment()) + geom_point(aes(color = factor(point_color)),size=4) + labs(x=paste("PC",axes[1]," (",round(eigenvalues[axes[1]],1),"%)",sep=""),y=paste("PC",axes[2]," (",round(eigenvalues[axes[2]],1),"%)",sep="")) + theme_classic() + theme(text=element_text(face="bold"),panel.background = element_rect(fill=NA,color="black"),aspect.ratio = 1) + scale_colour_manual(name=legend_title,breaks=levels(point_color),labels=names(colors),values=colors) + guides(colour = guide_legend(title.hjust=0.12))
  }
  else{
    ggplot(eigenvectors,aes(x=eigenvectors[,axes[1]],y=eigenvectors[,axes[2]]),environment = environment()) + geom_point(aes(color = factor(point_color)),size=4) + labs(x=paste("PC",axes[1]," (",round(eigenvalues[axes[1]],1),"%)",sep=""),y=paste("PC",axes[2]," (",round(eigenvalues[axes[2]],1),"%)",sep="")) + theme_classic() + theme(text=element_text(face="bold"),panel.background = element_rect(fill=NA,color="black"),aspect.ratio = 1) + scale_colour_manual(breaks=levels(point_color),labels=names(colors),values=colors,guide=FALSE) 
  }
}

filter_metadata = function(metadata,metric){
  if (metric == "WGBS") {
    metadata_matrix = metadata[which(!is.na(metadata$WGBS)),]
  } else if (metric == "DNase") {
    metadata_matrix = metadata[which(!is.na(metadata$DNase)),]
  } else if (metric == "H3K27ac") {
    metadata_matrix = metadata[which(!is.na(metadata$H3K27ac)),]
  } else{
    metadata_matrix = metadata
  }
}

plot_binary_heatmap = function(matrix,metric="chromHMM",state="none",min_sample=0,max_sample=128,enrichment_threshold=1.5,enrichment_column="Enrichment",subfamilies=NULL) 
{
  #Input matrix should be subfamily_state_sample_filter
  
  #Filter metadata based on metric
  metadata_matrix = filter_metadata(metadata,metric)
  
  #Column metadata
  column_metadata = data.frame(Group=metadata_matrix$Group,Anatomy=metadata_matrix$Anatomy,Age=metadata_matrix$Age,Cancer=metadata_matrix$Cancer,Germline=metadata_matrix$Germline,Type=metadata_matrix$Type)
  column_colors = list(Age=brewer.pal(4,"YlOrRd"),Cancer=c("white","red"),Germline=brewer.pal(6,"Dark2"),Type=brewer.pal(5,"Greens"),Group=group_colors,Anatomy=anatomy_colors,Class=class_colors[c(1:4,6:7)])
  
  #Create matrix
  if (state != "none"){
    matrix = matrix[which(matrix$State == state),] #Filter by state
  } 
  
  test = dcast(matrix,subfamily~Sample,value.var=enrichment_column)
  rownames(test) = test[,1]
  test[,setdiff(metadata_matrix$Sample,colnames(test))] = rep(NA,dim(test)[1])
  test = test[,2:dim(test)[2]]
  test = test[,as.vector(metadata_matrix$Sample)]
  
  #Convert to binary
  test[is.na(test) | test < enrichment_threshold] = 0
  test[test > enrichment_threshold] = 1
  
  #Filter by number of samples
  test = test[which(apply(test,1,sum) >= min_sample & apply(test,1,sum) <= max_sample),]
  
  #Filter by subfamilies of interest
  if (!is.null(subfamilies)){ 
    test = test[subfamilies,]
  }

  #Plot
  aheatmap(test,Rowv=FALSE,Colv=FALSE,distfun="binary",breaks=0.5,legend=FALSE,color=c("white","cornflowerblue"),border_color="NA",
           annRow=data.frame(Class=rmsk_TE_subfamily[match(rownames(test),rmsk_TE_subfamily$subfamily),]$class_update),annColors=column_colors,annCol=column_metadata,annLegend=FALSE)
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
  column_colors = list(Age=brewer.pal(4,"YlOrRd"),Cancer=c("white","red"),Germline=brewer.pal(6,"Dark2"),Type=brewer.pal(5,"Greens"),Group=group_colors,Anatomy=anatomy_colors,Class=class_colors[c(1:4,6:7)])
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
  if (state %in% chromHMM_states){
    test = merge(test,rmsk_TE_subfamily_ever[,c("subfamily",stateX)],by="subfamily")
  } else {
    test = merge(test,rmsk_TE_subfamily_ever[,c("subfamily",state)],by="subfamily")
  }
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

print_individual_TEs = function(subfamily_state_sample){
  input_matrix = read.table(file=subfamily_state_sample,sep='\t')
  test = list()
  for (i in 1:dim(input_matrix)[1]){
    subfamily = as.character(input_matrix[i,1])
    sample = as.character(input_matrix[i,2])
    state = as.character(input_matrix[i,3])
    if (state == "DNase"){
      matrix = TE_DNase_peaks
    } else if (state == "H3K27ac"){
      matrix = TE_H3K27ac_peaks
    } else {
      print("Error: not DNase or H3K27ac")
    }
    test[[i]] = matrix[which(matrix$subfamily == subfamily & matrix[[sample]] > 0),c(colnames(matrix)[1:7],sample)]
    test[[i]] = melt(test[[i]],id.vars=c("chromosome","start","stop","subfamily","family","class","strand"))
  }
  test = ldply(test)
  test_instance = aggregate(data=test,variable~chromosome+start+stop+subfamily+family+class+strand,length)
  lapply(unique(test_instance$subfamily),function(x) write.table(test_instance[which(test_instance$subfamily == x),c(1:4,8,7)],file=paste("enrichment/",x,"_",state,"_enriched.bed",sep=""),sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE))
}

get_subfamily_in_state = function(subfamily,state,metric){
  # Get subfamily members ever in state
  if (metric=="chromHMM"){
    subfamily_in_state = read.table(paste("chromHMM/subfamily/by_state/",subfamily,"_",state,".txt",sep=""),sep='\t')
    colnames(subfamily_in_state) = c("chromosome","start","stop","subfamily","class","family","strand","State","Overlap","Sample")
  } else if (metric=="WGBS") {
    subfamily_in_state = read.table(paste("WGBS/subfamily/by_state/",subfamily,"_",state,".txt",sep=""),sep='\t')
    colnames(subfamily_in_state) = c("chromosome","start","stop","subfamily","class","family","strand","Sample","Overlap")
  } else if (metric=="DNase" | metric=="H3K27ac"){
    subfamily_in_state = read.table(paste(metric,"/subfamily/",subfamily,"_",state,".txt",sep=""),sep='\t')
    colnames(subfamily_in_state) = c("chromosome","start","stop","subfamily","class","family","strand","Sample","Overlap")
  }
  print("Loaded matrix")
  
  #Calculate overlap
  if (metric == "WGBS") {
    subfamily_in_state = merge(subfamily_in_state,TE_meth_average[,c("chromosome","start","stop","subfamily","class","family","strand","CpGs")],by=c("chromosome","start","stop","subfamily","class","family","strand"))
    subfamily_in_state$Coverage = subfamily_in_state$Overlap/subfamily_in_state$CpGs
  } else {
    subfamily_in_state$Coverage = subfamily_in_state$Overlap/(subfamily_in_state$stop-subfamily_in_state$start)
  }
  print("Calculated overlap")
  
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
    png(paste("enrichment/heatmaps/Heatmap_",subfamily,"_",state,"_",enrichment,".png",sep=""))
  }
  
  plot_binary_heatmap_indv(subfamily_in_state,metric=metric,highlight_samples=samples)
  
  if (print_fig == "TRUE"){
    dev.off()
  }
}

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
  subfamily_matrix = melt(subfamily_matrix,id.vars=c("chromosome","start","stop","subfamily","class","family","strand"))
  colnames(subfamily_matrix)[8:9] = c("Sample","Coverage")
  print("Added missing samples")
  
  # Add metadata
  subfamily_matrix = merge(subfamily_matrix,metadata_matrix[,c("Sample","Age","Anatomy","Cancer","Germline","Group","Type")])
  subfamily_matrix = melt(subfamily_matrix,id.vars=c("chromosome","start","stop","subfamily","class","family","strand","Sample","Coverage"))
  colnames(subfamily_matrix)[10:11] = c("Category","Grouping")
  print("Added metadata")
  
  # Fraction of sample category each TE is in state
  TE_category_matrix = ddply(subfamily_matrix,.(chromosome,start,stop,subfamily,family,class,strand,Category,Grouping),here(summarise),Samples=sum(Coverage > coverage_threshold))
  print("Fraction of sample category each TE is in state")
  
  # Fraction of all samples each TE is in state
  TE_matrix = ddply(unique(subfamily_matrix[,1:9]),.(chromosome,start,stop,subfamily,family,class,strand),here(summarise),All=sum(Coverage > coverage_threshold))
  TE_category_matrix = merge(TE_category_matrix,TE_matrix,by=c("chromosome","start","stop","subfamily","class","family","strand"),all.x=TRUE)
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