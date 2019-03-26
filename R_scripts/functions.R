# General
count_na = function(data_frame){
  counts = apply(data_frame,2,function(x) sum(is.na(x)))
  return(counts)
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

format_pca = function(pca,metadata,metadata_col){
  pca$eigenvectors = cbind(as.data.frame(pca$x),metadata[match(rownames(pca$x),metadata[[metadata_col]]),])
  pca$eigenvalues = 100*pca$sdev^2/sum(pca$sdev^2)
  return(pca)
}

# Plotting
add_stars = function(positions=c(),heights=c(),symbol="*",size=2,color="red"){
  star_list =  mapply(function(x, y) {paste(" + annotate(\"text\",x=",x,",y=",y,",label=\"",symbol,"\",size=",size,",color=\"",color,"\")",sep="")},
                      positions,heights)
  star_list = paste(star_list,collapse=" ")
  return(star_list)
}

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

# Matrix processing
convert_class = function(class_vector){
  class_vector = factor(class_vector,levels=unique(c(levels(class_vector),"SVA","Other")))
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
  stats = as.data.frame(t(rbind(apply(distribution,2,function(x) sum(x[2:(samples+1)])/sum(x)),
                                apply(distribution,2,function(x) sum(x*seq(0,samples))/sum(x))/samples,
                                apply(distribution,2,function(x) sd(rep(seq(0,samples),x))/sqrt(sum(x))/samples),
                                apply(as.data.frame(distribution[2:(samples+1),]),2,function(x) sum(x*seq(1,samples))/sum(x))/samples)))
  colnames(stats) = c("Proportion_ever","Samples_avg_all","Samples_SE_all","Samples_avg_ever")
  
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

get_subfamily_in_state = function(subfamily,state,metric){
  # Get subfamily members ever in state
  if (metric=="chromHMM"){
    subfamily_in_state = read.table(paste("chromHMM/subfamily/by_state/",subfamily,"_",state,".txt",sep=""),sep='\t')[,1:10]
    colnames(subfamily_in_state) = c(TE_coordinates[c(1:4,6,5,7)],"Sample","Overlap","State")
  } else if (metric=="WGBS") {
    subfamily_in_state = read.table(paste("WGBS/subfamily/by_state/",subfamily,"_",state,".txt",sep=""),sep='\t')
    colnames(subfamily_in_state) = c(TE_coordinates[c(1:4,6,5,7)],"Sample","Overlap","State")
  } else if (metric=="DNase" | metric=="H3K27ac"){
    subfamily_in_state = read.table(paste(metric,"/subfamily/",subfamily,"_",state,".txt",sep=""),sep='\t')
    colnames(subfamily_in_state) = c(TE_coordinates[c(1:4,6,5,7)],"Sample","Overlap")
  }
  print("Loaded matrix")
  
  #Calculate overlap
  if (metric == "chromHMM") {
    subfamily_in_state$Coverage = subfamily_in_state$Overlap/(subfamily_in_state$stop-subfamily_in_state$start)
  } else {
    subfamily_in_state$Coverage = subfamily_in_state$Overlap
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

get_subfamily_enriched = function(subfamily,State){
  metric = ifelse(State %in% chromHMM_states,"chromHMM",ifelse(State %in% meth_states,"WGBS",State))
  
  subfamily_in_state = get_subfamily_in_state(subfamily,State,metric)
  samples = get_enriched_samples(subfamily,State,THRESHOLD_LOR)
  subfamily_enriched = subfamily_in_state[which(subfamily_in_state$Sample %in% samples),]
  
  subfamily_ubiq = aggregate(data=subfamily_in_state,Sample~chromosome+start+stop+subfamily+strand,length)
  colnames(subfamily_ubiq)[6] = "Total_samples"
  subfamily_bedfile = aggregate(data=subfamily_enriched,Sample~chromosome+start+stop+subfamily+strand,length)
  subfamily_bedfile = merge(subfamily_bedfile,subfamily_ubiq,by=c("chromosome","start","stop","subfamily","strand"))
  subfamily_bedfile$Score = (subfamily_bedfile$Sample/length(samples))-
    ((subfamily_bedfile$Total_samples-subfamily_bedfile$Sample)/(sample_counts["All",metric]-length(samples)))
  subfamily_bedfile = subfamily_bedfile[order(subfamily_bedfile$Score,subfamily_bedfile$chromosome,subfamily_bedfile$start),]
  
  write.table(subfamily_bedfile[,c(1:4,6,5)],
              file=paste("enrichment/bedfiles/",subfamily,"_",State,"_enriched.bed",sep=""),
              quote=FALSE,row.names = FALSE,col.names=FALSE,sep='\t')
  
  return(subfamily_bedfile)
}

permute_by_sample = function(matrix,metric,direction,threshold=0,filtering,threshold2=0){ # Should be by_sample_all, split by State
  print(head(matrix,n=1))
  
  metadata_matrix = filter_metadata(metadata,ifelse(unique(matrix$State) == "H3K27ac","H3K27ac",ifelse(unique(matrix$State) == "DNase","DNase",ifelse(unique(matrix$State) %in% meth_states,"WGBS","chromHMM"))))
  metadata_table = ldply(apply(metadata_matrix[,sample_categories],2,as.data.frame(table)))
  colnames(metadata_table) = c("Category","Grouping","Samples")
  print("Computing background")
  
  # Filter matrix to those samples above thresholds
  threshold_matrix = function(matrix,metric,threshold,filtering,threshold2){
    filter_matrix = merge(matrix,metadata_matrix[,c("Sample",sample_categories)],by="Sample")
    filter_matrix = filter_matrix[which(filter_matrix[[metric]] > threshold & filter_matrix[[filtering]] > threshold2),]
    table_matrix = ldply(apply(filter_matrix[,sample_categories],2,as.data.frame(table)))
    colnames(table_matrix) = c("Category","Grouping","Samples")
    table_matrix = merge(table_matrix,metadata_table[,1:2],all.y=TRUE,by=c("Category","Grouping"))
    table_matrix[is.na(table_matrix)] = 0
    return(table_matrix)
  }
  
  real = threshold_matrix(matrix,metric,threshold,filtering,threshold2)
  print("Computing real")
  
  permuted = rdply(1000,function(x) {permute_matrix = matrix; permute_matrix$Sample = sample(permute_matrix$Sample);threshold_matrix(permute_matrix,metric,threshold,filtering,threshold2)})  
  colnames(permuted) = c("Replicate","Category","Grouping","Replicate_samples")
  print("Computing permuted")
  
  permuted = merge(permuted,real,by=c("Category","Grouping"),all.x=TRUE)
  print("Merging")
  
  if(direction == "+"){
    over = ddply(permuted,.(Category,Grouping),summarise,Pvalue=sum(Replicate_samples >= Samples)/1000)
  } else if (direction == "-") {
    over = ddply(permuted,.(Category,Grouping),summarise,Pvalue=sum(Replicate_samples <= Samples)/1000)
  }
  
  print("Computing p-values")
  
  over = merge(over,real,by=c("Category","Grouping"),all.x=TRUE)
  over = merge(over,metadata_table,by=c("Category","Grouping"),all.x=TRUE)
  return(over)
}

# Mouse
tissue_matrix = function(x,matrix){
  c(All = dim(matrix[which(matrix$Human_samples < 5 & matrix[[paste(x,".x",sep="")]] == 2),])[1],
    Specific = dim(matrix[which(matrix$Human_samples < 5 & matrix[[paste(x,".x",sep="")]] == 2 & matrix$Mouse_samples < 5 & matrix[[paste(x,".y",sep="")]] == 2),])[1],
    On = dim(matrix[which(matrix$Human_samples < 5 & matrix[[paste(x,".x",sep="")]] == 2 & matrix$Mouse_samples > 8),])[1],
    Off = dim(matrix[which(matrix$Human_samples < 5 & matrix[[paste(x,".x",sep="")]] == 2 & matrix$Mouse_samples < 2),])[1])
}

# Replicates
transform_matrix = function(matrix){
  matrix[lower.tri(matrix)] = NA
  matrix = melt(as.matrix(matrix))
  colnames(matrix) = c("Sample 1","Sample 2","Correlation")
  matrix = matrix[which(matrix$`Sample 1` != matrix$`Sample 2` & !is.na(matrix$Correlation)),]
  matrix = matrix[order(matrix$Correlation,decreasing = TRUE),]
  matrix$Rank = 1:nrow(matrix)
  return(matrix)
}

load_state = function(state){
  # Load matrix
  print("Load matrix")
  state_sample = read.table(paste("chromHMM/chromHMM_",state,".txt",sep=""),sep='\t')[,c(1:8)]
  colnames(state_sample) = c(TE_coordinates[c(1:4,6,5,7)],"Sample")
  
  # Add state total 
  print("Adding state total")
  state_sample = merge(state_sample,chromHMM_TE_state[,c(TE_coordinates,state)],by=TE_coordinates)
  colnames(state_sample)[9] = "Total"
  
  return(state_sample)
}

load_state_group = function(state){
  # Make metadata table
  print("Create metadata table")
  metadata_table = ldply(apply(metadata[,sample_categories],2,as.data.frame(table)))
  colnames(metadata_table) = c("Category","Grouping","Samples")
  
  # Load matrix
  print("Load matrix")
  state_sample = read.table(paste("chromHMM/chromHMM_",state,".txt",sep=""),sep='\t')[,c(1:8)]
  colnames(state_sample) = c(TE_coordinates[c(1:4,6,5,7)],"Sample")
  
  # Add metadata
  print("Add metadata")
  state_sample = merge(state_sample,metadata[,c("Sample",sample_categories)],by="Sample")
  state_sample = melt(state_sample,id.vars=c(TE_coordinates[c(1:4,6,5,7)],"Sample"))
  colnames(state_sample)[9:10] = c("Category","Grouping")
  
  # Count samples in Group
  print("Count by group")
  state_group = ddply(state_sample,.(chromosome,start,stop,subfamily,family,class,strand,Category,Grouping),
                      summarise,Count=length(Grouping))
  
  # Add state total 
  print("Adding state total")
  state_group = merge(state_group,chromHMM_TE_state[,c(TE_coordinates,state)],by=TE_coordinates)
  colnames(state_group)[11] = "Total"
  
  # Add Group totals
  print("Add group totals")
  state_group = merge(state_group,metadata_table,by=c("Category","Grouping"))
  colnames(state_group)[12] = "Category.Total"
  
  # Save data
  save(state_group,file=paste("chromHMM/chromHMM_",state,".RData",sep=""))
  
  return(state_group)
}