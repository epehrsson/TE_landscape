load("R_datasets/shuffled_DNase.RData")
load("R_datasets/shuffled_H3K27ac.RData")
load("R_datasets/shuffled_WGBS.RData")

subfamily_chromHMM_sample_expand = expand.grid(subfamily = rmsk_TE_subfamily$subfamily,Sample = metadata$Sample,State = chromHMM_states)
subfamily_DNase_expand = expand.grid(subfamily = rmsk_TE_subfamily$subfamily,Sample = as.vector(metadata[which(!is.na(metadata$DNase)),]$Sample))
subfamily_H3K27ac_expand = expand.grid(subfamily = rmsk_TE_subfamily$subfamily,Sample = as.vector(metadata[which(!is.na(metadata$H3K27ac)),]$Sample))

shuffled_enrichment = lapply(seq(1,10,1),function(i) {
  print(i)
  
  ## Number of unique chromHMM bp/peaks/CpGs in state per sample per subfamily
  print("Load Length ijk")

  # chromHMM
  print("Load chromHMM ijk")
  subfamily_chromHMM_sample_shuffle = read.table(paste("chromHMM/shuffled_TEs/subfamily/subfamily_state_sample_",i,".txt",sep=""),sep='\t')
  colnames(subfamily_chromHMM_sample_shuffle) = c("subfamily","State","Sample","Length_ijk")
  subfamily_chromHMM_sample_shuffle = merge(subfamily_chromHMM_sample_shuffle,subfamily_chromHMM_sample_expand,by=c("subfamily","Sample","State"),all=TRUE)
   
  # DNase
  print("Load DNase ijk")
  subfamily_DNase_sample_shuffle = read.table(paste("DNase/shuffled/subfamily/true_summit/subfamily_DNase_sample_summit_",i,".txt",sep=""),sep='\t')
  colnames(subfamily_DNase_sample_shuffle) = c("subfamily","Sample","Length_ijk")
  subfamily_DNase_sample_shuffle = merge(subfamily_DNase_sample_shuffle,subfamily_DNase_expand,by=c("subfamily","Sample"),all=TRUE)
  subfamily_DNase_sample_shuffle$State = rep("DNase",dim(subfamily_DNase_sample_shuffle)[1])

  # H3K27ac
  print("Load H3K27ac ijk")
  subfamily_H3K27ac_sample_shuffle = read.table(paste("H3K27ac/shuffled/subfamily/true_summit/subfamily_H3K27ac_sample_summit_",i,".txt",sep=""),sep='\t')
  colnames(subfamily_H3K27ac_sample_shuffle) = c("subfamily","Sample","Length_ijk")
  subfamily_H3K27ac_sample_shuffle = merge(subfamily_H3K27ac_sample_shuffle,subfamily_H3K27ac_expand,by=c("subfamily","Sample"),all=TRUE)
  subfamily_H3K27ac_sample_shuffle$State = rep("H3K27ac",dim(subfamily_H3K27ac_sample_shuffle)[1])

  # WGBS (CpGs)
  print("Load WGBS ijk")
  subfamily_CpG_meth_shuffle = read.table(paste("WGBS/shuffled/subfamily/subfamily_CpG_Meth_states_",i,".txt",sep=""),sep='\t')
  subfamily_CpG_meth_shuffle$V1 = mapvalues(subfamily_CpG_meth_shuffle$V1,seq(8,44,1),as.vector(metadata[which(!is.na(metadata$WGBS)),]$Sample))
  colnames(subfamily_CpG_meth_shuffle) = c("Sample",meth_states,"subfamily")
  subfamily_CpG_meth_shuffle[is.na(subfamily_CpG_meth_shuffle)] = 0
  subfamily_CpG_meth_shuffle[,meth_states] = subfamily_CpG_meth_shuffle[,meth_states]/2
  subfamily_CpG_meth_shuffle = melt(subfamily_CpG_meth_shuffle,id.vars=c("Sample","subfamily"))
  colnames(subfamily_CpG_meth_shuffle)[3:4] = c("State","Length_ijk")

  # Combine
  print("Combine ijk")
  subfamily_state_sample_shuffle = rbind(subfamily_chromHMM_sample_shuffle,subfamily_DNase_sample_shuffle,subfamily_H3K27ac_sample_shuffle,subfamily_CpG_meth_shuffle)
  rm(list=c("subfamily_chromHMM_sample_shuffle","subfamily_DNase_sample_shuffle","subfamily_H3K27ac_sample_shuffle","subfamily_CpG_meth_shuffle"))

  ## Number of subfamily members in state
  print("Load members")
  
  # chromHMM
  print("Load chromHMM members")
  subfamily_chromHMM_sample_members_shuffle = read.table(paste("chromHMM/shuffled_TEs/subfamily/subfamily_state_sample_summit_",i,".txt",sep=""),sep='\t')
  colnames(subfamily_chromHMM_sample_members_shuffle) = c("subfamily","State","Sample","Members")

  # DNase
  print("Load DNase members")

  TE_DNase_peaks_members_shuffle = melt(aggregate(data=shuffled_DNase_potential[[i]][,c(4,8:60)],.~subfamily,function(x) sum(as.numeric(x) > 0)),id.var="subfamily")
  colnames(TE_DNase_peaks_members_shuffle)[2:3] = c("Sample","Members")
  TE_DNase_peaks_members_shuffle$State = rep("DNase",dim(TE_DNase_peaks_members_shuffle)[1])
  
  # H3K27ac
  print("Load H3K27ac members")

  TE_H3K27ac_peaks_members_shuffle = melt(aggregate(data=shuffled_H3K27ac_potential[[i]][,c(4,8:105)],.~subfamily,function(x) sum(as.numeric(x) > 0)),id.var="subfamily")
  colnames(TE_H3K27ac_peaks_members_shuffle)[2:3] = c("Sample","Members")
  TE_H3K27ac_peaks_members_shuffle$State = rep("H3K27ac",dim(TE_H3K27ac_peaks_members_shuffle)[1])
  
  # WGBS (average methylation)
  print("Load WGBS members")
    
  TE_meth_average_shuffle = melt(shuffled_WGBS_average[[i]][,c(4,8:44)],id.var="subfamily")
  colnames(TE_meth_average_shuffle)[2:3] = c("Sample","Methylation")
  
  TE_meth_subfamily_shuffle = ddply(TE_meth_average_shuffle,.(subfamily,Sample),summarise,Hypomethylated=sum(na.omit(Methylation) < 0.3),
                            Intermediate=sum(na.omit(Methylation) >= 0.3 & na.omit(Methylation) <= 0.7),Hypermethylated=sum(na.omit(Methylation) > 0.7),
                            Missing=sum(is.na(Methylation)))
  TE_meth_subfamily_shuffle = melt(TE_meth_subfamily_shuffle,id.vars=c("subfamily","Sample"))
  colnames(TE_meth_subfamily_shuffle)[3:4] = c("State","Members")
  
  rm(TE_meth_average_shuffle)

  # Combine
  print("Combine members")
  subfamily_state_sample_members_shuffle = rbind(subfamily_chromHMM_sample_members_shuffle,TE_DNase_peaks_members_shuffle,TE_H3K27ac_peaks_members_shuffle,TE_meth_subfamily_shuffle)
  rm(list=c("subfamily_chromHMM_sample_members_shuffle","TE_DNase_peaks_members_shuffle","TE_H3K27ac_peaks_members_shuffle","TE_meth_subfamily_shuffle"))

  # Combine both matrices
  print("Combine ijk and members")
  subfamily_state_sample_shuffle = merge(subfamily_state_sample_shuffle,subfamily_state_sample_members_shuffle,by=c("subfamily","State","Sample"),all=TRUE)
  rm(subfamily_state_sample_members_shuffle)

  # Fill in missing values
  subfamily_state_sample_shuffle[is.na(subfamily_state_sample_shuffle)] = 0
  
  return(subfamily_state_sample_shuffle)
})

rm(list=c("shuffled_DNase_potential","shuffled_H3K27ac_potential","shuffled_WGBS_average","subfamily_chromHMM_sample_expand","subfamily_DNase_expand","subfamily_H3K27ac_expand"))

# Subfamily potential
subfamily_state_potential_shuffle = shuffled_enrichment
names(subfamily_state_potential_shuffle) = seq(1,10,1)
subfamily_state_potential_shuffle = ldply(subfamily_state_potential_shuffle,.id="Iteration")
subfamily_state_potential_shuffle = ddply(subfamily_state_potential_shuffle,.(Iteration,subfamily,State),summarise,Samples=sum(Length_ijk > 0))
subfamily_state_potential_shuffle$State = factor(subfamily_state_potential_shuffle$State,levels=states[1:21])
subfamily_state_potential_shuffle$Metric = factor(ifelse(subfamily_state_potential_shuffle$State %in% chromHMM_states,"chromHMM",
                                                         ifelse(subfamily_state_potential_shuffle$State %in% meth_states,"WGBS",
                                                                as.character(subfamily_state_potential_shuffle$State))),levels=c("chromHMM","WGBS","DNase","H3K27ac"))
subfamily_state_potential_shuffle$Sample.Proportion = as.numeric(subfamily_state_potential_shuffle$Samples/sample_counts["All",subfamily_state_potential_shuffle$Metric])
subfamily_state_potential_shuffle$Iteration = factor(subfamily_state_potential_shuffle$Iteration,levels=seq(1,10,1))

# Combine with main matrix, each iteration 
shuffled_enrichment = lapply(shuffled_enrichment,function(x) {y = merge(x,subfamily_state_sample_combined,all.y=TRUE,by=c("State","Sample","subfamily")); y[is.na(y)] <- 0; return(y)})

# Length ik for methylation states
shuffled_subfam_CpG = lapply(seq(1,10,1),function(x) read.table(paste("WGBS/shuffled/subfamily/subfamily_CpG_count_",x,".txt",sep=""),sep="\t",col.names=c("subfamily","Length_ik")))
shuffled_enrichment = lapply(seq(1,10,1),function(x) merge(shuffled_enrichment[[x]],shuffled_subfam_CpG[[x]],by="subfamily",all.x=TRUE))

# LOR Enrichment for subfamily x state x sample, using new Length ijk
shuffled_enrichment = lapply(shuffled_enrichment,function(x) transform(x,Enrichment.Shuffled = ifelse(State %in% meth_states,
                                                                                                      log2(1e-20+((x$Length_ijk.x/x$Length_jk)/(x$Length_ik.y/x$Length_k))),
                                                                                                      log2(1e-20+((x$Length_ijk.x/x$Length_jk)/(x$Length_ik.x/x$Length_k))))))

# Filter
shuffled_enrichment_filter = lapply(shuffled_enrichment,function(x) x[which(x$Members.x > THRESHOLD_IJK_MEMBER & x$Count > THRESHOLD_IK_MEMBER),])

# Count instances
subfamily_state_expand = expand.grid(subfamily = levels(rmsk_TE_subfamily$subfamily),State = states[1:21])
subfamily_state_expand = join(subfamily_state_expand,rmsk_TE_subfamily[,c("subfamily","family","class_update")],by=c("subfamily"))

shuffled_enrichment_counts = lapply(shuffled_enrichment_filter,function(x) {y = ddply(x,.(class_update,family,subfamily,State),function(x) sum(x$Enrichment.Shuffled > THRESHOLD_LOR));
y = join(subfamily_state_expand,y,by=c("subfamily","family","class_update","State"),type="left")[,c(1,3,4,2,5)];
y[is.na(y)] <- 0; return(y)})
rm(subfamily_state_expand)
names(shuffled_enrichment_counts) = seq(1,10,1)
shuffled_enrichment_counts = ldply(shuffled_enrichment_counts,.id = "Iteration")

shuffled_enrichment_counts$Sample.Proportion = ifelse(shuffled_enrichment_counts$State %in% chromHMM_states,shuffled_enrichment_counts$V1/sample_counts["All","chromHMM"],
                                                         ifelse(shuffled_enrichment_counts$State %in% meth_states,shuffled_enrichment_counts$V1/sample_counts["All","WGBS"],
                                                                ifelse(shuffled_enrichment_counts$State == "DNase",shuffled_enrichment_counts$V1/sample_counts["All","DNase"],
                                                                       ifelse(shuffled_enrichment_counts$State == "H3K27ac",shuffled_enrichment_counts$V1/sample_counts["All","H3K27ac"],"NA"))))
shuffled_enrichment_counts$Sample.Proportion = as.numeric(shuffled_enrichment_counts$Sample.Proportion)

shuffled_enrichment_counts_combine = ddply(shuffled_enrichment_counts,.(Iteration,State),summarise,Subfamilies=sum(V1 > 0),Enrichments=sum(V1))
