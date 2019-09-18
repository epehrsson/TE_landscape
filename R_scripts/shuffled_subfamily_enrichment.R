# Calculates the log odds ratio enrichment of each shuffled TE subfamily in each epigenetic state in each sample
# then creates a table of the number of samples each subfamily is enriched LOR > 1.5 and the total number of enrichments, 
# for 10 iterations of shuffled TEs.
# Includes thresholds for considering a subfamily x state x sample combination

## subfamily_state_potential_shuffle - Number/proportion of samples in which each TE subfamily is annotated with each epigenetic state, 
## for each iteration of shuffled TEs
## shuffled_enrichment - LOR enrichment for each subfamily x state x sample combination, for each iteration of shuffled TEs
## shuffled_enrichment_filter - LOR enrichment for each subfamily x state x sample combination, for each iteration of shuffled TEs,
## limited to those that pass thresholds for the number of members overall and in the state
## shuffled_enrichment_counts - Number/proportion of samples with LOR enrichment > 1.5 per subfamily x state, for each iteration of shuffled TEs,
## including only those that pass subfamily member thresholds
## shuffled_enrichment_counts_combine - Number of subfamilies with at least one enrichment > 1.5 and total number of enrichments LOR > 1.5 per state, 
## for each iteration of shuffled TEs


# Load matrices of the epigenetic state of shuffled TEs in each sample, including
# WGBS: average methylation and DHS/H3K27ac: overlap with peak summits

load("R_datasets/shuffled_DNase.RData")
load("R_datasets/shuffled_H3K27ac.RData")
load("R_datasets/shuffled_WGBS.RData")

# Creates dataframes with all possible combinations of subfamily x sample x state for samples with chromHMM, DHS, and H3K27ac data, respectively
subfamily_chromHMM_sample_expand = expand.grid(subfamily = rmsk_TE_subfamily$subfamily,Sample = metadata$Sample,State = chromHMM_states)
subfamily_DNase_expand = expand.grid(subfamily = rmsk_TE_subfamily$subfamily,Sample = as.vector(metadata[which(!is.na(metadata$DNase)),]$Sample))
subfamily_H3K27ac_expand = expand.grid(subfamily = rmsk_TE_subfamily$subfamily,Sample = as.vector(metadata[which(!is.na(metadata$H3K27ac)),]$Sample))

# For each iteration of shuffled TEs, the length of each state in each subfamily per sample (ijk)
# and the number of subfamily members in each state per sample
shuffled_enrichment = lapply(seq(1,10,1),function(i) {
  print(i)
  
  # Length of each state in each subfamily in each sample (ijk)
  print("Load Length ijk")

  ## chromHMM: length of overlap (bp) between subfamily and chromHMM state per sample
  ## Including all possible combinations of subfamily x state x sample
  print("Load chromHMM ijk")
  subfamily_chromHMM_sample_shuffle = read.table(paste("chromHMM/shuffled_TEs/subfamily/subfamily_state_sample_",i,".txt",sep=""),sep='\t',
                                                 col.names=c("subfamily","State","Sample","Length_ijk"))
  subfamily_chromHMM_sample_shuffle = merge(subfamily_chromHMM_sample_shuffle,subfamily_chromHMM_sample_expand,by=c("subfamily","Sample","State"),all=TRUE)
   
  ## DHS: number of unique peaks whose summit overlaps each subfamily per sample
  ## Including all possible combinations of subfamily x state x sample
  print("Load DNase ijk")
  subfamily_DNase_sample_shuffle = read.table(paste("DNase/shuffled/subfamily/true_summit/subfamily_DNase_sample_summit_",i,".txt",sep=""),sep='\t',
                                              col.names=c("subfamily","Sample","Length_ijk"))
  subfamily_DNase_sample_shuffle = merge(subfamily_DNase_sample_shuffle,subfamily_DNase_expand,by=c("subfamily","Sample"),all=TRUE)
  subfamily_DNase_sample_shuffle$State = rep("DNase",dim(subfamily_DNase_sample_shuffle)[1])

  ## H3K27ac: number of unique peaks whose summit overlaps each subfamily per sample
  ## Including all possible combinations of subfamily x state x sample
  print("Load H3K27ac ijk")
  subfamily_H3K27ac_sample_shuffle = read.table(paste("H3K27ac/shuffled/subfamily/true_summit/subfamily_H3K27ac_sample_summit_",i,".txt",sep=""),sep='\t',
                                                col.names=c("subfamily","Sample","Length_ijk"))
  subfamily_H3K27ac_sample_shuffle = merge(subfamily_H3K27ac_sample_shuffle,subfamily_H3K27ac_expand,by=c("subfamily","Sample"),all=TRUE)
  subfamily_H3K27ac_sample_shuffle$State = rep("H3K27ac",dim(subfamily_H3K27ac_sample_shuffle)[1])

  ## WGBS: Number of unique CpGs overlapping each subfamily in each methylation state per sample
  print("Load WGBS ijk")
  subfamily_CpG_meth_shuffle = read.table(paste("WGBS/shuffled/subfamily/subfamily_CpG_Meth_states_",i,".txt",sep=""),sep='\t')
  ## Add sample names
  subfamily_CpG_meth_shuffle$V1 = mapvalues(subfamily_CpG_meth_shuffle$V1,seq(8,44,1),as.vector(metadata[which(!is.na(metadata$WGBS)),]$Sample))
  colnames(subfamily_CpG_meth_shuffle) = c("Sample",meth_states,"subfamily")
  subfamily_CpG_meth_shuffle[is.na(subfamily_CpG_meth_shuffle)] = 0
  subfamily_CpG_meth_shuffle[,meth_states] = subfamily_CpG_meth_shuffle[,meth_states]/2
  subfamily_CpG_meth_shuffle = melt(subfamily_CpG_meth_shuffle,id.vars=c("Sample","subfamily"),
                                    variable.name="State",value.name="Length_ijk")
  
  # Combine into a single dataframe
  print("Combine ijk")
  subfamily_state_sample_shuffle = rbind(subfamily_chromHMM_sample_shuffle,subfamily_DNase_sample_shuffle,subfamily_H3K27ac_sample_shuffle,subfamily_CpG_meth_shuffle)
  rm(list=c("subfamily_chromHMM_sample_shuffle","subfamily_DNase_sample_shuffle","subfamily_H3K27ac_sample_shuffle","subfamily_CpG_meth_shuffle"))

  # Number of subfamily members in each epigenetic state per sample
  print("Load members")
  
  ## chromHMM: number of subfamily members where 1) the TE overlaps the center of a 200bp chromHMM bin annotated with the state
  ## or 2) the state covers the majority of the TE, for TEs that do not overlap the center of any 200bp chromHMM bin,
  ## by chromHMM state and sample
  print("Load chromHMM members")
  subfamily_chromHMM_sample_members_shuffle = read.table(paste("chromHMM/shuffled_TEs/subfamily/subfamily_state_sample_summit_",i,".txt",sep=""),sep='\t',
                                                         col.names=c("subfamily","State","Sample","Members"))

  ## DHS: number of subfamily members overlapping the summit of a DHS peak per sample (column)
  print("Load DNase members")

  TE_DNase_peaks_members_shuffle = melt(aggregate(data=shuffled_DNase_potential[[i]][,c(4,8:60)],.~subfamily,function(x) sum(as.numeric(x) > 0)),
                                        id.var="subfamily",variable.name="Sample",value.name="Members")
  TE_DNase_peaks_members_shuffle$State = rep("DNase",dim(TE_DNase_peaks_members_shuffle)[1])
  
  ## H3K27ac: number of subfamily members overlapping the summit of an H3K27ac peak per sample (column)
  print("Load H3K27ac members")

  TE_H3K27ac_peaks_members_shuffle = melt(aggregate(data=shuffled_H3K27ac_potential[[i]][,c(4,8:105)],.~subfamily,function(x) sum(as.numeric(x) > 0)),
                                          id.var="subfamily",variable.name="Sample",value.name="Members")
  TE_H3K27ac_peaks_members_shuffle$State = rep("H3K27ac",dim(TE_H3K27ac_peaks_members_shuffle)[1])
  
  ## WGBS: Number of subfamily members in each methylation state per sample, based on TE average methylation
  print("Load WGBS members")
    
  TE_meth_average_shuffle = melt(shuffled_WGBS_average[[i]][,c(4,8:44)],
                                 id.var="subfamily",variable.name="Sample",value.name="Methylation")
  
  TE_meth_subfamily_shuffle = ddply(TE_meth_average_shuffle,.(subfamily,Sample),summarise,Hypomethylated=sum(na.omit(Methylation) < 0.3),
                            Intermediate=sum(na.omit(Methylation) >= 0.3 & na.omit(Methylation) <= 0.7),Hypermethylated=sum(na.omit(Methylation) > 0.7),
                            Missing=sum(is.na(Methylation)))
  TE_meth_subfamily_shuffle = melt(TE_meth_subfamily_shuffle,id.vars=c("subfamily","Sample"),
                                   variable.name="State",value.name="Members")
  
  rm(TE_meth_average_shuffle)

  # Combine into a single dataframe
  print("Combine members")
  subfamily_state_sample_members_shuffle = rbind(subfamily_chromHMM_sample_members_shuffle,TE_DNase_peaks_members_shuffle,TE_H3K27ac_peaks_members_shuffle,TE_meth_subfamily_shuffle)
  rm(list=c("subfamily_chromHMM_sample_members_shuffle","TE_DNase_peaks_members_shuffle","TE_H3K27ac_peaks_members_shuffle","TE_meth_subfamily_shuffle"))

  # Combine Length ijk and number of members in the state for each subfamily x state x sample combination
  print("Combine ijk and members")
  subfamily_state_sample_shuffle = merge(subfamily_state_sample_shuffle,subfamily_state_sample_members_shuffle,by=c("subfamily","State","Sample"),all=TRUE)
  rm(subfamily_state_sample_members_shuffle)

  # Fill in Length ijk and members for missing combinations
  subfamily_state_sample_shuffle[is.na(subfamily_state_sample_shuffle)] = 0
  
  return(subfamily_state_sample_shuffle)
})

rm(list=c("shuffled_DNase_potential","shuffled_H3K27ac_potential","shuffled_WGBS_average","subfamily_chromHMM_sample_expand","subfamily_DNase_expand","subfamily_H3K27ac_expand"))


# Number/proportion of samples in which each TE subfamily is annotated with each epigenetic state, for each iteration of shuffled TEs
## chromHMM: length of overlap >= 1bp, WGBS: number of CpGs >= 1, DHS/H3K27ac: number of peak summits >= 1
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


# Combine the dataframes of Length ijk for each subfamily x state x sample combination for each iteration of shuffled TEs 
# with the dataframe of LOR enrichment for each real subfamily x state x sample combination
shuffled_enrichment = lapply(shuffled_enrichment,function(x) {y = merge(x,subfamily_state_sample_combined,all.y=TRUE,by=c("State","Sample","subfamily")); y[is.na(y)] <- 0; return(y)})

# Add the length of the shuffled subfamilies (Length ik for all other states)
shuffled_ik = lapply(seq(1,10,1),function(x) read.table(paste("features/shuffled_TEs/subfamily/subfamily_merge_",x,"_length.txt",sep=""),sep='\t',
                                                        col.names=c("subfamily","Total_length","Total_length_noY")))
shuffled_enrichment = lapply(seq(1,10,1),function(x) {y <- merge(shuffled_enrichment[[x]],shuffled_ik[[x]],by="subfamily",all.x=TRUE); y[is.na(y)] = 0; return(y)})
rm(shuffled_ik)

# Add the number of CpGs per shuffled subfamily (Length ik for methylation states)
shuffled_subfam_CpG = lapply(seq(1,10,1),function(x) {y <-read.table(paste("WGBS/shuffled/subfamily/subfamily_CpG_count_",x,".txt",sep=""),sep="\t",
                                                                col.names=c("subfamily","CpGs")); y$CpGs = y$CpGs/2; return(y)})
shuffled_enrichment = lapply(seq(1,10,1),function(x) {y <- merge(shuffled_enrichment[[x]],shuffled_subfam_CpG[[x]],by="subfamily",all.x=TRUE); y[is.na(y)] = 0; return(y)})
rm(shuffled_subfam_CpG)

# Log odds ratio (LOR) enrichment for each shuffled subfamily x state x sample combination
# With pseudocount of 1e-20 to avoid undefined values
# Using the updated Length ijk/Length ik for shuffled subfamilies
shuffled_enrichment = lapply(shuffled_enrichment,function(x) transform(x,Enrichment.Shuffled = ifelse(State %in% meth_states,log2(1e-20+((x$Length_ijk.x/x$Length_jk)/(x$CpGs/x$Length_k))),
                                                                                                      ifelse(metadata[match(Sample,metadata$Sample),]$chrY == "Yes",
                                                                                                             log2(1e-20+((x$Length_ijk.x/x$Length_jk)/(x$Total_length/x$Length_k))),
                                                                                                             log2(1e-20+((x$Length_ijk.x/x$Length_jk)/(x$Total_length_noY/x$Length_k)))))))

# Filter LOR enrichments to only those where the subfamily has >10 members in the state and >30 overall
# For the methylation states, considers only TEs that overlap at least one CpG
shuffled_enrichment_filter = lapply(shuffled_enrichment,function(x) x[which(x$Members.x > THRESHOLD_IJK_MEMBER & x$Count > THRESHOLD_IK_MEMBER),])

# Number/proportion of samples with LOR enrichment > 1.5 per subfamily x state, including only those that pass subfamily member thresholds
## Including subfamily x state combinations without any enrichment LOR > 1.5
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

# Number of subfamilies with at least one enrichment > 1.5 and total number of enrichments LOR > 1.5 per state, by iteration
shuffled_enrichment_counts_combine = ddply(shuffled_enrichment_counts,.(Iteration,State),summarise,Subfamilies=sum(V1 > 0),Enrichments=sum(V1))
