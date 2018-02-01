subfamily = "AluSx"
state = "6_EnhG"
metric = "chromHMM"
enrichment=0
  
# Get members ever in state
subfamily_in_state = get_subfamily_in_state(subfamily,state,metric)
  
# Samples where subfamily is enriched in state
samples = get_enriched_samples(subfamily,state,enrichment)
  
# Samples of interest
samples_filter = as.vector(metadata[which(metadata$Sample %in% samples & metadata$Group == "Blood & T-cell"),]$Sample)
  
# Plot figure
#investigate_candidate_indv(subfamily,state,metric,enrichment,print_fig=FALSE)

# Subfamily members in the state in samples of interest, bed format with number of samples
subfamily_bedfile = aggregate(data=subfamily_in_state[which(subfamily_in_state$Sample %in% samples_filter),],
                              Sample~chromosome+start+stop+subfamily+strand,length)[,c(1:4,6,5)]
subfamily_bedfile[order(subfamily_bedfile$Sample,subfamily_bedfile$chromosome,subfamily_bedfile$start),]
  
write.table(subfamily_bedfile[which(subfamily_bedfile$Sample > 0),],
            file=paste("enrichment/",subfamily,"_",state,"_enriched_select.bed",sep=""),quote=FALSE,row.names = FALSE,col.names=FALSE,sep='\t')