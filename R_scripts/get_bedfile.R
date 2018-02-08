subfamily = "LFSINE_Vert"
state = "7_Enh"
metric = "chromHMM"
enrichment=0
  
# Get members ever in state
subfamily_in_state = get_subfamily_in_state(subfamily,state,metric)
  
# Samples where subfamily is enriched in state
samples = get_enriched_samples(subfamily,state,enrichment)
  
# Plot figure
investigate_candidate_indv(subfamily,state,metric,enrichment,print_fig=FALSE)

# Samples of interest
samples_filter = as.vector(metadata[which(metadata$Sample %in% samples & metadata$Anatomy == "BRAIN"),]$Sample)
length(samples_filter)
  
# Subfamily members in the state in samples of interest, bed format with number of samples
subfamily_ubiq = aggregate(data=subfamily_in_state,Sample~chromosome+start+stop+subfamily+strand,length)
colnames(subfamily_ubiq)[6] = "Total_samples"
subfamily_bedfile = aggregate(data=subfamily_in_state[which(subfamily_in_state$Sample %in% samples_filter),],
                              Sample~chromosome+start+stop+subfamily+strand,length)
subfamily_bedfile = merge(subfamily_bedfile,subfamily_ubiq,by=c("chromosome","start","stop","subfamily","strand"))
subfamily_bedfile[order(subfamily_bedfile$Sample,subfamily_bedfile$chromosome,subfamily_bedfile$start),]
  
write.table(subfamily_bedfile[which(subfamily_bedfile$Sample > 1 & subfamily_bedfile$Total_samples < 70),c(1:4,6,5)],
            file=paste("enrichment/",subfamily,"_",state,"_enrich",enrichment,"_Brain.bed",sep=""),quote=FALSE,row.names = FALSE,col.names=FALSE,sep='\t')
